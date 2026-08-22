# =============================================================================
# Nocturnal VMS Electrodermal Activity (EDA) Analysis Pipeline
# -----------------------------------------------------------------------------
# Reorganized for the public repo. Roughly follows how the analysis was run: 
#import data, clean it up, figure out event detection, validate it,
# pull metrics, then the actual clustering (fit -> drop the bad cluster ->
# refit), and last, the stuff that didn't make the final cut but is worth
# keeping around (rejected approaches, etc).
#
# Heads up: no raw or de-identified data lives in this repo. The file paths
# below are the original local ones (kept for transparency about the folder
# structure) - they won't resolve if you try to run this as-is.
# =============================================================================


################################################################################
# DATA IMPORT — raw EDA to 1-min summaries, attach REDCap start times
################################################################################

# Step 1 FIXED: Create 1-minute summaries from 32Hz EDA data
# ------------------------------------------------------------

library(tidyverse)

# Set up paths
base_path <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Michael Busa's files - Data"
output_path <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA Data Summaries"

# Create output directory if it doesn't exist
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
}

# Define participant IDs and visits
participants <- paste0("S", 37:295)
visits <- c("V2", "V3", "V4", "V5")

# Fixed function to create 1-minute summaries
create_minute_summary <- function(participant, visit, base_path, output_path) {
  
  # Construct path to EDA file - look specifically for eda.csv
  eda_file <- file.path(base_path, participant, visit, "movisens EDA", "csv", "eda.csv")
  
  # Check if file exists
  if (!file.exists(eda_file)) {
    return(NULL)  # Skip if file doesn't exist
  }
  
  # Read the EDA file with semicolon delimiter
  eda_data <- tryCatch({
    read.csv(eda_file, sep = ";", header = FALSE)
  }, error = function(e) {
    warning(paste("Error reading file for", participant, visit, ":", e$message))
    return(NULL)
  })
  
  if (is.null(eda_data)) {
    return(NULL)
  }
  
  # Take the first column (EDA values)
  eda_values <- as.numeric(as.character(eda_data[, 1]))
  
  # Check if we got valid data
  if (all(is.na(eda_values))) {
    warning(paste("All NA values for", participant, visit))
    return(NULL)
  }
  
  # Remove any leading/trailing NA values
  valid_indices <- which(!is.na(eda_values))
  if (length(valid_indices) == 0) {
    warning(paste("No valid data for", participant, visit))
    return(NULL)
  }
  
  first_valid <- min(valid_indices)
  last_valid <- max(valid_indices)
  eda_values <- eda_values[first_valid:last_valid]
  
  # Create 1-minute summaries (32Hz * 60 seconds = 1920 samples per minute)
  samples_per_minute <- 1920
  n_minutes <- floor(length(eda_values) / samples_per_minute)
  
  if (n_minutes == 0) {
    warning(paste("Not enough data for", participant, visit, "- only", length(eda_values), "samples"))
    return(NULL)
  }
  
  # Calculate minute-level means
  minute_summaries <- data.frame(
    minute = 1:n_minutes,
    EDA_mean = sapply(1:n_minutes, function(i) {
      start_idx <- (i - 1) * samples_per_minute + 1
      end_idx <- i * samples_per_minute
      mean(eda_values[start_idx:end_idx], na.rm = TRUE)
    })
  )
  
  # Save summary file
  output_file <- file.path(output_path, paste0(participant, visit, "EDAsummary.csv"))
  write.csv(minute_summaries, output_file, row.names = FALSE)
  
  return(paste("Processed:", participant, visit, "-", n_minutes, "minutes"))
}

# Process all participants and visits
results <- list()
errors <- list()

for (participant in participants) {
  for (visit in visits) {
    result <- create_minute_summary(participant, visit, base_path, output_path)
    if (!is.null(result)) {
      results <- c(results, result)
      print(result)
    }
  }
}

print("=====================================")
print(paste("Total files processed:", length(results)))
print("=====================================")
# Step 2: Add timestamps based on start times from REDCap
# --------------------------------------------------------

library(tidyverse)
library(lubridate)

# Set paths
summary_path <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA Data Summaries"
start_times_file <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/REDCap EDA Monitor Start Times.csv"

# Read start times
start_times <- read.csv(start_times_file)

# Check the structure of start times file
print("Start times file structure:")
print(head(start_times))
print(colnames(start_times))

# Function to add timestamps to summary file
add_timestamps <- function(summary_file, start_times) {
  
  # Extract participant and visit from filename
  # Example: S37V2EDAsummary.csv
  filename <- basename(summary_file)
  participant_num <- as.numeric(str_extract(filename, "(?<=S)\\d+"))
  visit_num <- str_extract(filename, "V\\d+")
  
  # Map visit to REDCap event name
  visit_mapping <- c(
    "V2" = "v2_arm_1",
    "V3" = "v3_arm_1",
    "V4" = "v4_arm_1",
    "V5" = "v5_for_resched_arm_1"
  )
  
  redcap_event <- visit_mapping[visit_num]
  
  # Find corresponding start time
  start_row <- start_times %>%
    filter(study_id == participant_num, 
           redcap_event_name == redcap_event)
  
  if (nrow(start_row) == 0) {
    warning(paste("No start time found for", filename))
    return(NULL)
  }
  
  # Parse start time (format: 4/23/2021 15:30)
  start_time <- mdy_hm(start_row$moveda_start[1])
  
  if (is.na(start_time)) {
    warning(paste("Invalid start time for", filename, "- start time:", start_row$moveda_start[1]))
    return(NULL)
  }
  
  # Read summary data
  summary_data <- read.csv(summary_file)
  
  # Add timestamps (1 minute intervals)
  summary_data$timestamp <- start_time + minutes(0:(nrow(summary_data) - 1))
  
  # Reorder columns to put timestamp first
  summary_data <- summary_data %>%
    select(timestamp, minute, EDA_mean)
  
  # Save updated file
  write.csv(summary_data, summary_file, row.names = FALSE)
  
  return(paste("Added timestamps to:", filename, "- Start:", start_time))
}

# Get all summary files
summary_files <- list.files(summary_path, pattern = "EDAsummary\\.csv$", full.names = TRUE)

print(paste("Found", length(summary_files), "summary files to process"))

# Add timestamps to all files
timestamp_results <- list()
for (file in summary_files) {
  result <- add_timestamps(file, start_times = start_times)
  if (!is.null(result)) {
    timestamp_results <- c(timestamp_results, result)
    print(result)
  }
}

print("=====================================")
print(paste("Timestamps added to", length(timestamp_results), "files"))
print("=====================================")


################################################################################
# PREPROCESSING — sleep trimming, exclusions, rolling CV, unit conversion
################################################################################


# trim to nighttime only, using GGIR sleep windows
### working on trimmed datasets for nighttime only ###
# Function to extract sleep times and trim EDA
process_participant_visit <- function(participant, visit) {
  
  # Construct file paths
  ggir_file <- paste0("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/GGIR_Output/",
                      participant, visit, "/output_Actigraphy/results/part4_nightsummary_sleep_cleaned.csv")
  
  eda_file <- paste0("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA Data Summaries/",
                     participant, visit, "EDAsummary.csv")
  
  output_file <- paste0("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA Data Summaries Trimmed/",
                        participant, visit, "EDAsummary_trimmed.csv")
  
  # Check if files exist
  if(!file.exists(ggir_file)) {
    cat("GGIR file not found for", participant, visit, "\n")
    return(NULL)
  }
  
  if(!file.exists(eda_file)) {
    cat("EDA file not found for", participant, visit, "\n")
    return(NULL)
  }
  
  # Read GGIR sleep data
  sleep_data <- read.csv(ggir_file, stringsAsFactors = FALSE)
  
  # Read EDA data
  eda_data <- read.csv(eda_file, stringsAsFactors = FALSE)
  
  # Check if timestamp column exists
  if(!"timestamp" %in% names(eda_data)) {
    cat("No timestamp column in EDA file for", participant, visit, "\n")
    return(NULL)
  }
  
  # Debug: Check what the timestamp looks like
  cat("First few timestamps:\n")
  print(head(eda_data$timestamp))
  
  # Parse EDA timestamps - with error handling
  eda_data$timestamp <- tryCatch({
    ymd_hms(eda_data$timestamp, quiet = TRUE)
  }, error = function(e) {
    # If ymd_hms fails, try other formats
    parse_date_time(eda_data$timestamp, orders = c("ymd HMS", "mdy HMS", "dmy HMS"))
  })
  
  # Remove rows with NA timestamps
  eda_data <- eda_data %>% filter(!is.na(timestamp))
  
  if(nrow(eda_data) == 0) {
    cat("No valid timestamps in EDA file for", participant, visit, "\n")
    return(NULL)
  }
  
  # Initialize trimmed data
  trimmed_data <- data.frame()
  
  # Process each night
  for(night_num in unique(sleep_data$night)) {
    night_info <- sleep_data %>% filter(night == night_num)
    
    if(nrow(night_info) == 0) next
    
    # Get calendar date and sleep times
    cal_date <- ymd(night_info$calendar_date[1])
    sleeponset_decimal <- as.numeric(night_info$sleeponset[1])
    wakeup_decimal <- as.numeric(night_info$wakeup[1])
    
    # Skip if missing data
    if(is.na(sleeponset_decimal) || is.na(wakeup_decimal)) {
      cat("Missing sleep data for", participant, visit, "night", night_num, "\n")
      next
    }
    
    # Convert decimal hours to datetime
    onset_hours <- floor(sleeponset_decimal)
    onset_minutes <- round((sleeponset_decimal - onset_hours) * 60)
    sleep_onset <- cal_date + hours(onset_hours) + minutes(onset_minutes)
    
    wake_hours <- floor(wakeup_decimal)
    wake_minutes <- round((wakeup_decimal - wake_hours) * 60)
    wake_time <- cal_date + hours(wake_hours) + minutes(wake_minutes)
    
    # If wake time is before sleep onset, it's the next day
    if(wake_time < sleep_onset) {
      wake_time <- wake_time + days(1)
    }
    
    cat("Night", night_num, "- Sleep:", sleep_onset, "Wake:", wake_time, "\n")
    
    # Filter EDA data for this sleep period
    night_eda <- eda_data %>%
      filter(timestamp >= sleep_onset & timestamp <= wake_time) %>%
      mutate(night = night_num,
             sleep_onset = as.character(sleep_onset),
             wake_time = as.character(wake_time))
    
    trimmed_data <- bind_rows(trimmed_data, night_eda)
  }
  
  # Save trimmed data
  if(nrow(trimmed_data) > 0) {
    write.csv(trimmed_data, output_file, row.names = FALSE)
    cat("Processed", participant, visit, "- Rows:", nrow(trimmed_data), "\n")
  } else {
    cat("No data extracted for", participant, visit, "\n")
  }
  
  return(trimmed_data)
}

# Test on one participant first
process_participant_visit("S70", "V2")





















##################
# Get all participant/visit combinations from the GGIR output folder
ggir_folders <- list.dirs("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/GGIR_Output", 
                          full.names = FALSE, recursive = FALSE)

# Process each one
results_summary <- data.frame()

for(folder in ggir_folders) {
  # Extract participant and visit from folder name (e.g., "S51V2")
  # Handle different ID lengths (S37 vs S101)
  if(grepl("V", folder)) {
    # Find where V starts
    v_position <- regexpr("V[2-5]", folder)
    participant <- substr(folder, 1, v_position - 1)  # Everything before V
    visit <- substr(folder, v_position, v_position + 1)  # V2, V3, V4, or V5
    
    result <- tryCatch({
      trimmed <- process_participant_visit(participant, visit)
      if(!is.null(trimmed)) {
        data.frame(
          participant = participant,
          visit = visit,
          rows = nrow(trimmed),
          nights = length(unique(trimmed$night)),
          status = "Success"
        )
      } else {
        data.frame(
          participant = participant,
          visit = visit,
          rows = 0,
          nights = 0,
          status = "Failed"
        )
      }
    }, error = function(e) {
      cat("ERROR for", participant, visit, ":", e$message, "\n")
      data.frame(
        participant = participant,
        visit = visit,
        rows = 0,
        nights = 0,
        status = paste0("Error: ", e$message)
      )
    })
    
    results_summary <- bind_rows(results_summary, result)
  }
}

print("All files processed!")
print(results_summary)

# Save the summary
write.csv(results_summary, 
          "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_trimming_summary.csv",
          row.names = FALSE)

# apply participant/visit/night exclusions, then compute 10-min rolling CV
# (full exclusion list -> Supplemental Table 1)
#### EXCLUDING PARTICIPANTS AND ROLLING WINDOW ANALYSIS #####

library(tidyverse)
library(lubridate)
library(zoo)

# ===== STEP 1: DEFINE EXCLUSIONS =====

# Complete visit exclusions (remove entire file)
full_exclusions <- c(
  "S57V2",
  "S75V2",
  "S79V2", "S79V3", "S79V4",
  "S91V2", "S91V3", "S91V4",
  "S101V4",
  "S103V3", "S103V4",
  "S131V2", "S131V3", "S131V4",
  "S154V4",
  "S176V2", "S176V3", "S176V4", "S176V5",
  "S198V3",
  "S242V2",
  "S266V2", "S266V3", "S266V4"
)

# Single night exclusions - format: list of c("participantVisit", night_to_EXCLUDE)
# The night listed is the one to DELETE, keeping the other
night_exclusions <- list(
  c("S69V3",  2),   # exclude night 2, keep night 1
  c("S67V4",  2),   # exclude night 2, keep night 1
  c("S86V3",  1),   # exclude night 1, keep night 2
  c("S86V4",  2),   # exclude night 2, keep night 1
  c("S88V4",  2),   # exclude night 2, keep night 1
  c("S101V2", 2),   # exclude night 2, keep night 1
  c("S101V3", 2),   # exclude night 2, keep night 1
  c("S104V3", 1),   # exclude night 1, keep night 2
  c("S104V4", 2),   # exclude night 2, keep night 1
  c("S146V2", 1),   # exclude night 1, keep night 2
  c("S154V3", 2),   # exclude night 2, keep night 1
  c("S226V3", 2),   # exclude night 2, keep night 1
  c("S226V4", 1),   # exclude night 1, keep night 2
  c("S250V4", 2)    # exclude night 2, keep night 1
)

# ===== STEP 2: LOAD AND FILTER ALL TRIMMED FILES =====

trimmed_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA Data Summaries Trimmed"

trimmed_files <- list.files(trimmed_folder, pattern = "EDAsummary_trimmed.csv", full.names = TRUE)

cat("Total trimmed files before exclusions:", length(trimmed_files), "\n\n")

# Load all files into one dataframe, applying exclusions as we go
all_eda <- data.frame()

for(file_path in trimmed_files) {
  
  # Extract participant and visit from filename
  file_name <- basename(file_path)
  file_base <- gsub("EDAsummary_trimmed.csv", "", file_name)
  v_position <- regexpr("V[2-5]", file_base)
  participant <- substr(file_base, 1, v_position - 1)
  visit <- substr(file_base, v_position, v_position + 1)
  pv_label <- paste0(participant, visit)  # e.g., "S70V2"
  
  # Skip full exclusions
  if(pv_label %in% full_exclusions) {
    cat("Excluding full visit:", pv_label, "\n")
    next
  }
  
  # Read the file
  eda_data <- read.csv(file_path, stringsAsFactors = FALSE)
  eda_data$timestamp <- ymd_hms(eda_data$timestamp)
  eda_data$participant <- participant
  eda_data$visit <- visit
  eda_data$pv_label <- pv_label
  
  # Apply single-night exclusions
  for(excl in night_exclusions) {
    if(pv_label == excl[1]) {
      night_to_remove <- as.numeric(excl[2])
      rows_before <- nrow(eda_data)
      eda_data <- eda_data %>% filter(night != night_to_remove)
      cat("Excluding", pv_label, "night", night_to_remove, 
          "- removed", rows_before - nrow(eda_data), "rows\n")
    }
  }
  
  # Add to master dataframe
  all_eda <- bind_rows(all_eda, eda_data)
}

cat("\n=== Loading complete ===\n")

# ===== STEP 3: DESCRIPTIVE STATISTICS =====

# Participants remaining
participants_remaining <- length(unique(all_eda$participant))

# Visits remaining
visits_remaining <- all_eda %>%
  distinct(participant, visit) %>%
  nrow()

# Nights remaining (each participant/visit/night combo)
nights_remaining <- all_eda %>%
  distinct(participant, visit, night) %>%
  nrow()

# Night duration (each night is a unique participant/visit/night combination)
night_durations <- all_eda %>%
  group_by(participant, visit, night) %>%
  summarise(
    start_time = min(timestamp),
    end_time = max(timestamp),
    duration_minutes = as.numeric(difftime(max(timestamp), min(timestamp), units = "mins")),
    .groups = "drop"
  )

avg_duration <- mean(night_durations$duration_minutes)
min_duration <- min(night_durations$duration_minutes)
max_duration <- max(night_durations$duration_minutes)
sd_duration  <- sd(night_durations$duration_minutes)

cat("\n========================================\n")
cat("      DESCRIPTIVE STATISTICS\n")
cat("========================================\n")
cat("Participants remaining:    ", participants_remaining, "\n")
cat("Visits remaining:          ", visits_remaining, "\n")
cat("Nights remaining:          ", nights_remaining, "\n")
cat("----------------------------------------\n")
cat("Avg night duration:        ", round(avg_duration, 1), "minutes (", round(avg_duration/60, 2), "hours)\n")
cat("SD night duration:         ", round(sd_duration, 1), "minutes\n")
cat("Min night duration:        ", round(min_duration, 1), "minutes (", round(min_duration/60, 2), "hours)\n")
cat("Max night duration:        ", round(max_duration, 1), "minutes (", round(max_duration/60, 2), "hours)\n")
cat("========================================\n")

# Breakdown by participant - how many visits and nights each has
participant_summary <- all_eda %>%
  group_by(participant) %>%
  summarise(
    n_visits = n_distinct(visit),
    n_nights = n_distinct(paste(visit, night)),
    .groups = "drop"
  )

cat("\nPer-participant breakdown:\n")
print(participant_summary)

# ===== STEP 4: ROLLING WINDOWS (10 MINUTE) WITH COEFFICIENT OF VARIATION =====

# Coefficient of variation = (SD / Mean) * 100
# This expresses variability RELATIVE to the mean
# High CV = lots of variability relative to average = potential event
# Low CV = stable EDA = quiet period

calculate_rolling_cv <- function(eda_data, window_size_minutes = 10) {
  
  eda_data <- eda_data %>% arrange(timestamp)
  
  eda_data <- eda_data %>%
    mutate(
      rolling_mean = rollmean(EDA_mean, k = window_size_minutes, fill = NA, align = "right"),
      rolling_sd   = rollapply(EDA_mean, width = window_size_minutes, FUN = sd, fill = NA, align = "right"),
      rolling_max  = rollapply(EDA_mean, width = window_size_minutes, FUN = max, fill = NA, align = "right"),
      rolling_min  = rollapply(EDA_mean, width = window_size_minutes, FUN = min, fill = NA, align = "right"),
      
      # Coefficient of variation (as percentage)
      rolling_cv   = (rolling_sd / rolling_mean) * 100,
      
      window_size  = window_size_minutes
    )
  
  # Remove incomplete windows
  eda_data <- eda_data %>% filter(!is.na(rolling_mean))
  
  return(eda_data)
}

# Apply rolling CV to entire dataset, processing each participant/visit/night separately
# (important to not let windows bleed across nights or participants)
all_eda_cv <- all_eda %>%
  group_by(participant, visit, night) %>%
  group_modify(~ calculate_rolling_cv(.x, window_size_minutes = 10)) %>%
  ungroup()

cat("\nRolling CV calculation complete!\n")
cat("Total rows with CV data:", nrow(all_eda_cv), "\n")

# Save the processed dataset
output_path <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min.csv"
write.csv(all_eda_cv, output_path, row.names = FALSE)
cat("Saved to:", output_path, "\n")

# build the final trimmed rolling-CV dataset (cutting anomalous/cutoff events)
#### creating a new rolling dataset that is trimmed of end cut off events and anomalies #####
library(tidyverse)
library(lubridate)

# Load the rolling CV dataset
rolling_cv <- read.csv("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min.csv",
                       stringsAsFactors = FALSE)

# Parse timestamp
rolling_cv$timestamp <- ymd_hms(rolling_cv$timestamp)

cat("Original dataset:\n")
cat("  Rows:", nrow(rolling_cv), "\n")
cat("  Participants:", n_distinct(rolling_cv$participant), "\n")
cat("  Visits:", n_distinct(paste0(rolling_cv$participant, rolling_cv$visit)), "\n\n")

# ===== STEP 1: FULL VISIT EXCLUSIONS =====
full_exclusions <- c(
  "S57V2",
  "S75V2",
  "S79V2", "S79V3", "S79V4",
  "S91V2", "S91V3", "S91V4",
  "S101V4",
  "S103V3", "S103V4",
  "S131V2", "S131V3", "S131V4",
  "S154V4",
  "S176V2", "S176V3", "S176V4", "S176V5",
  "S198V3",
  "S242V2",
  "S266V2", "S266V3", "S266V4",
  "S86V4",
  "S282V4"
)

rolling_cv <- rolling_cv %>%
  filter(!pv_label %in% full_exclusions)

cat("After full visit exclusions:", nrow(rolling_cv), "rows\n\n")

# ===== STEP 2: SINGLE NIGHT EXCLUSIONS =====
night_exclusions <- tribble(
  ~pv_label, ~night,
  "S69V3",   2,
  "S67V4",   2,
  "S86V3",   1,
  "S88V4",   2,
  "S101V2",  2,
  "S101V3",  2,
  "S104V3",  1,
  "S104V4",  2,
  "S146V2",  1,
  "S154V3",  2,
  "S226V3",  2,
  "S226V4",  1,
  "S250V4",  2,
  "S269V3",  2,
  "S244V3",  1,
  "S146V4",  1,
  "S143V2",  1
)

for(i in 1:nrow(night_exclusions)) {
  pv <- night_exclusions$pv_label[i]
  n  <- night_exclusions$night[i]
  
  before <- nrow(rolling_cv)
  rolling_cv <- rolling_cv %>%
    filter(!(pv_label == pv & night == n))
  after <- nrow(rolling_cv)
  
  if(before > after) {
    cat("Excluded", pv, "night", n, "- removed", before - after, "rows\n")
  }
}

cat("\nAfter night exclusions:", nrow(rolling_cv), "rows\n\n")

# ===== STEP 3: TIME-BASED TRIMMING =====
# Helper function to create time on the correct date
create_time_filter <- function(data, pv, night_num, before_time = NULL, after_time = NULL) {
  
  # Get data for this specific night
  night_data <- data %>% filter(pv_label == pv, night == night_num)
  
  if(nrow(night_data) == 0) return(data)
  
  # Get the date of this night (use the calendar date from sleep onset)
  night_date <- as.Date(min(night_data$timestamp))
  
  rows_before <- nrow(data)
  
  if(!is.null(before_time)) {
    # Keep data BEFORE this time (i.e., remove after)
    cutoff <- ymd_hm(paste(night_date, before_time))
    # If cutoff is before the start of the night, it's actually the next day
    if(cutoff < min(night_data$timestamp)) {
      cutoff <- cutoff + days(1)
    }
    data <- data %>%
      filter(!(pv_label == pv & night == night_num & timestamp >= cutoff))
  }
  
  if(!is.null(after_time)) {
    # Keep data AFTER this time (i.e., remove before)
    cutoff <- ymd_hm(paste(night_date, after_time))
    # If cutoff is way after, might be next day
    if(cutoff < min(night_data$timestamp)) {
      cutoff <- cutoff + days(1)
    }
    data <- data %>%
      filter(!(pv_label == pv & night == night_num & timestamp <= cutoff))
  }
  
  rows_after <- nrow(data)
  if(rows_before > rows_after) {
    cat("Trimmed", pv, "night", night_num, "- removed", rows_before - rows_after, "rows\n")
  }
  
  return(data)
}

# Apply all time trims
rolling_cv <- create_time_filter(rolling_cv, "S64V3", 1, before_time = "06:00")
rolling_cv <- create_time_filter(rolling_cv, "S94V3", 1, before_time = "07:30")
rolling_cv <- create_time_filter(rolling_cv, "S94V3", 2, before_time = "05:00")
rolling_cv <- create_time_filter(rolling_cv, "S88V3", 1, before_time = "06:00")
rolling_cv <- create_time_filter(rolling_cv, "S88V2", 1, before_time = "07:00")
rolling_cv <- create_time_filter(rolling_cv, "S86V3", 2, before_time = "06:00")
rolling_cv <- create_time_filter(rolling_cv, "S75V4", 2, before_time = "06:30")
rolling_cv <- create_time_filter(rolling_cv, "S70V3", 1, before_time = "07:30")
rolling_cv <- create_time_filter(rolling_cv, "S70V3", 2, before_time = "08:00")
rolling_cv <- create_time_filter(rolling_cv, "S70V2", 1, before_time = "06:15")
rolling_cv <- create_time_filter(rolling_cv, "S69V2", 2, before_time = "05:50")
rolling_cv <- create_time_filter(rolling_cv, "S295V5", 1, before_time = "07:30")
rolling_cv <- create_time_filter(rolling_cv, "S293V3", 1, before_time = "06:00")
rolling_cv <- create_time_filter(rolling_cv, "S292V2", 1, before_time = "08:25")
rolling_cv <- create_time_filter(rolling_cv, "S288V4", 1, before_time = "05:00")
rolling_cv <- create_time_filter(rolling_cv, "S288V3", 2, after_time = "22:30")
rolling_cv <- create_time_filter(rolling_cv, "S274V3", 1, before_time = "06:00")
rolling_cv <- create_time_filter(rolling_cv, "S270V2", 2, after_time = "12:30")
rolling_cv <- create_time_filter(rolling_cv, "S260V2", 1, before_time = "03:30")
rolling_cv <- create_time_filter(rolling_cv, "S260V2", 2, before_time = "04:00")
rolling_cv <- create_time_filter(rolling_cv, "S250V2", 1, after_time = "22:30")
rolling_cv <- create_time_filter(rolling_cv, "S250V2", 1, before_time = "06:00")
rolling_cv <- create_time_filter(rolling_cv, "S238V3", 1, before_time = "06:10")
rolling_cv <- create_time_filter(rolling_cv, "S226V2", 1, before_time = "10:00")
rolling_cv <- create_time_filter(rolling_cv, "S226V2", 2, before_time = "07:50")
rolling_cv <- create_time_filter(rolling_cv, "S225V3", 1, after_time = "00:00")
rolling_cv <- create_time_filter(rolling_cv, "S225V3", 2, before_time = "06:00")
rolling_cv <- create_time_filter(rolling_cv, "S213V4", 1, after_time = "23:15")
rolling_cv <- create_time_filter(rolling_cv, "S213V4", 1, before_time = "07:00")
rolling_cv <- create_time_filter(rolling_cv, "S213V3", 1, after_time = "23:00")
rolling_cv <- create_time_filter(rolling_cv, "S213V2", 1, after_time = "23:00")
rolling_cv <- create_time_filter(rolling_cv, "S213V2", 2, before_time = "07:15")
rolling_cv <- create_time_filter(rolling_cv, "S198V2", 1, before_time = "05:30")
rolling_cv <- create_time_filter(rolling_cv, "S196V3", 2, before_time = "07:30")
rolling_cv <- create_time_filter(rolling_cv, "S196V2", 2, before_time = "07:00")
rolling_cv <- create_time_filter(rolling_cv, "S193V3", 2, before_time = "04:00")
rolling_cv <- create_time_filter(rolling_cv, "S193V2", 1, after_time = "23:00")
rolling_cv <- create_time_filter(rolling_cv, "S193V2", 2, before_time = "05:00")
rolling_cv <- create_time_filter(rolling_cv, "S166V3", 1, before_time = "07:00")
rolling_cv <- create_time_filter(rolling_cv, "S154V2", 1, before_time = "10:50")
rolling_cv <- create_time_filter(rolling_cv, "S148V2", 1, after_time = "22:15")
rolling_cv <- create_time_filter(rolling_cv, "S148V2", 1, before_time = "06:00")
rolling_cv <- create_time_filter(rolling_cv, "S148V2", 2, before_time = "06:00")
rolling_cv <- create_time_filter(rolling_cv, "S146V4", 2, before_time = "05:15")
rolling_cv <- create_time_filter(rolling_cv, "S143V4", 2, before_time = "05:00")
rolling_cv <- create_time_filter(rolling_cv, "S143V3", 1, before_time = "06:30")
rolling_cv <- create_time_filter(rolling_cv, "S143V2", 2, before_time = "02:10")
rolling_cv <- create_time_filter(rolling_cv, "S122V4", 2, before_time = "05:00")
rolling_cv <- create_time_filter(rolling_cv, "S122V3", 2, before_time = "05:00")
rolling_cv <- create_time_filter(rolling_cv, "S122V2", 2, before_time = "04:30")
rolling_cv <- create_time_filter(rolling_cv, "S113V2", 2, after_time = "22:30")
rolling_cv <- create_time_filter(rolling_cv, "S104V4", 1, before_time = "05:00")
rolling_cv <- create_time_filter(rolling_cv, "S104V2", 2, before_time = "04:45")

cat("\nAfter time trimming:", nrow(rolling_cv), "rows\n\n")

# ===== SAVE TRIMMED DATASET =====
output_path <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min_trimmed.csv"
write.csv(rolling_cv, output_path, row.names = FALSE)
cat("Saved trimmed dataset to:", output_path, "\n\n")

# ===== DESCRIPTIVE STATISTICS =====
participants <- n_distinct(rolling_cv$participant)
visits <- n_distinct(paste0(rolling_cv$participant, rolling_cv$visit))
nights <- rolling_cv %>%
  distinct(participant, visit, night) %>%
  nrow()

# Calculate night durations
night_durations <- rolling_cv %>%
  group_by(participant, visit, night) %>%
  summarise(
    duration_minutes = as.numeric(difftime(max(timestamp), min(timestamp), units = "mins")),
    .groups = "drop"
  )

cat("========================================\n")
cat("      TRIMMED DATASET STATISTICS\n")
cat("========================================\n")
cat("Participants:          ", participants, "\n")
cat("Visits:                ", visits, "\n")
cat("Nights:                ", nights, "\n")
cat("----------------------------------------\n")
cat("Avg night duration:    ", round(mean(night_durations$duration_minutes), 1), " min (", 
    round(mean(night_durations$duration_minutes)/60, 2), " hrs)\n")
cat("SD night duration:     ", round(sd(night_durations$duration_minutes), 1), " min\n")
cat("Min night duration:    ", round(min(night_durations$duration_minutes), 1), " min (", 
    round(min(night_durations$duration_minutes)/60, 2), " hrs)\n")
cat("Max night duration:    ", round(max(night_durations$duration_minutes), 1), " min (", 
    round(max(night_durations$duration_minutes)/60, 2), " hrs)\n")
cat("========================================\n")

# convert EDA from microsiemens to micromhos
##########changing eda to uMhos and re-running visualizations############
library(tidyverse)
library(lubridate)

# Load trimmed dataset
rolling_cv_trimmed <- read.csv("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min_trimmed.csv",
                               stringsAsFactors = FALSE)
rolling_cv_trimmed$timestamp <- ymd_hms(rolling_cv_trimmed$timestamp)

# ===== ADD MHOS COLUMN =====
# Conversion: microsiemens to micromhos
rolling_cv_trimmed <- rolling_cv_trimmed %>%
  mutate(mhos = EDA_mean * 0.00305185)

cat("Added mhos column to dataset\n")
cat("Sample mhos values:\n")
print(head(rolling_cv_trimmed %>% select(EDA_mean, mhos)))

# Save updated dataset
output_path <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min_trimmed.csv"
write.csv(rolling_cv_trimmed, output_path, row.names = FALSE)
cat("\nUpdated dataset saved with mhos column\n\n")

# Create output folder
output_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA CV Events 10-10-75-5min MHOS"
dir.create(output_folder, showWarnings = FALSE)

# ===== DUAL CRITERIA EVENT DETECTION (using mhos) =====
identify_events_dual_5min <- function(cv_data, cv_start = 10, cv_end = 10, eda_percentile = 75) {
  cv_data <- cv_data %>% arrange(timestamp)
  # Calculate threshold based on MHOS not EDA_mean
  mhos_threshold <- quantile(cv_data$mhos, eda_percentile/100, na.rm = TRUE)
  cv_data$in_event <- FALSE
  cv_data$mhos_threshold <- mhos_threshold
  in_event <- FALSE
  both_conditions_start <- NULL
  
  for(i in 1:nrow(cv_data)) {
    if(!in_event) {
      if(!is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] > cv_start) {
        in_event <- TRUE
        both_conditions_start <- NULL
      }
    } else {
      cv_low <- !is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] < cv_end
      mhos_low <- !is.na(cv_data$mhos[i]) && cv_data$mhos[i] < mhos_threshold
      both_conditions_met <- cv_low && mhos_low
      
      if(both_conditions_met && is.null(both_conditions_start)) {
        both_conditions_start <- i
      }
      if(!both_conditions_met) {
        both_conditions_start <- NULL
      }
      if(both_conditions_met && !is.null(both_conditions_start)) {
        minutes_both_low <- i - both_conditions_start
        if(minutes_both_low >= 5) {
          in_event <- FALSE
          both_conditions_start <- NULL
        }
      }
    }
    cv_data$in_event[i] <- in_event
  }
  return(cv_data)
}

# ===== GET EVENT RECTANGLES =====
get_event_rects <- function(night_data) {
  night_data <- night_data %>% arrange(timestamp)
  rects <- list()
  event_start <- NULL
  
  for(i in 1:nrow(night_data)) {
    if(night_data$in_event[i] && (i == 1 || !night_data$in_event[i-1])) {
      event_start <- night_data$timestamp[i]
    }
    if(!is.null(event_start) && !night_data$in_event[i] && i > 1 && night_data$in_event[i-1]) {
      rects[[length(rects)+1]] <- data.frame(
        xmin = event_start,
        xmax = night_data$timestamp[i],
        night_label = night_data$night_label[i]
      )
      event_start <- NULL
    }
    if(i == nrow(night_data) && !is.null(event_start) && night_data$in_event[i]) {
      rects[[length(rects)+1]] <- data.frame(
        xmin = event_start,
        xmax = night_data$timestamp[i],
        night_label = night_data$night_label[i]
      )
    }
  }
  if(length(rects) == 0) return(data.frame())
  return(bind_rows(rects))
}

# ===== VISUALIZATION FUNCTION =====
visualize_mhos <- function(pv, data = rolling_cv_trimmed) {
  
  plot_data <- data %>%
    filter(pv_label == pv) %>%
    mutate(night_label = paste("Night", night))
  
  if(nrow(plot_data) == 0) {
    return(NULL)
  }
  
  pid <- plot_data$participant[1]
  visit <- plot_data$visit[1]
  nights <- unique(plot_data$night_label)
  
  # Process each night with 10/10/75
  processed <- map_dfr(nights, function(n) {
    identify_events_dual_5min(
      plot_data %>% filter(night_label == n),
      cv_start = 10,
      cv_end = 10,
      eda_percentile = 75
    )
  })
  
  # Get rectangles
  all_rects <- map_dfr(nights, function(n) {
    get_event_rects(processed %>% filter(night_label == n))
  })
  
  # Event counts per night
  if(nrow(all_rects) > 0) {
    event_counts <- all_rects %>%
      group_by(night_label) %>%
      summarise(n = n(), .groups = "drop")
  } else {
    event_counts <- data.frame(night_label = nights, n = 0)
  }
  
  event_counts <- data.frame(night_label = nights) %>%
    left_join(event_counts, by = "night_label") %>%
    mutate(n = ifelse(is.na(n), 0, n))
  
  subtitle_text <- paste(
    paste0(event_counts$night_label, ": ", event_counts$n, " events"),
    collapse = "   |   "
  )
  
  # Scale MHOS to CV axis
  cv_max <- max(processed$rolling_cv, na.rm = TRUE)
  mhos_min <- min(processed$mhos, na.rm = TRUE)
  mhos_max <- max(processed$mhos, na.rm = TRUE)
  
  processed <- processed %>%
    mutate(
      mhos_scaled = ((mhos - mhos_min) / (mhos_max - mhos_min)) * cv_max,
      thresh_scaled = ((mhos_threshold - mhos_min) / (mhos_max - mhos_min)) * cv_max
    )
  
  # Build plot
  p <- ggplot(processed, aes(x = timestamp)) +
    facet_wrap(~ night_label, scales = "free_x", ncol = 1)
  
  # Shading
  if(nrow(all_rects) > 0) {
    p <- p + geom_rect(
      data = all_rects,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = "#E63946", alpha = 0.2,
      inherit.aes = FALSE
    )
  }
  
  p <- p +
    # MHOS line - darker color (slate blue instead of grey)
    geom_line(aes(y = mhos_scaled), color = "#507eab", linewidth = 0.6, alpha = 0.9) +
    # CV line
    geom_line(aes(y = rolling_cv), color = "black", linewidth = 0.5) +
    # MHOS threshold line
    geom_line(aes(y = thresh_scaled), color = "orange", linewidth = 0.8, alpha = 0.7) +
    # CV threshold
    geom_hline(yintercept = 10, linetype = "dashed", color = "#E63946", linewidth = 0.8) +
    scale_x_datetime(date_breaks = "1 hour", date_labels = "%H:%M") +
    scale_y_continuous(
      name = "CV (%)",
      sec.axis = sec_axis(
        trans = ~ . * (mhos_max - mhos_min) / cv_max + mhos_min,
        name = "EDA (μmhos)"
      )
    ) +
    labs(
      title = paste0(pid, " ", visit, " — Event Detection"),
      subtitle = paste0(subtitle_text, "   |   Start: CV>10%   |   End: CV<10% AND EDA<75th percentile (orange) for 5+ min")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title.y.right = element_text(color = "#507eab", size = 10, face = "bold"),
      axis.text.y.right = element_text(color = "#507eab", size = 9),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "black"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1.5, "lines")
    )
  
  # Save
  output_path <- paste0(output_folder, "/", pid, "_", visit, "_mhos.png")
  ggsave(filename = output_path, plot = p, width = 14, height = 8, dpi = 300)
  
  total_events <- if(nrow(all_rects) > 0) nrow(all_rects) else 0
  cat("Saved:", pid, visit, "-", total_events, "events\n")
  
  return(total_events)
}

# ===== RUN FOR ALL PARTICIPANTS =====
pv_list <- rolling_cv_trimmed %>%
  distinct(pv_label) %>%
  arrange(pv_label) %>%
  pull(pv_label)

cat("Creating MHOS visualizations for", length(pv_list), "visits...\n\n")

for(pv in pv_list) {
  tryCatch({
    visualize_mhos(pv)
  }, error = function(e) {
    cat("ERROR for", pv, ":", e$message, "\n")
  })
}

cat("\n=====================================\n")
cat("All MHOS visualizations complete!\n")
cat("Saved to:", output_folder, "\n")
cat("=====================================\n")


################################################################################
# EVENT DETECTION — tried a few threshold combos before landing on the final one
################################################################################

# (this matches the Methods line about testing a range of onset/offset
# thresholds against participant-marked events - keeping all of them here,
# clearly labeled, since that's part of the actual method.)

# tried: CV on/off = 10%/10%, EDA offset = 75th percentile, no 5-min rule
####### UTILIZING DUAL METRIC ID FOR ALL PARTICIPANTS ####################
library(tidyverse)
library(lubridate)

rolling_cv_trimmed <- read.csv("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min_trimmed.csv",
                               stringsAsFactors = FALSE)
rolling_cv_trimmed$timestamp <- ymd_hms(rolling_cv_trimmed$timestamp)

# Create output folder
output_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA CV Events 10-10-75pct"
dir.create(output_folder, showWarnings = FALSE)

# ===== DUAL CRITERIA EVENT DETECTION =====
identify_events_dual <- function(cv_data, cv_start = 10, cv_end = 10, eda_percentile = 75) {
  cv_data <- cv_data %>% arrange(timestamp)
  eda_threshold <- quantile(cv_data$EDA_mean, eda_percentile/100, na.rm = TRUE)
  cv_data$in_event <- FALSE
  cv_data$eda_threshold <- eda_threshold
  in_event <- FALSE
  
  for(i in 1:nrow(cv_data)) {
    if(!in_event) {
      if(!is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] > cv_start) {
        in_event <- TRUE
      }
    } else {
      cv_low <- !is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] < cv_end
      eda_low <- !is.na(cv_data$EDA_mean[i]) && cv_data$EDA_mean[i] < eda_threshold
      if(cv_low && eda_low) {
        in_event <- FALSE
      }
    }
    cv_data$in_event[i] <- in_event
  }
  return(cv_data)
}

# ===== GET EVENT RECTANGLES =====
get_event_rects <- function(night_data) {
  night_data <- night_data %>% arrange(timestamp)
  rects <- list()
  event_start <- NULL
  
  for(i in 1:nrow(night_data)) {
    if(night_data$in_event[i] && (i == 1 || !night_data$in_event[i-1])) {
      event_start <- night_data$timestamp[i]
    }
    if(!is.null(event_start) && !night_data$in_event[i] && i > 1 && night_data$in_event[i-1]) {
      rects[[length(rects)+1]] <- data.frame(
        xmin = event_start,
        xmax = night_data$timestamp[i],
        night_label = night_data$night_label[i]
      )
      event_start <- NULL
    }
    if(i == nrow(night_data) && !is.null(event_start) && night_data$in_event[i]) {
      rects[[length(rects)+1]] <- data.frame(
        xmin = event_start,
        xmax = night_data$timestamp[i],
        night_label = night_data$night_label[i]
      )
    }
  }
  if(length(rects) == 0) return(data.frame())
  return(bind_rows(rects))
}

# ===== VISUALIZATION FUNCTION =====
visualize_dual_10_10_75 <- function(pv, data = rolling_cv_trimmed) {
  
  plot_data <- data %>%
    filter(pv_label == pv) %>%
    mutate(night_label = paste("Night", night))
  
  if(nrow(plot_data) == 0) {
    return(NULL)
  }
  
  pid <- plot_data$participant[1]
  visit <- plot_data$visit[1]
  nights <- unique(plot_data$night_label)
  
  # Process each night with 10/10/75
  processed <- map_dfr(nights, function(n) {
    identify_events_dual(
      plot_data %>% filter(night_label == n),
      cv_start = 10,
      cv_end = 10,
      eda_percentile = 75
    )
  })
  
  # Get rectangles
  all_rects <- map_dfr(nights, function(n) {
    get_event_rects(processed %>% filter(night_label == n))
  })
  
  # Event counts per night
  if(nrow(all_rects) > 0) {
    event_counts <- all_rects %>%
      group_by(night_label) %>%
      summarise(n = n(), .groups = "drop")
  } else {
    event_counts <- data.frame(night_label = nights, n = 0)
  }
  
  event_counts <- data.frame(night_label = nights) %>%
    left_join(event_counts, by = "night_label") %>%
    mutate(n = ifelse(is.na(n), 0, n))
  
  subtitle_text <- paste(
    paste0(event_counts$night_label, ": ", event_counts$n, " events"),
    collapse = "   |   "
  )
  
  # Scale EDA
  cv_max <- max(processed$rolling_cv, na.rm = TRUE)
  eda_min <- min(processed$EDA_mean, na.rm = TRUE)
  eda_max <- max(processed$EDA_mean, na.rm = TRUE)
  
  processed <- processed %>%
    mutate(
      eda_scaled = ((EDA_mean - eda_min) / (eda_max - eda_min)) * cv_max,
      thresh_scaled = ((eda_threshold - eda_min) / (eda_max - eda_min)) * cv_max
    )
  
  # Build plot
  p <- ggplot(processed, aes(x = timestamp)) +
    facet_wrap(~ night_label, scales = "free_x", ncol = 1)
  
  # Shading
  if(nrow(all_rects) > 0) {
    p <- p + geom_rect(
      data = all_rects,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = "#E63946", alpha = 0.2,
      inherit.aes = FALSE
    )
  }
  
  p <- p +
    geom_line(aes(y = eda_scaled), color = "grey75", linewidth = 0.4, alpha = 0.8) +
    geom_line(aes(y = rolling_cv), color = "gray30", linewidth = 0.5) +
    geom_line(aes(y = thresh_scaled), color = "orange", linewidth = 0.8, alpha = 0.7) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "#E63946", linewidth = 0.8) +
    scale_x_datetime(date_breaks = "1 hour", date_labels = "%H:%M") +
    scale_y_continuous(
      name = "CV (%)",
      sec.axis = sec_axis(
        trans = ~ . * (eda_max - eda_min) / cv_max + eda_min,
        name = "EDA (μS)"
      )
    ) +
    labs(
      title = paste0(pid, " ", visit, " — CV 10%/10% + EDA 75th Percentile"),
      subtitle = paste0(subtitle_text, "   |   Start: CV>10%   |   End: CV<10% AND EDA<75th %ile (orange line)")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title.y.right = element_text(color = "grey60", size = 9),
      axis.text.y.right = element_text(color = "grey60", size = 8),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1.5, "lines")
    )
  
  # Save
  output_path <- paste0(output_folder, "/", pid, "_", visit, "_10-10-75.png")
  ggsave(filename = output_path, plot = p, width = 14, height = 8, dpi = 300)
  
  total_events <- if(nrow(all_rects) > 0) nrow(all_rects) else 0
  cat("Saved:", pid, visit, "-", total_events, "events\n")
  
  return(total_events)
}

# ===== RUN FOR ALL PARTICIPANTS =====
pv_list <- rolling_cv_trimmed %>%
  distinct(pv_label) %>%
  arrange(pv_label) %>%
  pull(pv_label)

cat("Creating 10/10/75 dual criteria visualizations for", length(pv_list), "visits...\n\n")

for(pv in pv_list) {
  tryCatch({
    visualize_dual_10_10_75(pv)
  }, error = function(e) {
    cat("ERROR for", pv, ":", e$message, "\n")
  })
}

cat("\n=====================================\n")
cat("All visualizations complete!\n")
cat("Saved to:", output_folder, "\n")
cat("=====================================\n")

# tried: CV on/off = 10%/5%, EDA offset = 60th percentile, 5-min rule
#### trying 10/5 with 60th percentile and 5 minute time restraint ##########
library(tidyverse)
library(lubridate)

rolling_cv_trimmed <- read.csv("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min_trimmed.csv",
                               stringsAsFactors = FALSE)
rolling_cv_trimmed$timestamp <- ymd_hms(rolling_cv_trimmed$timestamp)

# Create output folder
output_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA CV Events 10-5-60pct-5min"
dir.create(output_folder, showWarnings = FALSE)

# ===== DUAL CRITERIA EVENT DETECTION WITH 5-MINUTE CV RULE =====
identify_events_dual_5min <- function(cv_data, cv_start = 10, cv_end = 5, eda_percentile = 60) {
  cv_data <- cv_data %>% arrange(timestamp)
  eda_threshold <- quantile(cv_data$EDA_mean, eda_percentile/100, na.rm = TRUE)
  cv_data$in_event <- FALSE
  cv_data$eda_threshold <- eda_threshold
  in_event <- FALSE
  cv_below_threshold_start <- NULL  # Track when CV first dropped below 5%
  
  for(i in 1:nrow(cv_data)) {
    if(!in_event) {
      # EVENT START
      if(!is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] > cv_start) {
        in_event <- TRUE
        cv_below_threshold_start <- NULL
      }
    } else {
      # Check conditions
      cv_low <- !is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] < cv_end
      eda_low <- !is.na(cv_data$EDA_mean[i]) && cv_data$EDA_mean[i] < eda_threshold
      
      # Track when CV first drops below threshold
      if(cv_low && is.null(cv_below_threshold_start)) {
        cv_below_threshold_start <- i
      }
      
      # If CV goes back above threshold, reset the counter
      if(!cv_low) {
        cv_below_threshold_start <- NULL
      }
      
      # EVENT ENDS: CV<5% for 5+ consecutive minutes AND EDA below threshold
      if(cv_low && eda_low && !is.null(cv_below_threshold_start)) {
        minutes_below <- i - cv_below_threshold_start
        if(minutes_below >= 5) {
          in_event <- FALSE
          cv_below_threshold_start <- NULL
        }
      }
    }
    
    cv_data$in_event[i] <- in_event
  }
  return(cv_data)
}

# ===== GET EVENT RECTANGLES =====
get_event_rects <- function(night_data) {
  night_data <- night_data %>% arrange(timestamp)
  rects <- list()
  event_start <- NULL
  
  for(i in 1:nrow(night_data)) {
    if(night_data$in_event[i] && (i == 1 || !night_data$in_event[i-1])) {
      event_start <- night_data$timestamp[i]
    }
    if(!is.null(event_start) && !night_data$in_event[i] && i > 1 && night_data$in_event[i-1]) {
      rects[[length(rects)+1]] <- data.frame(
        xmin = event_start,
        xmax = night_data$timestamp[i],
        night_label = night_data$night_label[i]
      )
      event_start <- NULL
    }
    if(i == nrow(night_data) && !is.null(event_start) && night_data$in_event[i]) {
      rects[[length(rects)+1]] <- data.frame(
        xmin = event_start,
        xmax = night_data$timestamp[i],
        night_label = night_data$night_label[i]
      )
    }
  }
  if(length(rects) == 0) return(data.frame())
  return(bind_rows(rects))
}

# ===== VISUALIZATION FUNCTION =====
visualize_dual_10_5_60_5min <- function(pv, data = rolling_cv_trimmed) {
  
  plot_data <- data %>%
    filter(pv_label == pv) %>%
    mutate(night_label = paste("Night", night))
  
  if(nrow(plot_data) == 0) {
    return(NULL)
  }
  
  pid <- plot_data$participant[1]
  visit <- plot_data$visit[1]
  nights <- unique(plot_data$night_label)
  
  # Process each night
  processed <- map_dfr(nights, function(n) {
    identify_events_dual_5min(
      plot_data %>% filter(night_label == n),
      cv_start = 10,
      cv_end = 5,
      eda_percentile = 60
    )
  })
  
  # Get rectangles
  all_rects <- map_dfr(nights, function(n) {
    get_event_rects(processed %>% filter(night_label == n))
  })
  
  # Event counts per night
  if(nrow(all_rects) > 0) {
    event_counts <- all_rects %>%
      group_by(night_label) %>%
      summarise(n = n(), .groups = "drop")
  } else {
    event_counts <- data.frame(night_label = nights, n = 0)
  }
  
  event_counts <- data.frame(night_label = nights) %>%
    left_join(event_counts, by = "night_label") %>%
    mutate(n = ifelse(is.na(n), 0, n))
  
  subtitle_text <- paste(
    paste0(event_counts$night_label, ": ", event_counts$n, " events"),
    collapse = "   |   "
  )
  
  # Scale EDA
  cv_max <- max(processed$rolling_cv, na.rm = TRUE)
  eda_min <- min(processed$EDA_mean, na.rm = TRUE)
  eda_max <- max(processed$EDA_mean, na.rm = TRUE)
  
  processed <- processed %>%
    mutate(
      eda_scaled = ((EDA_mean - eda_min) / (eda_max - eda_min)) * cv_max,
      thresh_scaled = ((eda_threshold - eda_min) / (eda_max - eda_min)) * cv_max
    )
  
  # Build plot
  p <- ggplot(processed, aes(x = timestamp)) +
    facet_wrap(~ night_label, scales = "free_x", ncol = 1)
  
  # Shading
  if(nrow(all_rects) > 0) {
    p <- p + geom_rect(
      data = all_rects,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = "#E63946", alpha = 0.2,
      inherit.aes = FALSE
    )
  }
  
  p <- p +
    geom_line(aes(y = eda_scaled), color = "grey75", linewidth = 0.4, alpha = 0.8) +
    geom_line(aes(y = rolling_cv), color = "gray30", linewidth = 0.5) +
    geom_line(aes(y = thresh_scaled), color = "orange", linewidth = 0.8, alpha = 0.7) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "#E63946", linewidth = 0.8) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "#2E86AB", linewidth = 0.8) +
    scale_x_datetime(date_breaks = "1 hour", date_labels = "%H:%M") +
    scale_y_continuous(
      name = "CV (%)",
      sec.axis = sec_axis(
        trans = ~ . * (eda_max - eda_min) / cv_max + eda_min,
        name = "EDA (μS)"
      )
    ) +
    labs(
      title = paste0(pid, " ", visit, " — CV 10%/5% + EDA 60th %ile + 5min Rule"),
      subtitle = paste0(subtitle_text, "   |   Start: CV>10%   |   End: CV<5% for 5+ min AND EDA<60th %ile (orange)")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title.y.right = element_text(color = "grey60", size = 9),
      axis.text.y.right = element_text(color = "grey60", size = 8),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1.5, "lines")
    )
  
  # Save
  output_path <- paste0(output_folder, "/", pid, "_", visit, "_10-5-60-5min.png")
  ggsave(filename = output_path, plot = p, width = 14, height = 8, dpi = 300)
  
  total_events <- if(nrow(all_rects) > 0) nrow(all_rects) else 0
  cat("Saved:", pid, visit, "-", total_events, "events\n")
  
  return(total_events)
}

# ===== RUN FOR ALL PARTICIPANTS =====
pv_list <- rolling_cv_trimmed %>%
  distinct(pv_label) %>%
  arrange(pv_label) %>%
  pull(pv_label)

cat("Creating 10/5/60/5min dual criteria visualizations for", length(pv_list), "visits...\n\n")

for(pv in pv_list) {
  tryCatch({
    visualize_dual_10_5_60_5min(pv)
  }, error = function(e) {
    cat("ERROR for", pv, ":", e$message, "\n")
  })
}

cat("\n=====================================\n")
cat("All visualizations complete!\n")
cat("Saved to:", output_folder, "\n")
cat("=====================================\n")

# tried: CV on/off = 10%/10%, EDA offset = 60th percentile, 5-min rule
############10/10 60 5min##################
library(tidyverse)
library(lubridate)

rolling_cv_trimmed <- read.csv("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min_trimmed.csv",
                               stringsAsFactors = FALSE)
rolling_cv_trimmed$timestamp <- ymd_hms(rolling_cv_trimmed$timestamp)

# Create output folder
output_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA CV Events 10-10-60pct-5min"
dir.create(output_folder, showWarnings = FALSE)

# ===== DUAL CRITERIA EVENT DETECTION WITH 5-MINUTE RULE =====
identify_events_dual_5min <- function(cv_data, cv_start = 10, cv_end = 10, eda_percentile = 60) {
  cv_data <- cv_data %>% arrange(timestamp)
  eda_threshold <- quantile(cv_data$EDA_mean, eda_percentile/100, na.rm = TRUE)
  cv_data$in_event <- FALSE
  cv_data$eda_threshold <- eda_threshold
  in_event <- FALSE
  both_conditions_start <- NULL  # Track when BOTH conditions became true
  
  for(i in 1:nrow(cv_data)) {
    if(!in_event) {
      # EVENT START
      if(!is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] > cv_start) {
        in_event <- TRUE
        both_conditions_start <- NULL
      }
    } else {
      # Check BOTH conditions
      cv_low <- !is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] < cv_end
      eda_low <- !is.na(cv_data$EDA_mean[i]) && cv_data$EDA_mean[i] < eda_threshold
      both_conditions_met <- cv_low && eda_low
      
      # Track when BOTH conditions first become true
      if(both_conditions_met && is.null(both_conditions_start)) {
        both_conditions_start <- i
      }
      
      # If either condition fails, reset the counter
      if(!both_conditions_met) {
        both_conditions_start <- NULL
      }
      
      # EVENT ENDS: BOTH conditions true for 5+ consecutive minutes
      if(both_conditions_met && !is.null(both_conditions_start)) {
        minutes_both_low <- i - both_conditions_start
        if(minutes_both_low >= 5) {
          in_event <- FALSE
          both_conditions_start <- NULL
        }
      }
    }
    
    cv_data$in_event[i] <- in_event
  }
  return(cv_data)
}

# ===== GET EVENT RECTANGLES =====
get_event_rects <- function(night_data) {
  night_data <- night_data %>% arrange(timestamp)
  rects <- list()
  event_start <- NULL
  
  for(i in 1:nrow(night_data)) {
    if(night_data$in_event[i] && (i == 1 || !night_data$in_event[i-1])) {
      event_start <- night_data$timestamp[i]
    }
    if(!is.null(event_start) && !night_data$in_event[i] && i > 1 && night_data$in_event[i-1]) {
      rects[[length(rects)+1]] <- data.frame(
        xmin = event_start,
        xmax = night_data$timestamp[i],
        night_label = night_data$night_label[i]
      )
      event_start <- NULL
    }
    if(i == nrow(night_data) && !is.null(event_start) && night_data$in_event[i]) {
      rects[[length(rects)+1]] <- data.frame(
        xmin = event_start,
        xmax = night_data$timestamp[i],
        night_label = night_data$night_label[i]
      )
    }
  }
  if(length(rects) == 0) return(data.frame())
  return(bind_rows(rects))
}

# ===== VISUALIZATION FUNCTION =====
visualize_dual_10_10_60_5min <- function(pv, data = rolling_cv_trimmed) {
  
  plot_data <- data %>%
    filter(pv_label == pv) %>%
    mutate(night_label = paste("Night", night))
  
  if(nrow(plot_data) == 0) {
    return(NULL)
  }
  
  pid <- plot_data$participant[1]
  visit <- plot_data$visit[1]
  nights <- unique(plot_data$night_label)
  
  # Process each night
  processed <- map_dfr(nights, function(n) {
    identify_events_dual_5min(
      plot_data %>% filter(night_label == n),
      cv_start = 10,
      cv_end = 10,
      eda_percentile = 60
    )
  })
  
  # Get rectangles
  all_rects <- map_dfr(nights, function(n) {
    get_event_rects(processed %>% filter(night_label == n))
  })
  
  # Event counts per night
  if(nrow(all_rects) > 0) {
    event_counts <- all_rects %>%
      group_by(night_label) %>%
      summarise(n = n(), .groups = "drop")
  } else {
    event_counts <- data.frame(night_label = nights, n = 0)
  }
  
  event_counts <- data.frame(night_label = nights) %>%
    left_join(event_counts, by = "night_label") %>%
    mutate(n = ifelse(is.na(n), 0, n))
  
  subtitle_text <- paste(
    paste0(event_counts$night_label, ": ", event_counts$n, " events"),
    collapse = "   |   "
  )
  
  # Scale EDA
  cv_max <- max(processed$rolling_cv, na.rm = TRUE)
  eda_min <- min(processed$EDA_mean, na.rm = TRUE)
  eda_max <- max(processed$EDA_mean, na.rm = TRUE)
  
  processed <- processed %>%
    mutate(
      eda_scaled = ((EDA_mean - eda_min) / (eda_max - eda_min)) * cv_max,
      thresh_scaled = ((eda_threshold - eda_min) / (eda_max - eda_min)) * cv_max
    )
  
  # Build plot
  p <- ggplot(processed, aes(x = timestamp)) +
    facet_wrap(~ night_label, scales = "free_x", ncol = 1)
  
  # Shading
  if(nrow(all_rects) > 0) {
    p <- p + geom_rect(
      data = all_rects,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = "#E63946", alpha = 0.2,
      inherit.aes = FALSE
    )
  }
  
  p <- p +
    geom_line(aes(y = eda_scaled), color = "grey75", linewidth = 0.4, alpha = 0.8) +
    geom_line(aes(y = rolling_cv), color = "gray30", linewidth = 0.5) +
    geom_line(aes(y = thresh_scaled), color = "orange", linewidth = 0.8, alpha = 0.7) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "#E63946", linewidth = 0.8) +
    scale_x_datetime(date_breaks = "1 hour", date_labels = "%H:%M") +
    scale_y_continuous(
      name = "CV (%)",
      sec.axis = sec_axis(
        trans = ~ . * (eda_max - eda_min) / cv_max + eda_min,
        name = "EDA (μS)"
      )
    ) +
    labs(
      title = paste0(pid, " ", visit, " — CV 10%/10% + EDA 60th %ile (5min sustained)"),
      subtitle = paste0(subtitle_text, "   |   Start: CV>10%   |   End: CV<10% AND EDA<60th %ile (orange) for 5+ min")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title.y.right = element_text(color = "grey60", size = 9),
      axis.text.y.right = element_text(color = "grey60", size = 8),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1.5, "lines")
    )
  
  # Save
  output_path <- paste0(output_folder, "/", pid, "_", visit, "_10-10-60-5min.png")
  ggsave(filename = output_path, plot = p, width = 14, height = 8, dpi = 300)
  
  total_events <- if(nrow(all_rects) > 0) nrow(all_rects) else 0
  cat("Saved:", pid, visit, "-", total_events, "events\n")
  
  return(total_events)
}

# ===== RUN FOR ALL PARTICIPANTS =====
pv_list <- rolling_cv_trimmed %>%
  distinct(pv_label) %>%
  arrange(pv_label) %>%
  pull(pv_label)

cat("Creating 10/10/60/5min dual criteria visualizations for", length(pv_list), "visits...\n\n")

for(pv in pv_list) {
  tryCatch({
    visualize_dual_10_10_60_5min(pv)
  }, error = function(e) {
    cat("ERROR for", pv, ":", e$message, "\n")
  })
}

cat("\n=====================================\n")
cat("All visualizations complete!\n")
cat("Saved to:", output_folder, "\n")
cat("=====================================\n")

# FINAL — CV on/off = 10%/10%, EDA offset = 75th percentile, 5-min rule.
###I tried many more but stupidly re-wrote the code rather than starting new blocks,
###List of all criteria tried is in the paper###
# this is what's in the manuscript. also generates the trace plots
# I picked a few of for Figure 1.
################10/10 75 5min##########
library(tidyverse)
library(lubridate)

rolling_cv_trimmed <- read.csv("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min_trimmed.csv",
                               stringsAsFactors = FALSE)
rolling_cv_trimmed$timestamp <- ymd_hms(rolling_cv_trimmed$timestamp)

# Create output folder
output_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA CV Events 10-10-75pct-5min"
dir.create(output_folder, showWarnings = FALSE)

# ===== DUAL CRITERIA EVENT DETECTION WITH 5-MINUTE RULE =====
identify_events_dual_5min <- function(cv_data, cv_start = 10, cv_end = 10, eda_percentile = 75) {
  cv_data <- cv_data %>% arrange(timestamp)
  eda_threshold <- quantile(cv_data$EDA_mean, eda_percentile/100, na.rm = TRUE)
  cv_data$in_event <- FALSE
  cv_data$eda_threshold <- eda_threshold
  in_event <- FALSE
  both_conditions_start <- NULL  # Track when BOTH conditions became true
  
  for(i in 1:nrow(cv_data)) {
    if(!in_event) {
      # EVENT START
      if(!is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] > cv_start) {
        in_event <- TRUE
        both_conditions_start <- NULL
      }
    } else {
      # Check BOTH conditions
      cv_low <- !is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] < cv_end
      eda_low <- !is.na(cv_data$EDA_mean[i]) && cv_data$EDA_mean[i] < eda_threshold
      both_conditions_met <- cv_low && eda_low
      
      # Track when BOTH conditions first become true
      if(both_conditions_met && is.null(both_conditions_start)) {
        both_conditions_start <- i
      }
      
      # If either condition fails, reset the counter
      if(!both_conditions_met) {
        both_conditions_start <- NULL
      }
      
      # EVENT ENDS: BOTH conditions true for 5+ consecutive minutes
      if(both_conditions_met && !is.null(both_conditions_start)) {
        minutes_both_low <- i - both_conditions_start
        if(minutes_both_low >= 5) {
          in_event <- FALSE
          both_conditions_start <- NULL
        }
      }
    }
    
    cv_data$in_event[i] <- in_event
  }
  return(cv_data)
}

# ===== GET EVENT RECTANGLES =====
get_event_rects <- function(night_data) {
  night_data <- night_data %>% arrange(timestamp)
  rects <- list()
  event_start <- NULL
  
  for(i in 1:nrow(night_data)) {
    if(night_data$in_event[i] && (i == 1 || !night_data$in_event[i-1])) {
      event_start <- night_data$timestamp[i]
    }
    if(!is.null(event_start) && !night_data$in_event[i] && i > 1 && night_data$in_event[i-1]) {
      rects[[length(rects)+1]] <- data.frame(
        xmin = event_start,
        xmax = night_data$timestamp[i],
        night_label = night_data$night_label[i]
      )
      event_start <- NULL
    }
    if(i == nrow(night_data) && !is.null(event_start) && night_data$in_event[i]) {
      rects[[length(rects)+1]] <- data.frame(
        xmin = event_start,
        xmax = night_data$timestamp[i],
        night_label = night_data$night_label[i]
      )
    }
  }
  if(length(rects) == 0) return(data.frame())
  return(bind_rows(rects))
}

# ===== VISUALIZATION FUNCTION =====
visualize_dual_10_10_75_5min <- function(pv, data = rolling_cv_trimmed) {
  
  plot_data <- data %>%
    filter(pv_label == pv) %>%
    mutate(night_label = paste("Night", night))
  
  if(nrow(plot_data) == 0) {
    return(NULL)
  }
  
  pid <- plot_data$participant[1]
  visit <- plot_data$visit[1]
  nights <- unique(plot_data$night_label)
  
  # Process each night
  processed <- map_dfr(nights, function(n) {
    identify_events_dual_5min(
      plot_data %>% filter(night_label == n),
      cv_start = 10,
      cv_end = 10,
      eda_percentile = 75
    )
  })
  
  # Get rectangles
  all_rects <- map_dfr(nights, function(n) {
    get_event_rects(processed %>% filter(night_label == n))
  })
  
  # Event counts per night
  if(nrow(all_rects) > 0) {
    event_counts <- all_rects %>%
      group_by(night_label) %>%
      summarise(n = n(), .groups = "drop")
  } else {
    event_counts <- data.frame(night_label = nights, n = 0)
  }
  
  event_counts <- data.frame(night_label = nights) %>%
    left_join(event_counts, by = "night_label") %>%
    mutate(n = ifelse(is.na(n), 0, n))
  
  subtitle_text <- paste(
    paste0(event_counts$night_label, ": ", event_counts$n, " events"),
    collapse = "   |   "
  )
  
  # Scale EDA
  cv_max <- max(processed$rolling_cv, na.rm = TRUE)
  eda_min <- min(processed$EDA_mean, na.rm = TRUE)
  eda_max <- max(processed$EDA_mean, na.rm = TRUE)
  
  processed <- processed %>%
    mutate(
      eda_scaled = ((EDA_mean - eda_min) / (eda_max - eda_min)) * cv_max,
      thresh_scaled = ((eda_threshold - eda_min) / (eda_max - eda_min)) * cv_max
    )
  
  # Build plot
  p <- ggplot(processed, aes(x = timestamp)) +
    facet_wrap(~ night_label, scales = "free_x", ncol = 1)
  
  # Shading
  if(nrow(all_rects) > 0) {
    p <- p + geom_rect(
      data = all_rects,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = "#E63946", alpha = 0.2,
      inherit.aes = FALSE
    )
  }
  
  p <- p +
    geom_line(aes(y = eda_scaled), color = "grey75", linewidth = 0.4, alpha = 0.8) +
    geom_line(aes(y = rolling_cv), color = "gray30", linewidth = 0.5) +
    geom_line(aes(y = thresh_scaled), color = "orange", linewidth = 0.8, alpha = 0.7) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "#E63946", linewidth = 0.8) +
    scale_x_datetime(date_breaks = "1 hour", date_labels = "%H:%M") +
    scale_y_continuous(
      name = "CV (%)",
      sec.axis = sec_axis(
        trans = ~ . * (eda_max - eda_min) / cv_max + eda_min,
        name = "EDA (μS)"
      )
    ) +
    labs(
      title = paste0(pid, " ", visit, " — CV 10%/10% + EDA 75th %ile (5min sustained)"),
      subtitle = paste0(subtitle_text, "   |   Start: CV>10%   |   End: CV<10% AND EDA<75th %ile (orange) for 5+ min")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title.y.right = element_text(color = "grey60", size = 9),
      axis.text.y.right = element_text(color = "grey60", size = 8),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1.5, "lines")
    )
  
  # Save
  output_path <- paste0(output_folder, "/", pid, "_", visit, "_10-10-75-5min.png")
  ggsave(filename = output_path, plot = p, width = 14, height = 8, dpi = 300)
  
  total_events <- if(nrow(all_rects) > 0) nrow(all_rects) else 0
  cat("Saved:", pid, visit, "-", total_events, "events\n")
  
  return(total_events)
}

# ===== RUN FOR ALL PARTICIPANTS =====
pv_list <- rolling_cv_trimmed %>%
  distinct(pv_label) %>%
  arrange(pv_label) %>%
  pull(pv_label)

cat("Creating 10/10/75/5min dual criteria visualizations for", length(pv_list), "visits...\n\n")

for(pv in pv_list) {
  tryCatch({
    visualize_dual_10_10_75_5min(pv)
  }, error = function(e) {
    cat("ERROR for", pv, ":", e$message, "\n")
  })
}

cat("\n=====================================\n")
cat("All visualizations complete!\n")
cat("Saved to:", output_folder, "\n")
cat("=====================================\n")


################################################################################
# VALIDATING AGAINST PARTICIPANT-MARKED EVENTS
################################################################################

# checking the algorithm against events participants actually flagged,
# false positives/negatives etc. - this is what backs up the accuracy
# claim in Methods.
############################# adding event markers#######################
library(tidyverse)
library(lubridate)
library(readxl)  # Need this for Excel files

rolling_cv_trimmed <- read.csv("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min_trimmed.csv",
                               stringsAsFactors = FALSE)
rolling_cv_trimmed$timestamp <- ymd_hms(rolling_cv_trimmed$timestamp)

# ===== STEP 1: READ AND ADD EXPERT MARKERS TO DATASET =====

read_all_expert_markers <- function() {
  
  all_markers <- data.frame()
  
  # Get all unique participant/visit combinations
  pv_list <- rolling_cv_trimmed %>%
    distinct(participant, visit)
  
  for(i in 1:nrow(pv_list)) {
    pid <- pv_list$participant[i]
    visit_label <- pv_list$visit[i]
    
    # Remove "S" from participant for folder path
    participant_num <- gsub("S", "", pid)
    
    # Construct path - CHANGE TO .xlsx
    event_marker_path <- paste0("Y:/Busa/Busa_AIMS/ShreyerDissertation/Event Marker Data/S",
                                participant_num, "/", visit_label, 
                                "/hf_subjective_objective.xlsx")
    
    # Check if file exists
    if(!file.exists(event_marker_path)) {
      # Silently skip - most won't have markers
      next
    }
    
    # Read the Excel file
    tryCatch({
      markers <- read_excel(event_marker_path)
      
      # Parse datetime with seconds
      markers$datetime <- mdy_hms(markers$datetime)
      
      # Round to nearest minute to match our dataset
      markers$timestamp_minute <- floor_date(markers$datetime, unit = "minute")
      
      # Filter for expert_identified = 1 or 3
      markers <- markers %>%
        filter(expert_identified %in% c(1, 3)) %>%
        select(timestamp_minute, expert_identified) %>%
        mutate(
          participant = pid,
          visit = visit_label,
          pv_label = paste0(pid, visit_label),
          event_type = ifelse(expert_identified == 1, "HF", "NS")
        )
      
      if(nrow(markers) > 0) {
        all_markers <- bind_rows(all_markers, markers)
        cat("Found", nrow(markers), "expert markers for", pid, visit_label, "\n")
      }
      
    }, error = function(e) {
      cat("Error reading", pid, visit_label, ":", e$message, "\n")
    })
  }
  
  return(all_markers)
}

# Read all expert markers
cat("Reading expert marker files...\n")
expert_markers <- read_all_expert_markers()

cat("\n=== Expert Markers Summary ===\n")
cat("Total expert markers found:", nrow(expert_markers), "\n")
cat("Hot flashes (type 1):", sum(expert_markers$expert_identified == 1), "\n")
cat("Night sweats (type 3):", sum(expert_markers$expert_identified == 3), "\n")
cat("Participants with markers:", n_distinct(expert_markers$pv_label), "\n")

# Add expert markers to the main dataset
rolling_cv_with_markers <- rolling_cv_trimmed %>%
  left_join(
    expert_markers %>% select(timestamp_minute, pv_label, expert_identified, event_type),
    by = c("timestamp" = "timestamp_minute", "pv_label" = "pv_label")
  )

cat("\nDataset rows:", nrow(rolling_cv_with_markers), "\n")
cat("Rows with expert markers:", sum(!is.na(rolling_cv_with_markers$expert_identified)), "\n")

# Save the enhanced dataset
output_path <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min_trimmed_with_markers.csv"
write.csv(rolling_cv_with_markers, output_path, row.names = FALSE)
cat("\nSaved dataset with markers to:", output_path, "\n")

# Show some examples
cat("\n=== Example Markers ===\n")
print(rolling_cv_with_markers %>% 
        filter(!is.na(expert_identified)) %>%
        select(pv_label, timestamp, EDA_mean, rolling_cv, expert_identified, event_type) %>%
        head(10))
#################adding event markers into visualuizations################################
library(tidyverse)
library(lubridate)

# Load the dataset with markers
rolling_cv_with_markers <- read.csv("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min_trimmed_with_markers.csv",
                                    stringsAsFactors = FALSE)
rolling_cv_with_markers$timestamp <- ymd_hms(rolling_cv_with_markers$timestamp)

# Create output folder
output_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA CV Events With Expert Markers"
dir.create(output_folder, showWarnings = FALSE)

# ===== DUAL CRITERIA EVENT DETECTION =====
identify_events_dual_5min <- function(cv_data, cv_start = 10, cv_end = 10, eda_percentile = 75) {
  cv_data <- cv_data %>% arrange(timestamp)
  eda_threshold <- quantile(cv_data$EDA_mean, eda_percentile/100, na.rm = TRUE)
  cv_data$in_event <- FALSE
  cv_data$eda_threshold <- eda_threshold
  in_event <- FALSE
  both_conditions_start <- NULL
  
  for(i in 1:nrow(cv_data)) {
    if(!in_event) {
      if(!is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] > cv_start) {
        in_event <- TRUE
        both_conditions_start <- NULL
      }
    } else {
      cv_low <- !is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] < cv_end
      eda_low <- !is.na(cv_data$EDA_mean[i]) && cv_data$EDA_mean[i] < eda_threshold
      both_conditions_met <- cv_low && eda_low
      
      if(both_conditions_met && is.null(both_conditions_start)) {
        both_conditions_start <- i
      }
      if(!both_conditions_met) {
        both_conditions_start <- NULL
      }
      if(both_conditions_met && !is.null(both_conditions_start)) {
        minutes_both_low <- i - both_conditions_start
        if(minutes_both_low >= 5) {
          in_event <- FALSE
          both_conditions_start <- NULL
        }
      }
    }
    cv_data$in_event[i] <- in_event
  }
  return(cv_data)
}

# ===== GET EVENT RECTANGLES =====
get_event_rects <- function(night_data) {
  night_data <- night_data %>% arrange(timestamp)
  rects <- list()
  event_start <- NULL
  
  for(i in 1:nrow(night_data)) {
    if(night_data$in_event[i] && (i == 1 || !night_data$in_event[i-1])) {
      event_start <- night_data$timestamp[i]
    }
    if(!is.null(event_start) && !night_data$in_event[i] && i > 1 && night_data$in_event[i-1]) {
      rects[[length(rects)+1]] <- data.frame(
        xmin = event_start,
        xmax = night_data$timestamp[i],
        night_label = night_data$night_label[i]
      )
      event_start <- NULL
    }
    if(i == nrow(night_data) && !is.null(event_start) && night_data$in_event[i]) {
      rects[[length(rects)+1]] <- data.frame(
        xmin = event_start,
        xmax = night_data$timestamp[i],
        night_label = night_data$night_label[i]
      )
    }
  }
  if(length(rects) == 0) return(data.frame())
  return(bind_rows(rects))
}

# ===== VISUALIZATION WITH EXPERT MARKERS =====
visualize_with_expert_markers <- function(pv, data = rolling_cv_with_markers) {
  
  plot_data <- data %>%
    filter(pv_label == pv) %>%
    mutate(night_label = paste("Night", night))
  
  if(nrow(plot_data) == 0) return(NULL)
  
  pid <- plot_data$participant[1]
  visit <- plot_data$visit[1]
  nights <- unique(plot_data$night_label)
  
  # Get expert markers for this participant
  expert_markers <- plot_data %>%
    filter(!is.na(expert_identified)) %>%
    select(timestamp, expert_identified, event_type, night_label)
  
  # Process each night
  processed <- map_dfr(nights, function(n) {
    identify_events_dual_5min(
      plot_data %>% filter(night_label == n),
      cv_start = 10, cv_end = 10, eda_percentile = 75
    )
  })
  
  # Get rectangles
  all_rects <- map_dfr(nights, function(n) {
    get_event_rects(processed %>% filter(night_label == n))
  })
  
  # Count events per night
  if(nrow(all_rects) > 0) {
    event_counts <- all_rects %>%
      group_by(night_label) %>%
      summarise(n_algo = n(), .groups = "drop")
  } else {
    event_counts <- data.frame(night_label = nights, n_algo = 0)
  }
  
  # Count expert markers per night
  if(nrow(expert_markers) > 0) {
    expert_counts <- expert_markers %>%
      group_by(night_label) %>%
      summarise(
        n_hf = sum(expert_identified == 1),
        n_ns = sum(expert_identified == 3),
        .groups = "drop"
      )
    
    event_counts <- data.frame(night_label = nights) %>%
      left_join(event_counts, by = "night_label") %>%
      left_join(expert_counts, by = "night_label") %>%
      mutate(
        n_algo = ifelse(is.na(n_algo), 0, n_algo),
        n_hf = ifelse(is.na(n_hf), 0, n_hf),
        n_ns = ifelse(is.na(n_ns), 0, n_ns)
      )
  } else {
    event_counts <- data.frame(night_label = nights) %>%
      left_join(event_counts, by = "night_label") %>%
      mutate(
        n_algo = ifelse(is.na(n_algo), 0, n_algo),
        n_hf = 0,
        n_ns = 0
      )
  }
  
  subtitle_text <- paste(
    sapply(1:nrow(event_counts), function(i) {
      paste0(event_counts$night_label[i], ": ", event_counts$n_algo[i], 
             " algo | ", event_counts$n_hf[i], " HF | ", event_counts$n_ns[i], " NS")
    }),
    collapse = "   |   "
  )
  
  # Scale EDA
  cv_max <- max(processed$rolling_cv, na.rm = TRUE)
  eda_min <- min(processed$EDA_mean, na.rm = TRUE)
  eda_max <- max(processed$EDA_mean, na.rm = TRUE)
  
  processed <- processed %>%
    mutate(
      eda_scaled = ((EDA_mean - eda_min) / (eda_max - eda_min)) * cv_max,
      thresh_scaled = ((eda_threshold - eda_min) / (eda_max - eda_min)) * cv_max
    )
  
  # Build plot
  p <- ggplot(processed, aes(x = timestamp)) +
    facet_wrap(~ night_label, scales = "free_x", ncol = 1)
  
  # Algorithm-detected events (pink shading)
  if(nrow(all_rects) > 0) {
    p <- p + geom_rect(
      data = all_rects,
      aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
      fill = "#E63946", alpha = 0.15,
      inherit.aes = FALSE
    )
  }
  
  # EDA/CV lines
  p <- p +
    geom_line(aes(y = eda_scaled), color = "grey75", linewidth = 0.4, alpha = 0.8) +
    geom_line(aes(y = rolling_cv), color = "gray30", linewidth = 0.5) +
    geom_line(aes(y = thresh_scaled), color = "orange", linewidth = 0.8, alpha = 0.7) +
    geom_hline(yintercept = 10, linetype = "dashed", color = "#E63946", linewidth = 0.8)
  
  # Add expert markers as vertical lines
  if(nrow(expert_markers) > 0) {
    # Hot flashes (type 1) - RED vertical lines
    hf_markers <- expert_markers %>% filter(expert_identified == 1)
    if(nrow(hf_markers) > 0) {
      p <- p + geom_vline(
        data = hf_markers,
        aes(xintercept = timestamp),
        color = "red", linewidth = 1.2, alpha = 0.9
      )
    }
    
    # Night sweats (type 3) - PURPLE vertical lines
    ns_markers <- expert_markers %>% filter(expert_identified == 3)
    if(nrow(ns_markers) > 0) {
      p <- p + geom_vline(
        data = ns_markers,
        aes(xintercept = timestamp),
        color = "purple", linewidth = 1.2, alpha = 0.9
      )
    }
  }
  
  p <- p +
    scale_x_datetime(date_breaks = "1 hour", date_labels = "%H:%M") +
    scale_y_continuous(
      name = "CV (%)",
      sec.axis = sec_axis(
        trans = ~ . * (eda_max - eda_min) / cv_max + eda_min,
        name = "EDA (μS)"
      )
    ) +
    labs(
      title = paste0(pid, " ", visit, " — Algorithm vs Expert Markers"),
      subtitle = paste0(subtitle_text, 
                        if(nrow(expert_markers) > 0) 
                          "\nRed line = Expert HF | Purple line = Expert NS | Pink shading = Algorithm detected" 
                        else "")
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.title.y.right = element_text(color = "grey60", size = 9),
      axis.text.y.right = element_text(color = "grey60", size = 8),
      strip.text = element_text(face = "bold", size = 11),
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 9, color = "gray40"),
      panel.grid.minor = element_blank(),
      panel.spacing = unit(1.5, "lines")
    )
  
  # Save
  output_path <- paste0(output_folder, "/", pid, "_", visit, "_with_expert.png")
  ggsave(filename = output_path, plot = p, width = 14, height = 8, dpi = 300)
  
  cat("Saved:", pid, visit, "- Algo:", sum(event_counts$n_algo), 
      "| HF:", sum(event_counts$n_hf), "| NS:", sum(event_counts$n_ns), "\n")
}

# ===== RUN FOR ALL PARTICIPANTS =====
pv_list <- rolling_cv_with_markers %>%
  distinct(pv_label) %>%
  arrange(pv_label) %>%
  pull(pv_label)

cat("Creating visualizations with expert markers for", length(pv_list), "visits...\n\n")

for(pv in pv_list) {
  tryCatch({
    visualize_with_expert_markers(pv)
  }, error = function(e) {
    cat("ERROR for", pv, ":", e$message, "\n")
  })
}

cat("\n=====================================\n")
cat("All visualizations complete!\n")
cat("=====================================\n")
###########comparing counts of algo vs expert ID ###############
library(tidyverse)
library(lubridate)

rolling_cv_with_markers <- read.csv("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min_trimmed_with_markers.csv",
                                    stringsAsFactors = FALSE)
rolling_cv_with_markers$timestamp <- ymd_hms(rolling_cv_with_markers$timestamp)

# ===== EVENT DETECTION FUNCTIONS =====
identify_events_dual_5min <- function(cv_data, cv_start = 10, cv_end = 10, eda_percentile = 75) {
  cv_data <- cv_data %>% arrange(timestamp)
  eda_threshold <- quantile(cv_data$EDA_mean, eda_percentile/100, na.rm = TRUE)
  cv_data$in_event <- FALSE
  cv_data$eda_threshold <- eda_threshold
  in_event <- FALSE
  both_conditions_start <- NULL
  
  for(i in 1:nrow(cv_data)) {
    if(!in_event) {
      if(!is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] > cv_start) {
        in_event <- TRUE
        both_conditions_start <- NULL
      }
    } else {
      cv_low <- !is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] < cv_end
      eda_low <- !is.na(cv_data$EDA_mean[i]) && cv_data$EDA_mean[i] < eda_threshold
      both_conditions_met <- cv_low && eda_low
      
      if(both_conditions_met && is.null(both_conditions_start)) {
        both_conditions_start <- i
      }
      if(!both_conditions_met) {
        both_conditions_start <- NULL
      }
      if(both_conditions_met && !is.null(both_conditions_start)) {
        minutes_both_low <- i - both_conditions_start
        if(minutes_both_low >= 5) {
          in_event <- FALSE
          both_conditions_start <- NULL
        }
      }
    }
    cv_data$in_event[i] <- in_event
  }
  return(cv_data)
}

get_event_rects <- function(night_data) {
  night_data <- night_data %>% arrange(timestamp)
  rects <- list()
  event_start <- NULL
  
  for(i in 1:nrow(night_data)) {
    if(night_data$in_event[i] && (i == 1 || !night_data$in_event[i-1])) {
      event_start <- night_data$timestamp[i]
    }
    if(!is.null(event_start) && !night_data$in_event[i] && i > 1 && night_data$in_event[i-1]) {
      rects[[length(rects)+1]] <- data.frame(
        xmin = event_start,
        xmax = night_data$timestamp[i]
      )
      event_start <- NULL
    }
    if(i == nrow(night_data) && !is.null(event_start) && night_data$in_event[i]) {
      rects[[length(rects)+1]] <- data.frame(
        xmin = event_start,
        xmax = night_data$timestamp[i]
      )
    }
  }
  if(length(rects) == 0) return(data.frame())
  return(bind_rows(rects))
}

# ===== SIMPLE COUNT COMPARISON =====
compare_counts <- function(pv, data = rolling_cv_with_markers) {
  
  plot_data <- data %>%
    filter(pv_label == pv) %>%
    mutate(night_label = paste("Night", night))
  
  if(nrow(plot_data) == 0) return(NULL)
  
  nights <- unique(plot_data$night_label)
  
  # Count expert markers (1 and 3 combined)
  n_expert_hf <- sum(!is.na(plot_data$expert_identified) & plot_data$expert_identified == 1)
  n_expert_ns <- sum(!is.na(plot_data$expert_identified) & plot_data$expert_identified == 3)
  n_expert_total <- n_expert_hf + n_expert_ns
  
  # Detect algorithm events
  processed <- map_dfr(nights, function(n) {
    identify_events_dual_5min(
      plot_data %>% filter(night_label == n),
      cv_start = 10, cv_end = 10, eda_percentile = 75
    )
  })
  
  all_rects <- get_event_rects(processed)
  n_algo <- if(nrow(all_rects) > 0) nrow(all_rects) else 0
  
  return(data.frame(
    pv_label = pv,
    n_algo_events = n_algo,
    n_expert_hf = n_expert_hf,
    n_expert_ns = n_expert_ns,
    n_expert_total = n_expert_total,
    difference = n_algo - n_expert_total
  ))
}

# ===== RUN FOR ALL PARTICIPANTS =====
pv_list <- rolling_cv_with_markers %>%
  distinct(pv_label) %>%
  arrange(pv_label) %>%
  pull(pv_label)

cat("Comparing algorithm vs expert counts for all participants...\n\n")

comparison_results <- map_dfr(pv_list, compare_counts)

# ===== SUMMARY STATISTICS =====
cat("=====================================\n")
cat("   ALGORITHM vs EXPERT COMPARISON\n")
cat("=====================================\n\n")

# Overall totals
cat("OVERALL TOTALS:\n")
cat("Algorithm events:", sum(comparison_results$n_algo_events), "\n")
cat("Expert HF (type 1):", sum(comparison_results$n_expert_hf), "\n")
cat("Expert NS (type 3):", sum(comparison_results$n_expert_ns), "\n")
cat("Expert total:", sum(comparison_results$n_expert_total), "\n")
cat("Difference (algo - expert):", sum(comparison_results$difference), "\n\n")

# Average per participant
cat("AVERAGE PER PARTICIPANT/VISIT:\n")
cat("Algo events:", round(mean(comparison_results$n_algo_events), 1), "\n")
cat("Expert events:", round(mean(comparison_results$n_expert_total), 1), "\n\n")

# Participants with expert markers
with_markers <- comparison_results %>% filter(n_expert_total > 0)
cat("PARTICIPANTS WITH EXPERT MARKERS (n=", nrow(with_markers), "):\n", sep="")
cat("Algo events:", sum(with_markers$n_algo_events), "\n")
cat("Expert events:", sum(with_markers$n_expert_total), "\n")
cat("Average algo per participant:", round(mean(with_markers$n_algo_events), 1), "\n")
cat("Average expert per participant:", round(mean(with_markers$n_expert_total), 1), "\n\n")

# Correlation
cat("CORRELATION:\n")
cat("Algo vs Expert counts: r =", round(cor(with_markers$n_algo_events, with_markers$n_expert_total), 3), "\n\n")

cat("=====================================\n\n")

# Save full comparison
write.csv(comparison_results, 
          "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/algo_vs_expert_counts.csv",
          row.names = FALSE)

cat("Full comparison saved to: algo_vs_expert_counts.csv\n\n")

# Show participants sorted by expert events
cat("=== Participants with Most Expert Events ===\n")
print(with_markers %>% 
        arrange(desc(n_expert_total)) %>%
        select(pv_label, n_algo_events, n_expert_hf, n_expert_ns, n_expert_total, difference) %>%
        head(10))

cat("\n=== Participants with Biggest Discrepancies ===\n")
print(with_markers %>% 
        arrange(desc(abs(difference))) %>%
        select(pv_label, n_algo_events, n_expert_total, difference) %>%
        head(10))


################################################################################
# PULLING EVENT METRICS — duration, AUC, time to peak, recovery time
################################################################################

########running metrics on events ############
library(tidyverse)
library(lubridate)

rolling_cv_trimmed <- read.csv("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min_trimmed.csv",
                               stringsAsFactors = FALSE)
rolling_cv_trimmed$timestamp <- ymd_hms(rolling_cv_trimmed$timestamp)

# ===== EVENT DETECTION =====
identify_events_dual_5min <- function(cv_data, cv_start = 10, cv_end = 10, eda_percentile = 75) {
  cv_data <- cv_data %>% arrange(timestamp)
  mhos_threshold <- quantile(cv_data$mhos, eda_percentile/100, na.rm = TRUE)
  cv_data$in_event <- FALSE
  cv_data$mhos_threshold <- mhos_threshold
  in_event <- FALSE
  both_conditions_start <- NULL
  
  for(i in 1:nrow(cv_data)) {
    if(!in_event) {
      if(!is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] > cv_start) {
        in_event <- TRUE
        both_conditions_start <- NULL
      }
    } else {
      cv_low <- !is.na(cv_data$rolling_cv[i]) && cv_data$rolling_cv[i] < cv_end
      mhos_low <- !is.na(cv_data$mhos[i]) && cv_data$mhos[i] < mhos_threshold
      both_conditions_met <- cv_low && mhos_low
      
      if(both_conditions_met && is.null(both_conditions_start)) {
        both_conditions_start <- i
      }
      if(!both_conditions_met) {
        both_conditions_start <- NULL
      }
      if(both_conditions_met && !is.null(both_conditions_start)) {
        minutes_both_low <- i - both_conditions_start
        if(minutes_both_low >= 5) {
          in_event <- FALSE
          both_conditions_start <- NULL
        }
      }
    }
    cv_data$in_event[i] <- in_event
  }
  return(cv_data)
}

# ===== EXTRACT EVENTS WITH ENHANCED METRICS =====
extract_events_with_metrics <- function(cv_data) {
  cv_data <- cv_data %>% arrange(timestamp)
  night_baseline <- quantile(cv_data$mhos, 0.25, na.rm = TRUE)
  
  events <- list()
  event_start <- NULL
  event_start_idx <- NULL
  
  for(i in 1:nrow(cv_data)) {
    if(cv_data$in_event[i] && (i == 1 || !cv_data$in_event[i-1])) {
      event_start <- cv_data$timestamp[i]
      event_start_idx <- i
    }
    
    if(!is.null(event_start) && !cv_data$in_event[i] && i > 1 && cv_data$in_event[i-1]) {
      event_end <- cv_data$timestamp[i-1]
      event_end_idx <- i - 1
      event_data <- cv_data[event_start_idx:event_end_idx, ]
      
      duration_min <- nrow(event_data)
      auc <- sum(event_data$mhos, na.rm = TRUE)
      peak_mhos <- max(event_data$mhos, na.rm = TRUE)
      mean_mhos <- mean(event_data$mhos, na.rm = TRUE)
      peak_cv <- max(event_data$rolling_cv, na.rm = TRUE)
      peak_idx <- which.max(event_data$mhos)
      time_to_peak <- peak_idx
      recovery_time <- duration_min - peak_idx
      
      if(nrow(event_data) > 1) {
        rises <- diff(event_data$mhos)
        max_rise_rate <- max(rises, na.rm = TRUE)
        falls <- -rises
        max_fall_rate <- max(falls, na.rm = TRUE)
      } else {
        max_rise_rate <- 0
        max_fall_rate <- 0
      }
      
      cv_during_event <- (sd(event_data$mhos, na.rm = TRUE) / mean(event_data$mhos, na.rm = TRUE)) * 100
      baseline_deviation <- peak_mhos - night_baseline
      peak_to_mean_ratio <- peak_mhos / mean_mhos
      asymmetry <- if(recovery_time > 0) time_to_peak / recovery_time else NA
      initial_mhos <- event_data$mhos[1]
      final_mhos <- event_data$mhos[nrow(event_data)]
      total_rise <- peak_mhos - initial_mhos
      total_fall <- peak_mhos - final_mhos
      
      events[[length(events) + 1]] <- data.frame(
        participant = cv_data$participant[i], visit = cv_data$visit[i], night = cv_data$night[i],
        event_start = event_start, event_end = event_end, duration_min = duration_min, auc = auc,
        peak_mhos = peak_mhos, mean_mhos = mean_mhos, peak_cv = peak_cv,
        time_to_peak = time_to_peak, recovery_time = recovery_time,
        max_rise_rate = max_rise_rate, max_fall_rate = max_fall_rate,
        cv_during_event = cv_during_event, baseline_deviation = baseline_deviation,
        peak_to_mean_ratio = peak_to_mean_ratio, asymmetry = asymmetry,
        initial_mhos = initial_mhos, final_mhos = final_mhos,
        total_rise = total_rise, total_fall = total_fall
      )
      
      event_start <- NULL
      event_start_idx <- NULL
    }
    
    if(i == nrow(cv_data) && !is.null(event_start) && cv_data$in_event[i]) {
      event_end <- cv_data$timestamp[i]
      event_end_idx <- i
      event_data <- cv_data[event_start_idx:event_end_idx, ]
      
      duration_min <- nrow(event_data)
      auc <- sum(event_data$mhos, na.rm = TRUE)
      peak_mhos <- max(event_data$mhos, na.rm = TRUE)
      mean_mhos <- mean(event_data$mhos, na.rm = TRUE)
      peak_cv <- max(event_data$rolling_cv, na.rm = TRUE)
      peak_idx <- which.max(event_data$mhos)
      time_to_peak <- peak_idx
      recovery_time <- duration_min - peak_idx
      
      if(nrow(event_data) > 1) {
        rises <- diff(event_data$mhos)
        max_rise_rate <- max(rises, na.rm = TRUE)
        falls <- -rises
        max_fall_rate <- max(falls, na.rm = TRUE)
      } else {
        max_rise_rate <- 0
        max_fall_rate <- 0
      }
      
      cv_during_event <- (sd(event_data$mhos, na.rm = TRUE) / mean(event_data$mhos, na.rm = TRUE)) * 100
      baseline_deviation <- peak_mhos - night_baseline
      peak_to_mean_ratio <- peak_mhos / mean_mhos
      asymmetry <- if(recovery_time > 0) time_to_peak / recovery_time else NA
      initial_mhos <- event_data$mhos[1]
      final_mhos <- event_data$mhos[nrow(event_data)]
      total_rise <- peak_mhos - initial_mhos
      total_fall <- peak_mhos - final_mhos
      
      events[[length(events) + 1]] <- data.frame(
        participant = cv_data$participant[i], visit = cv_data$visit[i], night = cv_data$night[i],
        event_start = event_start, event_end = event_end, duration_min = duration_min, auc = auc,
        peak_mhos = peak_mhos, mean_mhos = mean_mhos, peak_cv = peak_cv,
        time_to_peak = time_to_peak, recovery_time = recovery_time,
        max_rise_rate = max_rise_rate, max_fall_rate = max_fall_rate,
        cv_during_event = cv_during_event, baseline_deviation = baseline_deviation,
        peak_to_mean_ratio = peak_to_mean_ratio, asymmetry = asymmetry,
        initial_mhos = initial_mhos, final_mhos = final_mhos,
        total_rise = total_rise, total_fall = total_fall
      )
    }
  }
  
  if(length(events) == 0) return(data.frame())
  return(bind_rows(events))
}

# ===== PROCESS ALL PARTICIPANTS =====
cat("Extracting events with enhanced metrics...\n\n")

all_events <- data.frame()

pv_list <- rolling_cv_trimmed %>%
  distinct(pv_label, participant, visit) %>%
  arrange(pv_label)

for(i in 1:nrow(pv_list)) {
  pv <- pv_list$pv_label[i]
  plot_data <- rolling_cv_trimmed %>% filter(pv_label == pv) %>% mutate(night_label = paste("Night", night))
  nights <- unique(plot_data$night)
  
  for(n in nights) {
    night_data <- plot_data %>% filter(night == n)
    processed <- identify_events_dual_5min(night_data, cv_start = 10, cv_end = 10, eda_percentile = 75)
    events <- extract_events_with_metrics(processed)
    if(nrow(events) > 0) all_events <- bind_rows(all_events, events)
  }
  
  if(i %% 10 == 0) cat("Processed", i, "of", nrow(pv_list), "participants\n")
}

cat("\n=== Events Dataset Summary ===\n")
cat("Total events:", nrow(all_events), "\n")
cat("Participants:", n_distinct(all_events$participant), "\n\n")

# Save
events_path <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/events_dataset_full_metrics.csv"
write.csv(all_events, events_path, row.names = FALSE)
cat("Saved to:", events_path, "\n\n")

# ===== VISUALIZATIONS =====
if(!require(ggdist)) install.packages("ggdist")
library(ggdist)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"
dir.create(viz_folder, showWarnings = FALSE)

# Raincloud plots for key metrics
p1 <- ggplot(all_events, aes(x = 1, y = duration_min)) +
  stat_halfeye(adjust = 1, width = 0.6, .width = 0, fill = "#5A7A9B", alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "#E63946", alpha = 0.5) +
  geom_jitter(width = 0.05, alpha = 0.3, size = 0.8, color = "gray30") +
  coord_flip() +
  labs(title = "Event Duration Distribution", y = "Duration (minutes)", x = "") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title = element_text(face = "bold", size = 14), panel.grid.major.y = element_blank())
ggsave(paste0(viz_folder, "/duration_raincloud.png"), p1, width = 8, height = 4, dpi = 300)

p2 <- ggplot(all_events, aes(x = 1, y = max_rise_rate)) +
  stat_halfeye(adjust = 1, width = 0.6, .width = 0, fill = "#5A7A9B", alpha = 0.7) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "#E63946", alpha = 0.5) +
  geom_jitter(width = 0.05, alpha = 0.3, size = 0.8, color = "gray30") +
  coord_flip() +
  labs(title = "Max Rise Rate Distribution", y = "Max Rise Rate (μmhos/min)", x = "") +
  theme_minimal() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        plot.title = element_text(face = "bold", size = 14), panel.grid.major.y = element_blank())
ggsave(paste0(viz_folder, "/max_rise_rate_raincloud.png"), p2, width = 8, height = 4, dpi = 300)

# Scatterplots for clustering
p3 <- ggplot(all_events, aes(x = duration_min, y = max_rise_rate)) +
  geom_point(alpha = 0.5, size = 2, color = "#5A7A9B") +
  labs(title = "Duration vs Max Rise Rate", subtitle = "Looking for HF (short/sharp) vs NS (long/gradual) clusters",
       x = "Duration (minutes)", y = "Max Rise Rate (μmhos/min)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))
ggsave(paste0(viz_folder, "/duration_vs_rise_rate.png"), p3, width = 8, height = 6, dpi = 300)

p4 <- ggplot(all_events, aes(x = time_to_peak, y = recovery_time)) +
  geom_point(alpha = 0.5, size = 2, color = "#5A7A9B") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#E63946") +
  labs(title = "Time to Peak vs Recovery Time", subtitle = "Dashed line = symmetric events",
       x = "Time to Peak (min)", y = "Recovery Time (min)") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))
ggsave(paste0(viz_folder, "/time_to_peak_vs_recovery.png"), p4, width = 8, height = 6, dpi = 300)

p5 <- ggplot(all_events, aes(x = duration_min, y = auc)) +
  geom_point(alpha = 0.5, size = 2, color = "#5A7A9B") +
  geom_smooth(method = "lm", se = TRUE, color = "#E63946") +
  labs(title = "Duration vs Area Under Curve", x = "Duration (min)", y = "AUC") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14))
ggsave(paste0(viz_folder, "/duration_vs_auc.png"), p5, width = 8, height = 6, dpi = 300)

# Multi-panel comparison
p6 <- all_events %>%
  select(duration_min, max_rise_rate, time_to_peak, peak_to_mean_ratio, asymmetry) %>%
  pivot_longer(-duration_min, names_to = "metric", values_to = "value") %>%
  filter(!is.na(value) & !is.infinite(value)) %>%
  ggplot(aes(x = duration_min, y = value)) +
  geom_point(alpha = 0.4, size = 1.5, color = "#5A7A9B") +
  facet_wrap(~ metric, scales = "free_y", ncol = 2) +
  labs(title = "Event Characteristics by Duration", x = "Duration (minutes)", y = "Value") +
  theme_minimal() +
  theme(plot.title = element_text(face = "bold", size = 14), strip.text = element_text(face = "bold"))
ggsave(paste0(viz_folder, "/all_metrics_by_duration.png"), p6, width = 10, height = 8, dpi = 300)

cat("\n=== Visualizations Complete ===\n")
cat("Saved to:", viz_folder, "\n")


################################################################################
# VARIABLE SELECTION — correlation check before deciding on final variables
################################################################################

####################checking correlation##########################
library(tidyverse)
library(corrplot)

all_events <- read.csv("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/events_dataset_full_metrics.csv",
                       stringsAsFactors = FALSE)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"

# Calculate correlations
metrics_check <- all_events %>%
  select(duration_min, auc, time_to_peak, recovery_time, max_rise_rate, 
         peak_to_mean_ratio, cv_during_event) %>%
  na.omit()

cor_matrix <- cor(metrics_check)

cat("========================================\n")
cat("   CORRELATION MATRIX\n")
cat("========================================\n\n")

print(round(cor_matrix, 2))

# Visualize
png(paste0(viz_folder, "/correlation_heatmap.png"), width = 800, height = 800)
corrplot(cor_matrix, method = "color", type = "upper", 
         addCoef.col = "black", number.cex = 0.8,
         tl.col = "black", tl.srt = 45,
         col = colorRampPalette(c("#E63946", "white", "#5A7A9B"))(200),
         title = "Metric Correlations (Multicollinearity Check)",
         mar = c(0,0,2,0))
dev.off()

# Check specific correlations
cat("\n========================================\n")
cat("   KEY CORRELATIONS\n")
cat("========================================\n\n")

cat("Duration vs AUC:", round(cor(metrics_check$duration_min, metrics_check$auc), 3), "\n")
cat("Duration vs Time to Peak:", round(cor(metrics_check$duration_min, metrics_check$time_to_peak), 3), "\n")
cat("Duration vs Recovery:", round(cor(metrics_check$duration_min, metrics_check$recovery_time), 3), "\n")
cat("Time to Peak vs Recovery:", round(cor(metrics_check$time_to_peak, metrics_check$recovery_time), 3), "\n\n")

# Check if duration ≈ time_to_peak + recovery
all_events_check <- all_events %>%
  mutate(
    sum_components = time_to_peak + recovery_time,
    difference = duration_min - sum_components
  ) %>%
  filter(!is.na(sum_components))

cat("DURATION vs (TIME_TO_PEAK + RECOVERY):\n")
cat("Correlation:", round(cor(all_events_check$duration_min, all_events_check$sum_components, use="complete.obs"), 3), "\n")
cat("Mean difference:", round(mean(all_events_check$difference, na.rm=TRUE), 2), "minutes\n")
cat("SD of difference:", round(sd(all_events_check$difference, na.rm=TRUE), 2), "minutes\n\n")




################################################################################
# THE ACTUAL CLUSTERING — fit, drop the bad cluster, refit, characterize
################################################################################


# initial GMM fit on the 3 final variables (duration, AUC, time to peak), k=1-4
################re running with 3 variables ###################################
library(mclust)
library(cluster)
library(moments)
library(factoextra)
library(ggplot2)
library(dplyr)

# ===== SETUP =====

all_events <- read.csv("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/events_dataset_full_metrics.csv",
                       stringsAsFactors = FALSE)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"

metrics <- c("duration_min", "auc", "time_to_peak")

cluster_data <- na.omit(all_events[, metrics])

cat("N events:", nrow(cluster_data), "\n")
cat("Variables: duration_min, auc, time_to_peak\n")
cat("Recovery time dropped (mathematically dependent: duration = time_to_peak + recovery)\n\n")

# ===== LOG TRANSFORM =====

cluster_data_log <- data.frame(
  log_duration     = log(cluster_data$duration_min + 0.1),
  log_auc          = log(cluster_data$auc          + 0.1),
  log_time_to_peak = log(cluster_data$time_to_peak + 0.1)
)

cluster_data_log_scaled <- scale(cluster_data_log)

# ============================================================
# step 2: GMM fit, k = 1 to 4
# ============================================================

cat("\n\n========================================\n")
cat("   GMM: k = 1 to 4\n")
cat("========================================\n\n")

gmm_models  <- list()
results_table <- data.frame()

for (k in 1:4) {
  set.seed(123)
  gmm_k <- Mclust(cluster_data_log_scaled, G = k)
  gmm_models[[k]] <- gmm_k
  
  # Cluster sizes
  sizes <- table(factor(gmm_k$classification, levels = 1:4))
  
  results_table <- rbind(results_table, data.frame(
    k              = k,
    Total_N        = nrow(cluster_data_log_scaled),
    Model_Type     = gmm_k$modelName,
    Parameters     = gmm_k$df,
    Log_Likelihood = round(gmm_k$loglik, 2),
    BIC            = round(gmm_k$bic,    2),
    AIC            = round(2 * gmm_k$df - 2 * gmm_k$loglik, 2),
    Cluster_1_n    = as.integer(sizes[1]),
    Cluster_2_n    = ifelse(k >= 2, as.integer(sizes[2]), NA),
    Cluster_3_n    = ifelse(k >= 3, as.integer(sizes[3]), NA),
    Cluster_4_n    = ifelse(k >= 4, as.integer(sizes[4]), NA)
  ))
  
  cat("k=", k, ": BIC =", round(gmm_k$bic, 2),
      "| Model:", gmm_k$modelName,
      "| Sizes:", paste(as.integer(sizes[1:k]), collapse = "/"), "\n")
}

write.csv(results_table,
          paste0(viz_folder, "/GMM3var_Model_Fit_k1to4.csv"),
          row.names = FALSE)

# ---- BIC differences ----

bic_diffs <- data.frame(
  Comparison     = c("k=2 vs k=1", "k=3 vs k=2", "k=4 vs k=3"),
  BIC_Difference = c(
    results_table$BIC[2] - results_table$BIC[1],
    results_table$BIC[3] - results_table$BIC[2],
    results_table$BIC[4] - results_table$BIC[3]
  )
)

bic_diffs$Strength <- ifelse(abs(bic_diffs$BIC_Difference) < 2,  "Weak",
                             ifelse(abs(bic_diffs$BIC_Difference) < 6,  "Positive",
                                    ifelse(abs(bic_diffs$BIC_Difference) < 10, "Strong", "Very Strong")))

cat("\n--- BIC differences ---\n")
print(bic_diffs, row.names = FALSE)
write.csv(bic_diffs,
          paste0(viz_folder, "/GMM3var_BIC_Differences.csv"),
          row.names = FALSE)

# ---- Silhouette widths for k = 2, 3, 4 ----

sil_table <- data.frame()

for (k in 2:4) {
  sil <- silhouette(gmm_models[[k]]$classification, dist(cluster_data_log_scaled))
  sil_table <- rbind(sil_table, data.frame(
    k              = k,
    Avg_Silhouette = round(mean(sil[, 3]), 3)
  ))
}

cat("\n--- Average silhouette width ---\n")
print(sil_table, row.names = FALSE)
write.csv(sil_table,
          paste0(viz_folder, "/GMM3var_Silhouette_Width.csv"),
          row.names = FALSE)

# ============================================================
# step 3: cluster characterization (k = 4)
# ============================================================

cat("\n\n========================================\n")
cat("   CLUSTER CHARACTERIZATION: k = 4\n")
cat("========================================\n\n")

gmm_k4 <- gmm_models[[4]]

cluster_data$gmm_cluster <- gmm_k4$classification

# Order clusters by mean duration so labels are meaningful
cluster_order <- cluster_data %>%
  group_by(gmm_cluster) %>%
  summarise(mean_dur = mean(duration_min), .groups = "drop") %>%
  arrange(mean_dur) %>%
  mutate(rank = row_number(),
         cluster_label = c("Cluster 1 (shortest)",
                           "Cluster 2",
                           "Cluster 3",
                           "Cluster 4 (longest)"))

cluster_data <- cluster_data %>%
  left_join(cluster_order[, c("gmm_cluster", "cluster_label")], by = "gmm_cluster")

cat("Cluster sizes (ordered by duration):\n")
print(as.data.frame(cluster_order), row.names = FALSE)

# ---- Descriptive stats per cluster ----

test_vars <- c("duration_min", "auc", "time_to_peak")
test_labs <- c("Duration (min)", "AUC (mhos-min)", "Time to Peak (min)")

desc_stats <- data.frame()

for (i in seq_along(test_vars)) {
  var <- test_vars[i]
  lab <- test_labs[i]
  
  if (!var %in% colnames(cluster_data)) next
  
  for (cl in unique(cluster_data$cluster_label)) {
    vals <- cluster_data[[var]][cluster_data$cluster_label == cl]
    
    desc_stats <- rbind(desc_stats, data.frame(
      Variable      = lab,
      Cluster       = cl,
      N             = length(vals),
      Mean          = round(mean(vals, na.rm = TRUE), 2),
      SD            = round(sd(vals,   na.rm = TRUE), 2),
      Median        = round(median(vals, na.rm = TRUE), 2)
    ))
  }
}

cat("\n--- Descriptive statistics per cluster ---\n")
print(desc_stats, row.names = FALSE)
write.csv(desc_stats,
          paste0(viz_folder, "/GMM3var_k4_Descriptive_Stats.csv"),
          row.names = FALSE)

# ---- Save events with cluster labels ----

all_events$gmm_cluster    <- NA
all_events$cluster_label  <- NA

valid_rows <- which(complete.cases(all_events[, metrics]))
all_events$gmm_cluster[valid_rows]  <- cluster_data$gmm_cluster
all_events$cluster_label[valid_rows] <- cluster_data$cluster_label

write.csv(all_events,
          paste0(viz_folder, "/events_dataset_GMM3var_k4_clusters.csv"),
          row.names = FALSE)

# ============================================================
# step 4: PCA visualization (k = 4)
# ============================================================

cat("\n\nCreating PCA visualization...\n")

pca_result <- prcomp(cluster_data_log_scaled, scale. = FALSE)
var_exp    <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

pca_df <- data.frame(
  PC1     = pca_result$x[, 1],
  PC2     = pca_result$x[, 2],
  Cluster = cluster_data$cluster_label
)

colors <- c(
  "Cluster 1 (shortest)" = "#1D9E75",
  "Cluster 2"            = "#378ADD",
  "Cluster 3"            = "#EF9F27",
  "Cluster 4 (longest)"  = "#D85A30"
)

p_pca <- ggplot(pca_df, aes(x = PC1, y = PC2, color = Cluster, fill = Cluster)) +
  geom_point(alpha = 0.45, size = 1.8) +
  stat_ellipse(level = 0.95, geom = "polygon", alpha = 0.10, linewidth = 1) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(
    title    = "GMM k=4 Clusters — PCA Visualization (3-variable solution)",
    subtitle = paste0("PC1: ", var_exp[1], "%, PC2: ", var_exp[2],
                      "% | Clusters ordered by mean duration"),
    x        = paste0("PC1 (", var_exp[1], "% variance)"),
    y        = paste0("PC2 (", var_exp[2], "% variance)"),
    color    = NULL, fill = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", hjust = 0.5),
    plot.subtitle    = element_text(hjust = 0.5, size = 10)
  )

ggsave(paste0(viz_folder, "/GMM3var_k4_PCA.png"),
       p_pca, width = 9, height = 7, dpi = 300)

cat("Saved: GMM3var_k4_PCA.png\n")

# ============================================================
# step 5: summary
# ============================================================

cat("\n\n========================================\n")
cat("   QUICK RESULTS PREVIEW\n")
cat("========================================\n\n")

cat("BIC (higher = better):\n")
for (i in 1:4) cat("  k=", i, ": ", results_table$BIC[i], "\n", sep="")

cat("\nBIC differences:\n")
print(bic_diffs[, c("Comparison", "BIC_Difference", "Strength")], row.names = FALSE)

cat("\nSilhouette width:\n")
for (i in 1:nrow(sil_table)) cat("  k=", sil_table$k[i], ": ", sil_table$Avg_Silhouette[i], "\n", sep="")

cat("\nCluster sizes (k=4, ordered by duration):\n")
for (cl in cluster_order$cluster_label) {
  cat("  ", cl, ": n =", sum(cluster_data$cluster_label == cl), "\n")
}

cat("\n========================================\n")
cat("   FILES SAVED\n")
cat("========================================\n\n")
cat("Assumption checks:\n")
cat("  GMM3var_Assumption_Normality_Original.csv\n")
cat("  GMM3var_Assumption_Normality_Transformed.csv\n")
cat("  GMM3var_Correlation_Matrix.csv\n")
cat("  GMM3var_Assumption_Multicollinearity.csv\n")
cat("  GMM3var_Assumption_Outliers.csv\n\n")
cat("GMM results:\n")
cat("  GMM3var_Model_Fit_k1to4.csv\n")
cat("  GMM3var_BIC_Differences.csv\n")
cat("  GMM3var_Silhouette_Width.csv\n")
cat("  GMM3var_k4_Descriptive_Stats.csv\n")
cat("  GMM3var_k4_Pairwise_Tests.csv\n")
cat("  events_dataset_GMM3var_k4_clusters.csv\n\n")
cat("Visualization:\n")
cat("  GMM3var_k4_PCA.png\n")

# fit both k=2 and k=4 on the same data, label everything, save it out -

library(mclust)
library(cluster)
library(ggplot2)
library(patchwork)
library(dplyr)

# ===== SETUP =====

all_events <- read.csv("C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/events_dataset_full_metrics.csv",
                       stringsAsFactors = FALSE)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"

metrics      <- c("duration_min", "auc", "time_to_peak")
cluster_data <- na.omit(all_events[, metrics])

cluster_data_log <- data.frame(
  log_duration     = log(cluster_data$duration_min + 0.1),
  log_auc          = log(cluster_data$auc          + 0.1),
  log_time_to_peak = log(cluster_data$time_to_peak + 0.1)
)

cluster_data_log_scaled <- scale(cluster_data_log)

# ===== FIT k=2 AND k=4 =====

set.seed(123)
gmm_k2 <- Mclust(cluster_data_log_scaled, G = 2)

set.seed(123)
gmm_k4 <- Mclust(cluster_data_log_scaled, G = 4)

cat("k=2: BIC =", round(gmm_k2$bic, 2), "\n")
cat("k=4: BIC =", round(gmm_k4$bic, 2), "\n\n")

# ===== PCA =====

pca_result <- prcomp(cluster_data_log_scaled, scale. = FALSE)
var_exp    <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

pc_scores <- data.frame(
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2]
)

# ---- Order k=2 clusters by duration ----

cluster_data$k2 <- gmm_k2$classification

k2_order <- cluster_data %>%
  group_by(k2) %>%
  summarise(mean_dur = mean(duration_min), .groups = "drop") %>%
  arrange(mean_dur) %>%
  mutate(k2_label = c("Short-duration\n(n=XXX)", "Long-duration\n(n=XXX)"))

# Replace XXX with actual counts
k2_order$k2_label <- paste0(
  c("Short-duration", "Long-duration"),
  "\n(n=", table(cluster_data$k2)[k2_order$k2], ")"
)

cluster_data <- cluster_data %>%
  left_join(k2_order[, c("k2", "k2_label")], by = "k2")

# ---- Order k=4 clusters by duration ----

cluster_data$k4 <- gmm_k4$classification

k4_order <- cluster_data %>%
  group_by(k4) %>%
  summarise(mean_dur = mean(duration_min), .groups = "drop") %>%
  arrange(mean_dur) %>%
  mutate(k4_label = paste0(
    c("Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4"),
    "\n(n=", table(cluster_data$k4)[k4[order(mean_dur)]], ")"
  ))

# Safer approach for k4 labels
k4_sizes <- as.data.frame(table(cluster_data$k4))
colnames(k4_sizes) <- c("k4", "n")
k4_sizes$k4 <- as.integer(as.character(k4_sizes$k4))

k4_order <- cluster_data %>%
  group_by(k4) %>%
  summarise(mean_dur = mean(duration_min), .groups = "drop") %>%
  arrange(mean_dur) %>%
  left_join(k4_sizes, by = "k4") %>%
  mutate(rank = row_number(),
         k4_label = paste0("Cluster ", rank, "\n(n=", n, ")"))

cluster_data <- cluster_data %>%
  left_join(k4_order[, c("k4", "k4_label")], by = "k4")

# ---- Add labels to PC scores ----

pc_scores$k2_label <- cluster_data$k2_label
pc_scores$k4_label <- cluster_data$k4_label

# ===== SIDE-BY-SIDE PCA PLOT =====

colors_k2 <- c(
  grep("Short", unique(pc_scores$k2_label), value = TRUE) |> setNames("#1D9E75", .),
  grep("Long",  unique(pc_scores$k2_label), value = TRUE) |> setNames("#D85A30", .)
)

# Simpler color assignment
k2_unique <- unique(pc_scores$k2_label)
k2_unique_sorted <- k2_unique[order(nchar(k2_unique))]  # short label first

colors_k2 <- c("#1D9E75", "#D85A30")
names(colors_k2) <- k2_unique[order(
  sapply(k2_unique, function(x) {
    cluster_data$mean_dur[match(
      cluster_data$k2_label, x
    )] |> mean(na.rm = TRUE)
  })
)]

# Cleaner approach
k2_labels_ordered <- k2_order$k2_label
k4_labels_ordered <- k4_order$k4_label

colors_k2 <- setNames(c("#1D9E75", "#D85A30"), k2_labels_ordered)
colors_k4 <- setNames(c("#1D9E75", "#378ADD", "#EF9F27", "#D85A30"), k4_labels_ordered)

# k=2 plot
p_k2 <- ggplot(pc_scores, aes(x = PC1, y = PC2,
                              color = k2_label, fill = k2_label)) +
  geom_point(alpha = 0.45, size = 1.5) +
  stat_ellipse(level = 0.95, geom = "polygon",
               alpha = 0.12, linewidth = 0.9) +
  scale_color_manual(values = colors_k2, name = NULL) +
  scale_fill_manual(values  = colors_k2, name = NULL) +
  labs(
    title    = "k = 2",
    subtitle = paste0("BIC = ", round(gmm_k2$bic, 0),
                      " | Silhouette = 0.113"),
    x        = paste0("PC1 (", var_exp[1], "%)"),
    y        = paste0("PC2 (", var_exp[2], "%)")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle    = element_text(hjust = 0.5, size = 9,
                                    color = "gray40")
  )

# k=4 plot
p_k4 <- ggplot(pc_scores, aes(x = PC1, y = PC2,
                              color = k4_label, fill = k4_label)) +
  geom_point(alpha = 0.45, size = 1.5) +
  stat_ellipse(level = 0.95, geom = "polygon",
               alpha = 0.10, linewidth = 0.9) +
  scale_color_manual(values = colors_k4, name = NULL) +
  scale_fill_manual(values  = colors_k4, name = NULL) +
  labs(
    title    = "k = 4",
    subtitle = paste0("BIC = ", round(gmm_k4$bic, 0),
                      " | Silhouette = 0.140"),
    x        = paste0("PC1 (", var_exp[1], "%)"),
    y        = paste0("PC2 (", var_exp[2], "%)")
  ) +
  theme_minimal(base_size = 11) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 9),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 14),
    plot.subtitle    = element_text(hjust = 0.5, size = 9,
                                    color = "gray40")
  )

# Combined
p_combined <- p_k2 + p_k4 +
  plot_annotation(
    title    = "Gaussian Mixture Model Cluster Solutions",
    subtitle = paste0("PCA captures ", var_exp[1] + var_exp[2],
                      "% of variance | 95% confidence ellipses"),
    theme = theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 15),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40")
    )
  )

ggsave(paste0(viz_folder, "/GMM3var_k2_vs_k4_PCA.png"),
       p_combined, width = 13, height = 6, dpi = 300)

cat("Saved: GMM3var_k2_vs_k4_PCA.png\n\n")

# --- cluster characterization: k=2 ---

cat("========================================\n")
cat("   CHARACTERIZATION: k = 2\n")
cat("========================================\n\n")

valid_rows <- which(complete.cases(all_events[, metrics]))

test_vars <- c("duration_min", "auc", "time_to_peak")
test_labs <- c("Duration (min)", "AUC (mhos-min)", "Time to Peak (min)")

desc_k2 <- data.frame()

for (i in seq_along(test_vars)) {
  var <- test_vars[i]
  lab <- test_labs[i]
  if (!var %in% colnames(cluster_data)) next
  
  for (cl in k2_labels_ordered) {
    vals <- cluster_data[[var]][cluster_data$k2_label == cl]
    desc_k2 <- rbind(desc_k2, data.frame(
      Variable = lab,
      Cluster  = cl,
      N        = length(vals),
      Mean     = round(mean(vals, na.rm = TRUE), 2),
      SD       = round(sd(vals,   na.rm = TRUE), 2),
      Median   = round(median(vals, na.rm = TRUE), 2)
    ))
  }
}

print(desc_k2, row.names = FALSE)
write.csv(desc_k2,
          paste0(viz_folder, "/GMM3var_k2_Descriptive_Stats.csv"),
          row.names = FALSE)

# --- cluster characterization: k=4 ---

cat("\n========================================\n")
cat("   CHARACTERIZATION: k = 4\n")
cat("========================================\n\n")

desc_k4 <- data.frame()

for (i in seq_along(test_vars)) {
  var <- test_vars[i]
  lab <- test_labs[i]
  if (!var %in% colnames(cluster_data)) next
  
  for (cl in k4_labels_ordered) {
    vals <- cluster_data[[var]][cluster_data$k4_label == cl]
    desc_k4 <- rbind(desc_k4, data.frame(
      Variable = lab,
      Cluster  = cl,
      N        = length(vals),
      Mean     = round(mean(vals, na.rm = TRUE), 2),
      SD       = round(sd(vals,   na.rm = TRUE), 2),
      Median   = round(median(vals, na.rm = TRUE), 2)
    ))
  }
}

print(desc_k4, row.names = FALSE)
write.csv(desc_k4,
          paste0(viz_folder, "/GMM3var_k4_Descriptive_Stats.csv"),
          row.names = FALSE)

# ===== SAVE EVENTS WITH BOTH CLUSTER LABELS =====

all_events$gmm_k2_cluster <- NA
all_events$gmm_k2_label   <- NA
all_events$gmm_k4_cluster <- NA
all_events$gmm_k4_label   <- NA

all_events$gmm_k2_cluster[valid_rows] <- cluster_data$k2
all_events$gmm_k2_label[valid_rows]   <- cluster_data$k2_label
all_events$gmm_k4_cluster[valid_rows] <- cluster_data$k4
all_events$gmm_k4_label[valid_rows]   <- cluster_data$k4_label

write.csv(all_events,
          paste0(viz_folder, "/events_dataset_GMM3var_k2_k4_clusters.csv"),
          row.names = FALSE)

cat("Saved: events_dataset_GMM3var_k2_k4_clusters.csv\n")

# cluster 3 (n=63) small and scattered all over the PCA space instead of forming a real cluster,
# so it got cut and everything got refit on the remaining 1027 events.
# this is the step described in the manuscript's Methods section.
# removing cluster 3 events and re-fitting
library(mclust)
library(cluster)
library(ggplot2)
library(ggdist)
library(patchwork)
library(dplyr)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"

# ===== LOAD AND FILTER =====

events_labeled <- read.csv(paste0(viz_folder, "/events_dataset_GMM3var_k2_k4_clusters.csv"),
                           stringsAsFactors = FALSE)

cat("Original N:", nrow(events_labeled), "\n")

# Identify and remove Cluster 3 events
# Cluster 3 is the 3rd shortest by duration in the k=4 solution
cluster3_label <- events_labeled %>%
  filter(!is.na(gmm_k4_label)) %>%
  group_by(gmm_k4_label) %>%
  summarise(mean_dur = mean(duration_min, na.rm = TRUE), .groups = "drop") %>%
  arrange(mean_dur) %>%
  slice(3) %>%
  pull(gmm_k4_label)

cat("Removing cluster label:", cluster3_label, "\n")
cat("N in Cluster 3:", sum(events_labeled$gmm_k4_label == cluster3_label, na.rm = TRUE), "\n")

events_filtered <- events_labeled %>%
  filter(is.na(gmm_k4_label) | gmm_k4_label != cluster3_label)

cat("N after removal:", nrow(events_filtered), "\n\n")

# Save filtered dataset
write.csv(events_filtered,
          paste0(viz_folder, "/events_dataset_cluster3_removed.csv"),
          row.names = FALSE)

cat("Saved: events_dataset_cluster3_removed.csv\n\n")

# ===== PREPARE DATA FOR RE-CLUSTERING =====

metrics      <- c("duration_min", "auc", "time_to_peak")
cluster_data <- events_filtered %>%
  filter(!is.na(duration_min), !is.na(auc), !is.na(time_to_peak))

cat("N events for re-clustering:", nrow(cluster_data), "\n\n")

cluster_data_log <- data.frame(
  log_duration     = log(cluster_data$duration_min + 0.1),
  log_auc          = log(cluster_data$auc          + 0.1),
  log_time_to_peak = log(cluster_data$time_to_peak + 0.1)
)

cluster_data_log_scaled <- scale(cluster_data_log)

# ===== RUN GMM k=1 TO k=4 =====

cat("========================================\n")
cat("   GMM: k = 1 to 4 (Cluster 3 removed)\n")
cat("========================================\n\n")

gmm_models    <- list()
results_table <- data.frame()

for (k in 1:4) {
  set.seed(123)
  gmm_k <- Mclust(cluster_data_log_scaled, G = k)
  gmm_models[[k]] <- gmm_k
  
  sizes <- table(factor(gmm_k$classification, levels = 1:4))
  
  results_table <- rbind(results_table, data.frame(
    k              = k,
    Total_N        = nrow(cluster_data_log_scaled),
    Model_Type     = gmm_k$modelName,
    Parameters     = gmm_k$df,
    Log_Likelihood = round(gmm_k$loglik, 2),
    BIC            = round(gmm_k$bic,    2),
    AIC            = round(2 * gmm_k$df - 2 * gmm_k$loglik, 2),
    Cluster_1_n    = as.integer(sizes[1]),
    Cluster_2_n    = ifelse(k >= 2, as.integer(sizes[2]), NA),
    Cluster_3_n    = ifelse(k >= 3, as.integer(sizes[3]), NA),
    Cluster_4_n    = ifelse(k >= 4, as.integer(sizes[4]), NA)
  ))
  
  cat("k=", k, ": BIC =", round(gmm_k$bic, 2),
      "| Model:", gmm_k$modelName,
      "| Sizes:", paste(as.integer(sizes[1:k]), collapse = "/"), "\n")
}

write.csv(results_table,
          paste0(viz_folder, "/GMM3var_noC3_Model_Fit_k1to4.csv"),
          row.names = FALSE)

# ===== BIC DIFFERENCES =====

bic_diffs <- data.frame(
  Comparison     = c("k=2 vs k=1", "k=3 vs k=2", "k=4 vs k=3"),
  BIC_Difference = c(
    results_table$BIC[2] - results_table$BIC[1],
    results_table$BIC[3] - results_table$BIC[2],
    results_table$BIC[4] - results_table$BIC[3]
  )
)

bic_diffs$Strength <- ifelse(abs(bic_diffs$BIC_Difference) < 2,  "Weak",
                             ifelse(abs(bic_diffs$BIC_Difference) < 6,  "Positive",
                                    ifelse(abs(bic_diffs$BIC_Difference) < 10, "Strong",
                                           "Very Strong")))

cat("\n--- BIC differences ---\n")
print(bic_diffs, row.names = FALSE)

write.csv(bic_diffs,
          paste0(viz_folder, "/GMM3var_noC3_BIC_Differences.csv"),
          row.names = FALSE)

# ===== SILHOUETTE WIDTHS =====

sil_table <- data.frame()

for (k in 2:4) {
  sil <- silhouette(gmm_models[[k]]$classification,
                    dist(cluster_data_log_scaled))
  sil_table <- rbind(sil_table, data.frame(
    k              = k,
    Avg_Silhouette = round(mean(sil[, 3]), 3)
  ))
}

cat("\n--- Silhouette widths ---\n")
print(sil_table, row.names = FALSE)

write.csv(sil_table,
          paste0(viz_folder, "/GMM3var_noC3_Silhouette_Width.csv"),
          row.names = FALSE)

# ===== CHARACTERIZE k=4 SOLUTION =====

cat("\n\n========================================\n")
cat("   CLUSTER CHARACTERIZATION: k = 4\n")
cat("========================================\n\n")

gmm_k4 <- gmm_models[[4]]
cluster_data$k4 <- gmm_k4$classification

k4_order <- cluster_data %>%
  group_by(k4) %>%
  summarise(mean_dur = mean(duration_min), .groups = "drop") %>%
  arrange(mean_dur) %>%
  mutate(rank = row_number(),
         k4_label = paste0("Cluster ", rank))

cluster_data <- cluster_data %>%
  left_join(k4_order[, c("k4", "k4_label")], by = "k4")

k4_labels_ordered <- paste0("Cluster ", 1:4)

cat("Cluster sizes (ordered by duration):\n")
for (cl in k4_labels_ordered) {
  n <- sum(cluster_data$k4_label == cl)
  cat(" ", cl, ": n =", n, "(", round(100*n/nrow(cluster_data), 1), "%)\n")
}

# ===== DESCRIPTIVE STATS =====

test_vars <- c("duration_min", "auc", "time_to_peak")
test_labs <- c("Duration (min)", "AUC (mhos-min)", "Time to Peak (min)")

cat("\n--- Descriptive stats per cluster ---\n\n")

for (i in seq_along(test_vars)) {
  var <- test_vars[i]
  lab <- test_labs[i]
  if (!var %in% colnames(cluster_data)) next
  cat(lab, ":\n")
  for (cl in k4_labels_ordered) {
    vals <- cluster_data[[var]][cluster_data$k4_label == cl]
    cat("  ", cl, ": M =", round(mean(vals, na.rm=TRUE), 2),
        ", SD =", round(sd(vals, na.rm=TRUE), 2), "\n")
  }
}

# Save descriptive stats
desc_rows <- data.frame()

for (i in seq_along(test_vars)) {
  var <- test_vars[i]
  lab <- test_labs[i]
  if (!var %in% colnames(cluster_data)) next
  
  row_vals <- sapply(k4_labels_ordered, function(cl) {
    vals <- cluster_data[[var]][cluster_data$k4_label == cl]
    paste0(round(mean(vals, na.rm=TRUE), 2),
           " (", round(sd(vals, na.rm=TRUE), 2), ")")
  })
  
  desc_rows <- rbind(desc_rows, c(Variable = lab, row_vals))
}

sizes_row <- c("n (%)", sapply(k4_labels_ordered, function(cl) {
  n <- sum(cluster_data$k4_label == cl)
  paste0(n, " (", round(100 * n / nrow(cluster_data), 1), "%)")
}))

table2 <- rbind(sizes_row, desc_rows)
colnames(table2) <- c("Variable",
                      "Cluster 1 (shortest)",
                      "Cluster 2",
                      "Cluster 3",
                      "Cluster 4 (longest)")

write.csv(table2,
          paste0(viz_folder, "/GMM3var_noC3_k4_Descriptive_Stats.csv"),
          row.names = FALSE)

cat("\nSaved: GMM3var_noC3_k4_Descriptive_Stats.csv\n")

# ===== PCA VISUALIZATION =====

pca_result <- prcomp(cluster_data_log_scaled, scale. = FALSE)
var_exp    <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

colors_k4 <- c("Cluster 1" = "#1D9E75",
               "Cluster 2" = "#378ADD",
               "Cluster 3" = "#EF9F27",
               "Cluster 4" = "#D85A30")

k4_sizes_vec <- sapply(k4_labels_ordered, function(cl)
  sum(cluster_data$k4_label == cl))

k4_display <- paste0(k4_labels_ordered, " (n=", k4_sizes_vec, ")")
names(k4_display) <- k4_labels_ordered

cluster_data$k4_display <- k4_display[cluster_data$k4_label]
colors_display <- setNames(unname(colors_k4), k4_display)

pc_df <- data.frame(
  PC1        = pca_result$x[, 1],
  PC2        = pca_result$x[, 2],
  k4_display = factor(cluster_data$k4_display, levels = k4_display)
)

p_pca <- ggplot(pc_df, aes(x = PC1, y = PC2,
                           color = k4_display,
                           fill  = k4_display)) +
  geom_point(alpha = 0.45, size = 1.8) +
  stat_ellipse(level = 0.95, geom = "polygon",
               alpha = 0.10, linewidth = 0.9) +
  scale_color_manual(values = colors_display, name = NULL) +
  scale_fill_manual(values  = colors_display, name = NULL) +
  labs(
    title    = "GMM k=4 cluster solution (Cluster 3 removed)",
    subtitle = paste0("PC1: ", var_exp[1], "%, PC2: ", var_exp[2],
                      "% (", sum(var_exp), "% total variance) | 95% confidence ellipses"),
    x = paste0("PC1 (", var_exp[1], "%)"),
    y = paste0("PC2 (", var_exp[2], "%)")
  ) +
  theme_minimal(base_size = 12) +
  theme(
    legend.position  = "bottom",
    legend.text      = element_text(size = 10),
    panel.grid.minor = element_blank(),
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "gray50")
  )

ggsave(paste0(viz_folder, "/GMM3var_noC3_k4_PCA.png"),
       p_pca, width = 8, height = 6, dpi = 300)

cat("Saved: GMM3var_noC3_k4_PCA.png\n")

# just the simple stats (sizes, mixing weights, means) - used these to
# build the decision tree in Canva, so no need for the full covariance
# matrices here.
# ===== EXTRACT PARAMETERS FOR DECISION TREE =====

# k=2 model
gmm_k2 <- gmm_models[[2]]

cat("\n=== K2 CLUSTER SIZES ===\n")
print(table(gmm_k2$classification))
cat("Mixing weights (pi):", round(gmm_k2$parameters$pro, 4), "\n")

cat("\n=== K4 CLUSTER SIZES ===\n")
print(table(gmm_k4$classification))  # gmm_k4 already defined in this section
cat("Mixing weights (pi):", round(gmm_k4$parameters$pro, 4), "\n")

cat("\n=== CROSSTAB K2 vs K4 ===\n")
print(table(k2 = gmm_k2$classification, k4 = cluster_data$k4_label))

cat("\n=== K2 MEANS (log scale) ===\n")
print(round(gmm_k2$parameters$mean, 4))

cat("\n=== K4 MEANS (log scale) ===\n")
print(round(gmm_k4$parameters$mean, 4))

# Note: full covariance matrices (gmm_k2$parameters$variance$sigma,
# gmm_k4$parameters$variance$sigma) are available via mclust's output object
# but are omitted here since only the simple summary stats (n, pi, means)
# were used to build the decision tree diagram.

# assumption checks on the final (post-removal) sample -> Supplemental Table 2
#####assumption checks###########
library(dplyr)
library(moments)
library(mclust)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"

events_filtered <- read.csv(
  paste0(viz_folder, "/events_dataset_cluster3_removed.csv"),
  stringsAsFactors = FALSE
)

metrics      <- c("duration_min", "auc", "time_to_peak")
cluster_data <- events_filtered %>%
  filter(!is.na(duration_min), !is.na(auc), !is.na(time_to_peak))

cluster_data_log <- data.frame(
  log_duration     = log(cluster_data$duration_min + 0.1),
  log_auc          = log(cluster_data$auc          + 0.1),
  log_time_to_peak = log(cluster_data$time_to_peak + 0.1)
)

cluster_data_log_scaled <- scale(cluster_data_log)

# ===== NORMALITY: ORIGINAL =====

norm_orig <- data.frame()

for (var in metrics) {
  sw   <- shapiro.test(cluster_data[[var]])
  norm_orig <- rbind(norm_orig, data.frame(
    Variable  = var,
    Transform = "Original",
    Shapiro_W = round(sw$statistic, 4),
    Shapiro_p = ifelse(sw$p.value < .001, "<.001", round(sw$p.value, 3)),
    Skewness  = round(skewness(cluster_data[[var]]), 3),
    Kurtosis  = round(kurtosis(cluster_data[[var]]), 3)
  ))
}

# ===== NORMALITY: LOG-TRANSFORMED =====

log_vars  <- c("log_duration", "log_auc", "log_time_to_peak")
log_labs  <- c("duration_min", "auc", "time_to_peak")

norm_log <- data.frame()

for (i in seq_along(log_vars)) {
  var  <- log_vars[i]
  sw   <- shapiro.test(cluster_data_log[[var]])
  norm_log <- rbind(norm_log, data.frame(
    Variable  = log_labs[i],
    Transform = "Log",
    Shapiro_W = round(sw$statistic, 4),
    Shapiro_p = ifelse(sw$p.value < .001, "<.001", round(sw$p.value, 3)),
    Skewness  = round(skewness(cluster_data_log[[var]]), 3),
    Kurtosis  = round(kurtosis(cluster_data_log[[var]]), 3)
  ))
}

normality_table <- rbind(norm_orig, norm_log) %>%
  arrange(Variable, Transform)

# ===== CORRELATION MATRIX =====

cor_matrix <- cor(cluster_data[, metrics], use = "complete.obs")
cor_df <- as.data.frame(round(cor_matrix, 3))
cor_df$Variable <- rownames(cor_df)
cor_df <- cor_df %>% select(Variable, everything())
rownames(cor_df) <- NULL

# ===== OUTLIERS (IQR) =====

outlier_table <- data.frame()

for (var in metrics) {
  q1    <- quantile(cluster_data[[var]], .25)
  q3    <- quantile(cluster_data[[var]], .75)
  iqr   <- q3 - q1
  n_out <- sum(cluster_data[[var]] < q1 - 1.5*iqr |
                 cluster_data[[var]] > q3 + 1.5*iqr)
  
  outlier_table <- rbind(outlier_table, data.frame(
    Variable     = var,
    n_outliers   = n_out,
    pct_outliers = round(100 * n_out / nrow(cluster_data), 2)
  ))
}

# ===== SAVE ALL TABLES =====

write.csv(normality_table,
          paste0(viz_folder, "/Supp_Assumptions_Normality.csv"),
          row.names = FALSE)

write.csv(cor_df,
          paste0(viz_folder, "/Supp_Assumptions_Correlations.csv"),
          row.names = FALSE)

write.csv(outlier_table,
          paste0(viz_folder, "/Supp_Assumptions_Outliers.csv"),
          row.names = FALSE)

cat("Saved:\n")
cat("  Supp_Assumptions_Normality.csv\n")
cat("  Supp_Assumptions_Correlations.csv\n")
cat("  Supp_Assumptions_Outliers.csv\n\n")

cat("=== NORMALITY ===\n")
print(normality_table, row.names = FALSE)

cat("\n=== CORRELATIONS ===\n")
print(cor_df, row.names = FALSE)

cat("\n=== OUTLIERS ===\n")
print(outlier_table, row.names = FALSE)

# PCA biplots, k=2 / k=3 / k=4 side by side -> Figure 2
############looking at k3###############
library(mclust)
library(cluster)
library(ggplot2)
library(ggdist)
library(patchwork)
library(dplyr)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"

# ===== LOAD AND PREP =====

events_filtered <- read.csv(
  paste0(viz_folder, "/events_dataset_cluster3_removed.csv"),
  stringsAsFactors = FALSE
)

metrics      <- c("duration_min", "auc", "time_to_peak")
cluster_data <- events_filtered %>%
  filter(!is.na(duration_min), !is.na(auc), !is.na(time_to_peak))

cluster_data_log <- data.frame(
  log_duration     = log(cluster_data$duration_min + 0.1),
  log_auc          = log(cluster_data$auc          + 0.1),
  log_time_to_peak = log(cluster_data$time_to_peak + 0.1)
)

cluster_data_log_scaled <- scale(cluster_data_log)

# ===== FIT k=2, k=3, k=4 =====

set.seed(123)
gmm_k2 <- Mclust(cluster_data_log_scaled, G = 2)
set.seed(123)
gmm_k3 <- Mclust(cluster_data_log_scaled, G = 3)
set.seed(123)
gmm_k4 <- Mclust(cluster_data_log_scaled, G = 4)

# ===== ORDER CLUSTERS BY DURATION =====

for (k_val in c(2, 3, 4)) {
  gmm_obj <- get(paste0("gmm_k", k_val))
  cluster_data[[paste0("k", k_val)]] <- gmm_obj$classification
  
  k_order <- cluster_data %>%
    group_by(.data[[paste0("k", k_val)]]) %>%
    summarise(mean_dur = mean(duration_min), .groups = "drop") %>%
    arrange(mean_dur) %>%
    mutate(label = paste0("Cluster ", row_number()))
  
  cluster_data <- cluster_data %>%
    left_join(
      k_order %>% rename(!!paste0("k", k_val, "_label") := label),
      by = paste0("k", k_val)
    ) %>%
    select(-mean_dur)
}

# ===== PCA =====

pca_result <- prcomp(cluster_data_log_scaled, scale. = FALSE)
var_exp    <- round(100 * summary(pca_result)$importance[2, 1:2], 1)

# Build display labels with n for each solution
make_display_labels <- function(label_col, colors) {
  sizes  <- table(cluster_data[[label_col]])
  labels <- paste0(names(sizes), " (n=", as.integer(sizes), ")")
  names(labels)  <- names(sizes)
  named_colors   <- setNames(colors, labels)
  list(labels = labels, colors = named_colors)
}

k2_info <- make_display_labels("k2_label",
                               c("#1D9E75", "#D85A30"))
k3_info <- make_display_labels("k3_label",
                               c("#1D9E75", "#378ADD", "#D85A30"))
k4_info <- make_display_labels("k4_label",
                               c("#1D9E75", "#378ADD", "#EF9F27", "#D85A30"))

cluster_data$k2_display <- k2_info$labels[cluster_data$k2_label]
cluster_data$k3_display <- k3_info$labels[cluster_data$k3_label]
cluster_data$k4_display <- k4_info$labels[cluster_data$k4_label]

pc_df <- data.frame(
  PC1        = pca_result$x[, 1],
  PC2        = pca_result$x[, 2],
  k2_display = factor(cluster_data$k2_display,
                      levels = k2_info$labels),
  k3_display = factor(cluster_data$k3_display,
                      levels = k3_info$labels),
  k4_display = factor(cluster_data$k4_display,
                      levels = k4_info$labels)
)

# ===== THREE-PANEL PCA =====

make_pca_plot <- function(pc_df, display_col, colors,
                          title, subtitle) {
  ggplot(pc_df, aes(x = PC1, y = PC2,
                    color = .data[[display_col]],
                    fill  = .data[[display_col]])) +
    geom_point(alpha = 0.45, size = 1.8) +
    stat_ellipse(level = 0.95, geom = "polygon",
                 alpha = 0.12, linewidth = 0.9) +
    scale_color_manual(values = colors, name = NULL) +
    scale_fill_manual(values  = colors, name = NULL) +
    labs(
      title    = title,
      subtitle = subtitle,
      x = paste0("PC1 (", var_exp[1], "%)"),
      y = paste0("PC2 (", var_exp[2], "%)")
    ) +
    theme_minimal(base_size = 11) +
    theme(
      legend.position  = "bottom",
      legend.text      = element_text(size = 8),
      panel.grid.minor = element_blank(),
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
      plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "gray50")
    )
}

p_k2 <- make_pca_plot(pc_df, "k2_display", k2_info$colors,
                      "k = 2",
                      paste0("BIC = ", round(gmm_k2$bic, 0), " | Silhouette = 0.417"))

p_k3 <- make_pca_plot(pc_df, "k3_display", k3_info$colors,
                      "k = 3",
                      paste0("BIC = ", round(gmm_k3$bic, 0), " | Silhouette = 0.152"))

p_k4 <- make_pca_plot(pc_df, "k4_display", k4_info$colors,
                      "k = 4",
                      paste0("BIC = ", round(gmm_k4$bic, 0), " | Silhouette = 0.192"))

p_combined <- p_k2 + p_k3 + p_k4 +
  plot_layout(ncol = 3) +
  plot_annotation(
    title    = "Gaussian Mixture Model cluster solutions",
    subtitle = paste0("PC1: ", var_exp[1], "%, PC2: ", var_exp[2],
                      "% (", sum(var_exp), "% total variance) | 95% confidence ellipses"),
    theme = theme(
      plot.title    = element_text(face = "bold", hjust = 0.5, size = 14),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray50")
    )
  )

ggsave(paste0(viz_folder, "/Fig_PCA_k2_k3_k4_noC3.png"),
       p_combined, width = 18, height = 6, dpi = 300)

cat("Saved: Fig_PCA_k2_k3_k4_noC3.png\n\n")

# ===== DESCRIPTIVE STATS FOR k=3 =====

cat("========================================\n")
cat("   k=3 CLUSTER SIZES AND DESCRIPTIVES\n")
cat("========================================\n\n")

k3_labels_ordered <- paste0("Cluster ", 1:3)

cat("Cluster sizes:\n")
for (cl in k3_labels_ordered) {
  n <- sum(cluster_data$k3_label == cl)
  cat(" ", cl, ": n =", n,
      "(", round(100*n/nrow(cluster_data), 1), "%)\n")
}

test_vars <- c("duration_min", "auc", "time_to_peak")
test_labs <- c("Duration (min)", "AUC (mhos-min)", "Time to Peak (min)")

cat("\n")
for (i in seq_along(test_vars)) {
  var <- test_vars[i]
  lab <- test_labs[i]
  if (!var %in% colnames(cluster_data)) next
  cat(lab, ":\n")
  for (cl in k3_labels_ordered) {
    vals <- cluster_data[[var]][cluster_data$k3_label == cl]
    cat("  ", cl, ": M =", round(mean(vals, na.rm=TRUE), 2),
        ", SD =", round(sd(vals, na.rm=TRUE), 2), "\n")
  }
}

# recalculated cluster stats on the final sample -> Table 2
##########recalculate clustering characteristics################################
library(mclust)
library(dplyr)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"

# ===== LOAD AND PREP =====

events_filtered <- read.csv(
  paste0(viz_folder, "/events_dataset_cluster3_removed.csv"),
  stringsAsFactors = FALSE
)

metrics      <- c("duration_min", "auc", "time_to_peak")
cluster_data <- events_filtered %>%
  filter(!is.na(duration_min), !is.na(auc), !is.na(time_to_peak))

cat("N events:", nrow(cluster_data), "\n\n")

cluster_data_log <- data.frame(
  log_duration     = log(cluster_data$duration_min + 0.1),
  log_auc          = log(cluster_data$auc          + 0.1),
  log_time_to_peak = log(cluster_data$time_to_peak + 0.1)
)

cluster_data_log_scaled <- scale(cluster_data_log)

# ===== REFIT GMM k=4 =====

set.seed(123)
gmm_k4 <- Mclust(cluster_data_log_scaled, G = 4)

cluster_data$k4 <- gmm_k4$classification

k4_order <- cluster_data %>%
  group_by(k4) %>%
  summarise(mean_dur = mean(duration_min), .groups = "drop") %>%
  arrange(mean_dur) %>%
  mutate(k4_label = paste0("Cluster ", row_number()))

cluster_data <- cluster_data %>%
  left_join(k4_order[, c("k4", "k4_label")], by = "k4")

k4_labels_ordered <- paste0("Cluster ", 1:4)

# ===== DESCRIPTIVE STATS (3 clustering metrics only) =====

test_vars <- c("duration_min", "auc", "time_to_peak")
test_labs <- c("Duration (min)", "AUC (mhos-min)", "Time to Peak (min)")

cat("========================================\n")
cat("   DESCRIPTIVE STATS BY CLUSTER\n")
cat("   (duration, AUC, time to peak only)\n")
cat("========================================\n\n")

cat("Cluster sizes:\n")
for (cl in k4_labels_ordered) {
  n <- sum(cluster_data$k4_label == cl)
  cat(" ", cl, ": n =", n,
      "(", round(100 * n / nrow(cluster_data), 1), "%)\n")
}

cat("\n")

for (i in seq_along(test_vars)) {
  var <- test_vars[i]
  lab <- test_labs[i]
  cat(lab, ":\n")
  for (cl in k4_labels_ordered) {
    vals <- cluster_data[[var]][cluster_data$k4_label == cl]
    cat("  ", cl, ": M =", round(mean(vals, na.rm = TRUE), 2),
        ", SD =", round(sd(vals, na.rm = TRUE), 2), "\n")
  }
  cat("\n")
}

# ===== SAVE =====

# Descriptive stats table
desc_rows <- data.frame()

sizes_row <- c("n (%)", sapply(k4_labels_ordered, function(cl) {
  n <- sum(cluster_data$k4_label == cl)
  paste0(n, " (", round(100 * n / nrow(cluster_data), 1), "%)")
}))

for (i in seq_along(test_vars)) {
  var <- test_vars[i]
  lab <- test_labs[i]
  
  row_vals <- sapply(k4_labels_ordered, function(cl) {
    vals <- cluster_data[[var]][cluster_data$k4_label == cl]
    paste0(round(mean(vals, na.rm = TRUE), 2),
           " (", round(sd(vals, na.rm = TRUE), 2), ")")
  })
  
  desc_rows <- rbind(desc_rows, c(Variable = lab, row_vals))
}

desc_table <- rbind(sizes_row, desc_rows)
colnames(desc_table) <- c("Variable",
                          "Cluster 1 (shortest)",
                          "Cluster 2",
                          "Cluster 3",
                          "Cluster 4 (longest)")

write.csv(desc_table,
          paste0(viz_folder, "/GMM_noC3_k4_DescStats_3metrics.csv"),
          row.names = FALSE)

cat("\nSaved: GMM_noC3_k4_DescStats_3metrics.csv\n")

# grabbing a few representative events per cluster to plot as example traces
###########choosing k4 after all. selecting representative events for visualization##############
library(dplyr)
library(ggplot2)
library(patchwork)
library(lubridate)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"
eda_path   <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min_trimmed.csv"

# ===== LOAD DATA =====

events_clustered <- read.csv(
  paste0(viz_folder, "/events_dataset_cluster3_removed.csv"),
  stringsAsFactors = FALSE
)

eda_data <- read.csv(eda_path, stringsAsFactors = FALSE)
eda_data$timestamp <- ymd_hms(eda_data$timestamp)

# ===== REFIT GMM AND ASSIGN k=4 LABELS =====

library(mclust)

metrics      <- c("duration_min", "auc", "time_to_peak")
cluster_data <- events_clustered %>%
  filter(!is.na(duration_min), !is.na(auc), !is.na(time_to_peak))

cluster_data_log <- data.frame(
  log_duration     = log(cluster_data$duration_min + 0.1),
  log_auc          = log(cluster_data$auc          + 0.1),
  log_time_to_peak = log(cluster_data$time_to_peak + 0.1)
)

cluster_data_log_scaled <- scale(cluster_data_log)

set.seed(123)
gmm_k4 <- Mclust(cluster_data_log_scaled, G = 4)

cluster_data$k4 <- gmm_k4$classification

k4_order <- cluster_data %>%
  group_by(k4) %>%
  summarise(mean_dur = mean(duration_min), .groups = "drop") %>%
  arrange(mean_dur) %>%
  mutate(k4_label = paste0("Cluster ", row_number()))

cluster_data <- cluster_data %>%
  left_join(k4_order[, c("k4", "k4_label")], by = "k4")

# Cluster means for reference
cat("Cluster means (duration / time to peak):\n")
cluster_data %>%
  group_by(k4_label) %>%
  summarise(
    n        = n(),
    mean_dur = round(mean(duration_min), 2),
    mean_ttp = round(mean(time_to_peak), 2),
    mean_auc = round(mean(auc), 2),
    .groups  = "drop"
  ) %>%
  print()

# ===== SELECT REPRESENTATIVE EVENTS PER CLUSTER =====
# Find events closest to cluster mean on duration AND time to peak

n_examples <- 3  # number of example traces per cluster

representative_events <- cluster_data %>%
  filter(!is.na(event_start)) %>%
  group_by(k4_label) %>%
  mutate(
    mean_dur = mean(duration_min),
    mean_ttp = mean(time_to_peak),
    dist_to_mean = sqrt(
      ((duration_min  - mean_dur) / sd(duration_min))^2 +
        ((time_to_peak  - mean_ttp) / sd(time_to_peak))^2
    )
  ) %>%
  arrange(dist_to_mean) %>%
  slice_head(n = n_examples) %>%
  ungroup() %>%
  select(k4_label, participant, visit, night,
         event_start, event_end,
         duration_min, time_to_peak, auc)

cat("\nSelected representative events:\n")
print(as.data.frame(representative_events), row.names = FALSE)

# ===== PULL RAW EDA SIGNAL FOR EACH EVENT =====

pull_event_signal <- function(participant_id, visit_id, night_id,
                              start_str, end_str, eda_df) {
  
  start_dt <- ymd_hms(start_str)
  end_dt   <- ymd_hms(end_str)
  
  # Add small buffer around event for context
  buffer_mins <- 5
  window_start <- start_dt - minutes(buffer_mins)
  window_end   <- end_dt   + minutes(buffer_mins)
  
  signal <- eda_df %>%
    filter(
      participant == participant_id,
      visit       == visit_id,
      night       == night_id,
      timestamp   >= window_start,
      timestamp   <= window_end
    ) %>%
    mutate(
      minutes_from_onset = as.numeric(
        difftime(timestamp, start_dt, units = "mins")
      ),
      in_event = timestamp >= start_dt & timestamp <= end_dt
    )
  
  return(signal)
}

# Build signal dataset for all representative events
all_signals <- data.frame()

for (i in 1:nrow(representative_events)) {
  ev <- representative_events[i, ]
  
  signal <- pull_event_signal(
    participant_id = ev$participant,
    visit_id       = ev$visit,
    night_id       = ev$night,
    start_str      = ev$event_start,
    end_str        = ev$event_end,
    eda_df         = eda_data
  )
  
  if (nrow(signal) == 0) {
    cat("No signal found for event", i, "-", ev$participant, ev$visit, "\n")
    next
  }
  
  signal$cluster_label  <- ev$k4_label
  signal$event_id       <- i
  signal$example_num    <- which(
    representative_events$k4_label == ev$k4_label &
      representative_events$event_start == ev$event_start
  )
  signal$duration_label <- paste0("dur=", round(ev$duration_min, 0), "min")
  signal$ttp_label      <- paste0("TTP=", round(ev$time_to_peak, 1), "min")
  
  all_signals <- rbind(all_signals, signal)
}

cat("\nSignal rows pulled:", nrow(all_signals), "\n")
cat("Events with signal:", length(unique(all_signals$event_id)), "\n\n")

# ===== CLUSTER LABELS FOR PLOT TITLES =====

# NOTE: relabeled to match final manuscript terminology (Brief/Longer HF,
# Camelback HF, Night Sweat) - an earlier pass through this section had
# "Camelback" and "Prolonged HF" on the wrong clusters (2 and 3, respectively)
# relative to what ended up in the paper. Cluster rank (by mean duration)
# was always correct; only the display labels below were stale.
cluster_descriptions <- c(
  "Cluster 1" = "Cluster 1: Brief Hot Flash\n(very short TTP)",
  "Cluster 2" = "Cluster 2: Longer Hot Flash\n(brief initial rise + main peak)",
  "Cluster 3" = "Cluster 3: Camelback Hot Flash\n(double-peak morphology)",
  "Cluster 4" = "Cluster 4: Night Sweat\n(slow sustained onset)"
)

all_signals$cluster_description <- cluster_descriptions[
  all_signals$cluster_label
]

# ===== PLOT =====

colors_k4 <- c(
  "Cluster 1: Brief Hot Flash\n(very short TTP)"                = "#1D9E75",
  "Cluster 2: Longer Hot Flash\n(brief initial rise + main peak)" = "#378ADD",
  "Cluster 3: Camelback Hot Flash\n(double-peak morphology)"      = "#EF9F27",
  "Cluster 4: Night Sweat\n(slow sustained onset)"                = "#D85A30"
)

# One panel per cluster, examples overlaid
p_examples <- ggplot(all_signals,
                     aes(x = minutes_from_onset, y = mhos,
                         group = event_id,
                         color = cluster_description,
                         alpha = in_event)) +
  geom_rect(
    aes(xmin = 0, xmax = duration_label |> gsub("dur=", "", x = _) |>
          gsub("min", "", x = _) |> as.numeric(),
        ymin = -Inf, ymax = Inf),
    fill = "gray90", alpha = 0.3,
    inherit.aes = FALSE,
    data = all_signals %>%
      distinct(cluster_label, event_id, duration_label, cluster_description) %>%
      mutate(dur_num = as.numeric(gsub("min", "",
                                       gsub("dur=", "", duration_label))))
  ) +
  geom_line(linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed",
             color = "gray40", linewidth = 0.6) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4),
                     guide = "none") +
  scale_color_manual(values = colors_k4, guide = "none") +
  facet_wrap(~ cluster_description, ncol = 2, scales = "free_y") +
  labs(
    title    = "Representative event traces by GMM cluster",
    subtitle = "Dashed line = event onset | Shaded = event window | Each line = one event",
    x        = "Minutes from event onset",
    y        = "EDA (µmhos)"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", hjust = 0.5, size = 13),
    plot.subtitle    = element_text(hjust = 0.5, size = 9, color = "gray50"),
    strip.text       = element_text(face = "bold", size = 10),
    panel.grid.minor = element_blank(),
    panel.spacing    = unit(1.5, "lines")
  )

ggsave(paste0(viz_folder, "/Fig_Example_Event_Traces.png"),
       p_examples, width = 12, height = 9, dpi = 300)

cat("Saved: Fig_Example_Event_Traces.png\n")

# ===== ALSO SAVE THE REPRESENTATIVE EVENTS TABLE =====

write.csv(representative_events,
          paste0(viz_folder, "/Representative_Events_per_Cluster.csv"),
          row.names = FALSE)

cat("Saved: Representative_Events_per_Cluster.csv\n")

# individual trace plots for candidate events - this is where I picked the
# actual events that ended up in Figure 3
##################candidate event visualizations, separate######################
library(dplyr)
library(ggplot2)
library(lubridate)
library(mclust)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"
eda_path   <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/EDA_rolling_cv_10min_trimmed.csv"

# Create subfolder for individual plots
traces_folder <- paste0(viz_folder, "/Event_Traces")
dir.create(traces_folder, showWarnings = FALSE)

# ===== LOAD DATA =====

events_clustered <- read.csv(
  paste0(viz_folder, "/events_dataset_cluster3_removed.csv"),
  stringsAsFactors = FALSE
)

eda_data <- read.csv(eda_path, stringsAsFactors = FALSE)
eda_data$timestamp <- ymd_hms(eda_data$timestamp)

# ===== REFIT GMM k=4 =====

metrics      <- c("duration_min", "auc", "time_to_peak")
cluster_data <- events_clustered %>%
  filter(!is.na(duration_min), !is.na(auc), !is.na(time_to_peak))

cluster_data_log <- data.frame(
  log_duration     = log(cluster_data$duration_min + 0.1),
  log_auc          = log(cluster_data$auc          + 0.1),
  log_time_to_peak = log(cluster_data$time_to_peak + 0.1)
)

cluster_data_log_scaled <- scale(cluster_data_log)

set.seed(123)
gmm_k4 <- Mclust(cluster_data_log_scaled, G = 4)

cluster_data$k4 <- gmm_k4$classification

k4_order <- cluster_data %>%
  group_by(k4) %>%
  summarise(mean_dur = mean(duration_min), .groups = "drop") %>%
  arrange(mean_dur) %>%
  mutate(k4_label = paste0("Cluster ", row_number()))

cluster_data <- cluster_data %>%
  left_join(k4_order[, c("k4", "k4_label")], by = "k4")

# ===== SELECT 5 REPRESENTATIVE EVENTS PER CLUSTER =====
# Closest to cluster mean on duration AND time to peak (z-score distance)

n_examples <- 5

representative_events <- cluster_data %>%
  filter(!is.na(event_start), !is.na(participant), !is.na(visit)) %>%
  group_by(k4_label) %>%
  mutate(
    mean_dur  = mean(duration_min,  na.rm = TRUE),
    sd_dur    = sd(duration_min,    na.rm = TRUE),
    mean_ttp  = mean(time_to_peak,  na.rm = TRUE),
    sd_ttp    = sd(time_to_peak,    na.rm = TRUE),
    dist_to_mean = sqrt(
      ((duration_min - mean_dur) / sd_dur)^2 +
        ((time_to_peak - mean_ttp) / sd_ttp)^2
    )
  ) %>%
  arrange(dist_to_mean) %>%
  slice_head(n = n_examples) %>%
  ungroup() %>%
  mutate(example_num = ave(seq_along(k4_label),
                           k4_label, FUN = seq_along)) %>%
  select(k4_label, example_num, participant, visit, night,
         event_start, event_end,
         duration_min, time_to_peak, auc)

cat("Representative events selected:\n")
print(as.data.frame(
  representative_events %>%
    select(k4_label, example_num, participant, visit,
           duration_min, time_to_peak)
), row.names = FALSE)

# ===== CLUSTER COLORS AND DESCRIPTIONS =====

colors_k4 <- c(
  "Cluster 1" = "#1D9E75",
  "Cluster 2" = "#378ADD",
  "Cluster 3" = "#EF9F27",
  "Cluster 4" = "#D85A30"
)

# ===== PULL SIGNAL AND PLOT EACH EVENT =====

buffer_mins <- 5
saved_count <- 0
failed_count <- 0

for (i in 1:nrow(representative_events)) {
  
  ev <- representative_events[i, ]
  
  start_dt     <- ymd_hms(ev$event_start)
  end_dt       <- ymd_hms(ev$event_end)
  window_start <- start_dt - minutes(buffer_mins)
  window_end   <- end_dt   + minutes(buffer_mins)
  
  signal <- eda_data %>%
    filter(
      participant == ev$participant,
      visit       == ev$visit,
      night       == ev$night,
      timestamp   >= window_start,
      timestamp   <= window_end
    ) %>%
    mutate(
      minutes_from_onset = as.numeric(
        difftime(timestamp, start_dt, units = "mins")
      ),
      in_event = timestamp >= start_dt & timestamp <= end_dt
    )
  
  if (nrow(signal) == 0) {
    cat("WARNING: No signal found for",
        ev$participant, ev$visit, "night", ev$night,
        "| event", i, "\n")
    failed_count <- failed_count + 1
    next
  }
  
  # Y axis range — zoom to signal, add 10% padding
  y_min <- min(signal$mhos, na.rm = TRUE)
  y_max <- max(signal$mhos, na.rm = TRUE)
  y_pad <- (y_max - y_min) * 0.12
  y_lo  <- max(0, y_min - y_pad)
  y_hi  <- y_max + y_pad
  
  # X axis: from -buffer to event end + buffer (in minutes)
  x_lo <- -buffer_mins - 0.5
  x_hi <- ev$duration_min + buffer_mins + 0.5
  
  cluster_color <- colors_k4[ev$k4_label]
  
  p <- ggplot(signal, aes(x = minutes_from_onset, y = mhos)) +
    
    # Shaded event window
    annotate("rect",
             xmin = 0, xmax = ev$duration_min,
             ymin = y_lo, ymax = y_hi,
             fill = cluster_color, alpha = 0.10) +
    
    # Pre/post buffer shading
    annotate("rect",
             xmin = x_lo, xmax = 0,
             ymin = y_lo, ymax = y_hi,
             fill = "gray80", alpha = 0.20) +
    annotate("rect",
             xmin = ev$duration_min, xmax = x_hi,
             ymin = y_lo, ymax = y_hi,
             fill = "gray80", alpha = 0.20) +
    
    # EDA line — lighter outside event, full color inside
    geom_line(data = signal %>% filter(!in_event),
              color = "gray60", linewidth = 0.8) +
    geom_line(data = signal %>% filter(in_event),
              color = cluster_color, linewidth = 1.2) +
    
    # Event onset and offset lines
    geom_vline(xintercept = 0,
               linetype = "dashed", color = "gray30",
               linewidth = 0.7) +
    geom_vline(xintercept = ev$duration_min,
               linetype = "dashed", color = "gray30",
               linewidth = 0.7) +
    
    # Peak marker
    geom_vline(xintercept = ev$time_to_peak,
               linetype = "dotted", color = cluster_color,
               linewidth = 0.9) +
    
    # Annotations
    annotate("text",
             x = 0.3, y = y_hi - y_pad * 0.3,
             label = "onset", size = 3,
             color = "gray40", hjust = 0, fontface = "italic") +
    annotate("text",
             x = ev$time_to_peak + 0.3, y = y_hi - y_pad * 0.3,
             label = "peak", size = 3,
             color = cluster_color, hjust = 0, fontface = "italic") +
    annotate("text",
             x = ev$duration_min + 0.3, y = y_hi - y_pad * 0.3,
             label = "offset", size = 3,
             color = "gray40", hjust = 0, fontface = "italic") +
    
    # Stats in corner
    annotate("text",
             x = x_hi - 0.5, y = y_lo + y_pad * 0.5,
             label = paste0("Duration: ", round(ev$duration_min, 1), " min\n",
                            "Time to peak: ", round(ev$time_to_peak, 1), " min\n",
                            "AUC: ", round(ev$auc, 1), " µmhos·min"),
             size = 3, hjust = 1, vjust = 0,
             color = "gray30", fontface = "italic") +
    
    scale_x_continuous(
      breaks = seq(floor(x_lo), ceiling(x_hi), by = 5),
      limits = c(x_lo, x_hi)
    ) +
    scale_y_continuous(limits = c(y_lo, y_hi)) +
    
    labs(
      title    = paste0(ev$k4_label,
                        "  |  ", ev$participant,
                        "  |  ", ev$visit),
      subtitle = paste0("Example ", ev$example_num, " of ", n_examples,
                        " representative events"),
      x        = "Minutes from event onset",
      y        = "EDA (µmhos)"
    ) +
    
    theme_minimal(base_size = 12) +
    theme(
      plot.title       = element_text(
        face = "bold", hjust = 0,
        size = 13,
        color = cluster_color),
      plot.subtitle    = element_text(hjust = 0, size = 9,
                                      color = "gray50"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "gray90"),
      axis.title       = element_text(size = 11),
      axis.text        = element_text(size = 10)
    )
  
  # Save with informative filename
  filename <- paste0(
    traces_folder, "/",
    gsub(" ", "_", ev$k4_label), "_",
    sprintf("%02d", as.integer(ev$example_num)), "_",
    ev$participant, "_", ev$visit, ".png"
  )
  
  ggsave(filename, p, width = 9, height = 5, dpi = 300)
  saved_count <- saved_count + 1
  cat("Saved:", basename(filename), "\n")
}

cat("\n========================================\n")
cat("Saved:", saved_count, "plots\n")
cat("Failed:", failed_count, "events (no signal found)\n")
cat("Location:", traces_folder, "\n")
cat("========================================\n")

# sanity check on what got dropped when cluster 3 got removed
###############quick analysis of what i deleted when i deleted c3###########
library(dplyr)
library(mclust)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"

# Load the original labeled dataset (before removal)
events_original <- read.csv(
  paste0(viz_folder, "/events_dataset_GMM3var_k2_k4_clusters.csv"),
  stringsAsFactors = FALSE
)

# Load the filtered dataset
events_filtered <- read.csv(
  paste0(viz_folder, "/events_dataset_cluster3_removed.csv"),
  stringsAsFactors = FALSE
)

cat("Original N:", nrow(events_original), "\n")
cat("Filtered N:", nrow(events_filtered), "\n")
cat("Removed N:", nrow(events_original) - nrow(events_filtered), "\n\n")

# Get the removed events
removed_events <- events_original %>%
  filter(!is.na(gmm_k4_label)) %>%
  anti_join(events_filtered %>% select(participant, visit, night, event_start),
            by = c("participant", "visit", "night", "event_start"))

cat("Confirmed removed N:", nrow(removed_events), "\n\n")

# Characterize removed events
test_vars <- c("duration_min", "auc", "time_to_peak")
test_labs <- c("Duration (min)", "AUC (mhos-min)", "Time to Peak (min)")

cat("========================================\n")
cat("   REMOVED EVENTS CHARACTERISTICS\n")
cat("========================================\n\n")

removed_stats <- data.frame()

for (i in seq_along(test_vars)) {
  var <- test_vars[i]
  lab <- test_labs[i]
  if (!var %in% colnames(removed_events)) next
  
  vals <- removed_events[[var]]
  removed_stats <- rbind(removed_stats, data.frame(
    Variable = lab,
    N        = sum(!is.na(vals)),
    Mean     = round(mean(vals, na.rm = TRUE), 2),
    SD       = round(sd(vals,   na.rm = TRUE), 2),
    Median   = round(median(vals, na.rm = TRUE), 2),
    Min      = round(min(vals, na.rm = TRUE), 2),
    Max      = round(max(vals, na.rm = TRUE), 2)
  ))
  
  cat(lab, ": M =", round(mean(vals, na.rm=TRUE), 2),
      ", SD =", round(sd(vals, na.rm=TRUE), 2),
      ", Median =", round(median(vals, na.rm=TRUE), 2),
      ", Range =", round(min(vals, na.rm=TRUE), 2),
      "-", round(max(vals, na.rm=TRUE), 2), "\n")
}

write.csv(removed_stats,
          paste0(viz_folder, "/Removed_Cluster3_Characteristics.csv"),
          row.names = FALSE)

# ===== COMPARE REMOVED TO RETAINED =====
#  see the removed group looks different

cat("\n\n========================================\n")
cat("   REMOVED vs RETAINED MEANS\n")
cat("========================================\n\n")

events_original$removed <- !(events_original$event_start %in%
                               events_filtered$event_start)

comparison <- data.frame()

for (i in seq_along(test_vars)) {
  var <- test_vars[i]
  lab <- test_labs[i]
  if (!var %in% colnames(events_original)) next
  
  removed_vals  <- events_original[[var]][events_original$removed == TRUE  & !is.na(events_original[[var]])]
  retained_vals <- events_original[[var]][events_original$removed == FALSE & !is.na(events_original[[var]])]
  
  comparison <- rbind(comparison, data.frame(
    Variable       = lab,
    Removed_Mean   = round(mean(removed_vals), 2),
    Removed_SD     = round(sd(removed_vals),   2),
    Retained_Mean  = round(mean(retained_vals), 2),
    Retained_SD    = round(sd(retained_vals),   2)
  ))
  
  cat(lab, ":\n")
  cat("  Removed:  M =", round(mean(removed_vals), 2),
      ", SD =", round(sd(removed_vals), 2), "\n")
  cat("  Retained: M =", round(mean(retained_vals), 2),
      ", SD =", round(sd(retained_vals), 2), "\n\n")
}

write.csv(comparison,
          paste0(viz_folder, "/Removed_vs_Retained_Comparison.csv"),
          row.names = FALSE)

# ===== PARTICIPANT DISTRIBUTION =====

cat("========================================\n")
cat("   HOW MANY PARTICIPANTS CONTRIBUTED\n")
cat("   REMOVED EVENTS?\n")
cat("========================================\n\n")

participant_summary <- removed_events %>%
  group_by(participant) %>%
  summarise(n_removed = n(), .groups = "drop") %>%
  arrange(desc(n_removed))

cat("Total participants with removed events:",
    nrow(participant_summary), "\n")
cat("Mean removed events per participant:",
    round(mean(participant_summary$n_removed), 1), "\n")
cat("Range:", min(participant_summary$n_removed),
    "-", max(participant_summary$n_removed), "\n\n")

print(as.data.frame(participant_summary), row.names = FALSE)

write.csv(participant_summary,
          paste0(viz_folder, "/Removed_Events_By_Participant.csv"),
          row.names = FALSE)

# full participant/visit/night exclusion table -> Supplemental Table 1
#####exclusions CSV##########
library(dplyr)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"

# ===== BUILD EXCLUSION TABLE =====

# ---- FULL VISIT EXCLUSIONS ----

full_exclusions <- data.frame(
  Participant = c(
    "S57", "S75", "S79", "S79", "S79",
    "S91", "S91", "S91",
    "S101", "S103", "S103",
    "S131", "S131", "S131",
    "S154", "S176", "S176", "S176", "S176",
    "S198", "S242",
    "S266", "S266", "S266",
    "S86", "S282"
  ),
  Visit = c(
    "V2", "V2", "V2", "V3", "V4",
    "V2", "V3", "V4",
    "V4", "V3", "V4",
    "V2", "V3", "V4",
    "V4", "V2", "V3", "V4", "V5",
    "V3", "V2",
    "V2", "V3", "V4",
    "V4", "V4"
  ),
  Night = rep("All", 26),
  Type  = rep("Full visit excluded", 26),
  Reason = rep("Signal quality: amplitude changes exceeding 1000 µmhos or insufficient recording duration", 26)
)

# ---- SINGLE NIGHT EXCLUSIONS ----

night_exclusions <- data.frame(
  Participant = c(
    "S69", "S67", "S86", "S86", "S88",
    "S101", "S101", "S104", "S104",
    "S146", "S154", "S226", "S226",
    "S250", "S269", "S244", "S146", "S143"
  ),
  Visit = c(
    "V3", "V4", "V3", "V4", "V4",
    "V2", "V3", "V3", "V4",
    "V2", "V3", "V3", "V4",
    "V4", "V3", "V3", "V4", "V2"
  ),
  Night = c(
    "2", "2", "1", "2", "2",
    "2", "2", "1", "2",
    "1", "2", "2", "1",
    "2", "2", "1", "1", "1"
  ),
  Type   = rep("Single night excluded", 18),
  Reason = rep("Signal quality: amplitude changes exceeding 1000 µmhos, insufficient recording duration, or signal dropout", 18)
)

# ---- MISSING START TIME ----

missing_start <- data.frame(
  Participant = c("S51", "S37", "S37"),
  Visit       = c("V2", "V3", "V2"),
  Night       = rep("All", 3),
  Type        = rep("Missing EDA start time", 3),
  Reason      = rep("No EDA monitor initialization time recorded in REDCap; timestamps could not be assigned", 3)
)

# ---- TIME TRIMS ----
# Each row = one trim applied to one night

time_trims <- data.frame(
  Participant = c(
    "S64",  "S94",  "S94",  "S88",  "S88",
    "S86",  "S75",  "S70",  "S70",  "S70",
    "S69",  "S295", "S293", "S292", "S288",
    "S288", "S274", "S270", "S260", "S260",
    "S250", "S250", "S238", "S226", "S226",
    "S225", "S225", "S213", "S213", "S213",
    "S213", "S198", "S196", "S196", "S193",
    "S193", "S193", "S166", "S154", "S148",
    "S148", "S148", "S146", "S143", "S143",
    "S143", "S122", "S122", "S122", "S113",
    "S104", "S104"
  ),
  Visit = c(
    "V3",  "V3",  "V3",  "V3",  "V2",
    "V3",  "V4",  "V3",  "V3",  "V2",
    "V2",  "V5",  "V3",  "V2",  "V4",
    "V3",  "V3",  "V2",  "V2",  "V2",
    "V2",  "V2",  "V3",  "V2",  "V2",
    "V3",  "V3",  "V4",  "V4",  "V3",
    "V2",  "V2",  "V3",  "V2",  "V3",
    "V2",  "V2",  "V3",  "V2",  "V2",
    "V2",  "V2",  "V4",  "V4",  "V3",
    "V2",  "V4",  "V3",  "V2",  "V2",
    "V4",  "V2"
  ),
  Night = c(
    "1", "1", "2", "1", "1",
    "2", "2", "1", "2", "1",
    "2", "1", "1", "1", "1",
    "2", "1", "2", "1", "2",
    "1", "1", "1", "1", "2",
    "1", "2", "1", "1", "1",
    "2", "1", "2", "2", "2",
    "1", "2", "1", "1", "1",
    "1", "2", "2", "2", "1",
    "2", "2", "2", "2", "2",
    "1", "2"
  ),
  Type = rep("Time-based trim", 52),
  Reason = c(
    "Trimmed after 06:00 (wake artifact)",
    "Trimmed after 07:30 (wake artifact)",
    "Trimmed after 05:00 (wake artifact)",
    "Trimmed after 06:00 (wake artifact)",
    "Trimmed after 07:00 (wake artifact)",
    "Trimmed after 06:00 (wake artifact)",
    "Trimmed after 06:30 (wake artifact)",
    "Trimmed after 07:30 (wake artifact)",
    "Trimmed after 08:00 (wake artifact)",
    "Trimmed after 06:15 (wake artifact)",
    "Trimmed after 05:50 (wake artifact)",
    "Trimmed after 07:30 (wake artifact)",
    "Trimmed after 06:00 (wake artifact)",
    "Trimmed after 08:25 (wake artifact)",
    "Trimmed after 05:00 (wake artifact)",
    "Trimmed before 22:30 (pre-sleep artifact)",
    "Trimmed after 06:00 (wake artifact)",
    "Trimmed before 12:30 (mid-day artifact)",
    "Trimmed after 03:30 (early wake artifact)",
    "Trimmed after 04:00 (early wake artifact)",
    "Trimmed before 22:30 (pre-sleep artifact)",
    "Trimmed after 06:00 (wake artifact)",
    "Trimmed after 06:10 (wake artifact)",
    "Trimmed after 10:00 (wake artifact)",
    "Trimmed after 07:50 (wake artifact)",
    "Trimmed before 00:00 (pre-sleep artifact)",
    "Trimmed after 06:00 (wake artifact)",
    "Trimmed before 23:15 (pre-sleep artifact)",
    "Trimmed after 07:00 (wake artifact)",
    "Trimmed before 23:00 (pre-sleep artifact)",
    "Trimmed after 07:15 (wake artifact)",
    "Trimmed after 05:30 (wake artifact)",
    "Trimmed after 07:30 (wake artifact)",
    "Trimmed after 07:00 (wake artifact)",
    "Trimmed after 04:00 (early wake artifact)",
    "Trimmed before 23:00 (pre-sleep artifact)",
    "Trimmed after 05:00 (wake artifact)",
    "Trimmed after 07:00 (wake artifact)",
    "Trimmed after 10:50 (wake artifact)",
    "Trimmed before 22:15 (pre-sleep artifact)",
    "Trimmed after 06:00 (wake artifact)",
    "Trimmed after 06:00 (wake artifact)",
    "Trimmed after 05:15 (wake artifact)",
    "Trimmed after 05:00 (wake artifact)",
    "Trimmed after 06:30 (wake artifact)",
    "Trimmed after 02:10 (early wake artifact)",
    "Trimmed after 05:00 (wake artifact)",
    "Trimmed after 05:00 (wake artifact)",
    "Trimmed after 04:30 (early wake artifact)",
    "Trimmed before 22:30 (pre-sleep artifact)",
    "Trimmed after 05:00 (wake artifact)",
    "Trimmed after 04:45 (early wake artifact)"
  )
)

# ---- COMBINE ALL ----

exclusion_table <- bind_rows(
  missing_start,
  full_exclusions,
  night_exclusions,
  time_trims
) %>%
  arrange(Type, Participant, Visit, Night)

# ---- SUMMARY COUNTS ----

cat("========================================\n")
cat("   EXCLUSION SUMMARY\n")
cat("========================================\n\n")

cat("Missing EDA start time:", nrow(missing_start), "visits\n")
cat("Full visit exclusions: ", nrow(full_exclusions), "visits\n")
cat("Single night exclusions:", nrow(night_exclusions), "nights\n")
cat("Time-based trims:      ", nrow(time_trims), "nights\n")
cat("Total rows in table:   ", nrow(exclusion_table), "\n\n")

# ---- SAVE ----

write.csv(exclusion_table,
          paste0(viz_folder, "/Supplemental_Exclusions_Table.csv"),
          row.names = FALSE)

cat("Saved: Supplemental_Exclusions_Table.csv\n\n")

# ---- PRINT SUMMARY TABLE ----

cat("=== SUMMARY BY TYPE ===\n")
exclusion_table %>%
  group_by(Type) %>%
  summarise(N = n(), .groups = "drop") %>%
  print(row.names = FALSE)

# final remaining participant/visit/night counts
###remaining peeps and visits#########
library(dplyr)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"

events_filtered <- read.csv(
  paste0(viz_folder, "/events_dataset_cluster3_removed.csv"),
  stringsAsFactors = FALSE
)

cat("Total events:", nrow(events_filtered), "\n")
cat("Participants:", n_distinct(events_filtered$participant), "\n")
cat("Visits:", n_distinct(paste(events_filtered$participant, events_filtered$visit)), "\n")
cat("Nights:", n_distinct(paste(events_filtered$participant, events_filtered$visit, events_filtered$night)), "\n")

# final sample counts
#################counts#########################
library(dplyr)
library(mclust)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"

events_filtered <- read.csv(
  paste0(viz_folder, "/events_dataset_cluster3_removed.csv"),
  stringsAsFactors = FALSE
)

metrics      <- c("duration_min", "auc", "time_to_peak")
cluster_data <- events_filtered %>%
  filter(!is.na(duration_min), !is.na(auc), !is.na(time_to_peak))

cluster_data_log <- data.frame(
  log_duration     = log(cluster_data$duration_min + 0.1),
  log_auc          = log(cluster_data$auc          + 0.1),
  log_time_to_peak = log(cluster_data$time_to_peak + 0.1)
)

cluster_data_log_scaled <- scale(cluster_data_log)

set.seed(123)
gmm_k4 <- Mclust(cluster_data_log_scaled, G = 4)

cluster_data$k4 <- gmm_k4$classification

k4_order <- cluster_data %>%
  group_by(k4) %>%
  summarise(mean_dur = mean(duration_min), .groups = "drop") %>%
  arrange(mean_dur) %>%
  mutate(k4_label = paste0("Cluster ", row_number()))

cluster_data <- cluster_data %>%
  left_join(k4_order[, c("k4", "k4_label")], by = "k4")

cluster_data %>%
  group_by(k4_label) %>%
  summarise(
    n   = n(),
    pct = round(100 * n() / nrow(cluster_data), 1)
  ) %>%
  mutate(n_pct = paste0(n, " (", pct, "%)")) %>%
  select(k4_label, n_pct) %>%
  print(row.names = FALSE)















library(dplyr)
library(mclust)
library(cluster)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"

events_filtered <- read.csv(
  paste0(viz_folder, "/events_dataset_cluster3_removed.csv"),
  stringsAsFactors = FALSE
)

metrics      <- c("duration_min", "auc", "time_to_peak")
cluster_data <- events_filtered %>%
  filter(!is.na(duration_min), !is.na(auc), !is.na(time_to_peak))

cluster_data_log <- data.frame(
  log_duration     = log(cluster_data$duration_min + 0.1),
  log_auc          = log(cluster_data$auc          + 0.1),
  log_time_to_peak = log(cluster_data$time_to_peak + 0.1)
)

cluster_data_log_scaled <- scale(cluster_data_log)

# Fit k=1 to k=4
gmm_models <- list()

for (k in 1:4) {
  set.seed(123)
  gmm_models[[k]] <- Mclust(cluster_data_log_scaled, G = k)
}

# Silhouette for k=2 to k=4
sil_vals <- c(NA, sapply(2:4, function(k) {
  sil <- silhouette(gmm_models[[k]]$classification,
                    dist(cluster_data_log_scaled))
  round(mean(sil[, 3]), 3)
}))

# Build table
model_table <- data.frame(
  k                = 1:4,
  Model_Type       = sapply(1:4, function(k) gmm_models[[k]]$modelName),
  Parameters       = sapply(1:4, function(k) gmm_models[[k]]$df),
  Log_Likelihood   = sapply(1:4, function(k) round(gmm_models[[k]]$loglik, 2)),
  BIC              = sapply(1:4, function(k) round(gmm_models[[k]]$bic, 2)),
  AIC              = sapply(1:4, function(k) round(2 * gmm_models[[k]]$df - 2 * gmm_models[[k]]$loglik, 2)),
  Avg_Silhouette   = ifelse(is.na(sil_vals), "—", as.character(sil_vals)),
  ΔBIC             = c("—",
                       round(gmm_models[[2]]$bic - gmm_models[[1]]$bic, 2),
                       round(gmm_models[[3]]$bic - gmm_models[[2]]$bic, 2),
                       round(gmm_models[[4]]$bic - gmm_models[[3]]$bic, 2))
)

print(model_table, row.names = FALSE)

write.csv(model_table,
          paste0(viz_folder, "/Table1_Model_Fit_noC3.csv"),
          row.names = FALSE)

cat("\nSaved: Table1_Model_Fit_noC3.csv\n")


################################################################################
##deriving the formulas for HF and NS so others can use ########################
# 
#Pulls a portable classification rule out of the final k=2 GMM solution
# (the one reported in the manuscript: BIC/silhouette on the 1027-event,
# cluster-3-removed sample) so someone without R or mclust can classify a
# single detected EDA event as HF-like (short, rapid-onset) or NS-like
# (long, slow-rising) from just its three metrics: duration, AUC, time to
# peak.
#
# Then checks that the formula actually reproduces what the real GMM says,
# using the same 1027 events it was derived from.
# =============================================================================

library(dplyr)
library(mclust)

viz_folder <- "C:/Users/sofiy/OneDrive - University of Massachusetts/Desktop/Dissertation/Chapter 2/Data/Event Clustering Visualizations"

# ===== LOAD FINAL (POST-REMOVAL) DATASET =====

events_filtered <- read.csv(
  paste0(viz_folder, "/events_dataset_cluster3_removed.csv"),
  stringsAsFactors = FALSE
)

metrics <- c("duration_min", "auc", "time_to_peak")
cluster_data <- events_filtered %>%
  filter(!is.na(duration_min), !is.na(auc), !is.na(time_to_peak))

cluster_data_log <- data.frame(
  log_duration     = log(cluster_data$duration_min + 0.1),
  log_auc          = log(cluster_data$auc          + 0.1),
  log_time_to_peak = log(cluster_data$time_to_peak + 0.1)
)

cluster_data_log_scaled <- scale(cluster_data_log)

# refit k=2 with the same seed used throughout - reproduces the exact
# published solution
set.seed(123)
gmm_k2 <- Mclust(cluster_data_log_scaled, G = 2)

cat("Sanity check - should match Table 1 (k=2 row):\n")
cat("  BIC:", round(gmm_k2$bic, 2), "\n")
cat("  Cluster sizes:", paste(table(gmm_k2$classification), collapse = " / "), "\n\n")

# ===== PULL OUT EVERYTHING NEEDED TO CLASSIFY A NEW EVENT =====

log_center <- attr(cluster_data_log_scaled, "scaled:center")
log_scale  <- attr(cluster_data_log_scaled, "scaled:scale")

pi_raw    <- gmm_k2$parameters$pro
mu_raw    <- gmm_k2$parameters$mean
sigma_raw <- gmm_k2$parameters$variance$sigma

# order components by mean duration, so component 1 = HF-like (shorter),
# component 2 = NS-like (longer) - matches how clusters are labeled
# everywhere else in the pipeline
comp_order <- order(mu_raw["log_duration", ])

pi_k    <- pi_raw[comp_order]
mu_k    <- mu_raw[, comp_order]
sigma_k <- sigma_raw[, , comp_order]

# ===== PRINT THE FORMULA (PLAIN LANGUAGE + EXACT NUMBERS) =====

cat("========================================\n")
cat("   HF vs NS CLASSIFICATION FORMULA\n")
cat("========================================\n\n")

cat("Step 1 - log-transform the raw event metrics:\n")
cat("  log_duration     = log(duration_min + 0.1)\n")
cat("  log_auc          = log(auc + 0.1)\n")
cat("  log_time_to_peak = log(time_to_peak + 0.1)\n\n")

cat("Step 2 - standardize using these fixed reference values\n")
cat("(from the original 1027-event training sample):\n\n")
for (i in seq_along(log_center)) {
  cat("  z_", names(log_center)[i], " = (log_", names(log_center)[i],
      " - ", round(log_center[i], 6), ") / ", round(log_scale[i], 6), "\n", sep = "")
}

cat("\nStep 3 - mixing weights (pi):\n")
cat("  pi_HF =", round(pi_k[1], 6), "\n")
cat("  pi_NS =", round(pi_k[2], 6), "\n\n")

cat("Step 4 - cluster means (mu), in standardized log-space:\n")
print(round(mu_k, 6))

cat("\nStep 5 - cluster covariance matrices (Sigma), in standardized log-space:\n")
cat("\n-- HF component --\n")
print(round(sigma_k[, , 1], 6))
cat("\n-- NS component --\n")
print(round(sigma_k[, , 2], 6))

cat("\nStep 6 - for a new event, compute each cluster's density:\n")
cat("  density_HF = pi_HF * MVN_density(z, mean = mu_HF, cov = Sigma_HF)\n")
cat("  density_NS = pi_NS * MVN_density(z, mean = mu_NS, cov = Sigma_NS)\n\n")
cat("  predicted_type = whichever density is larger\n")
cat("  posterior_prob = density_k / (density_HF + density_NS)\n\n")

# ===== PORTABLE CLASSIFICATION FUNCTION =====
# only needs base R (solve/det) - no mclust or mvtnorm required to USE it,
# just to derive it above.

dmvnorm_manual <- function(x, mean, sigma) {
  k        <- length(x)
  diff     <- x - mean
  exponent <- -0.5 * t(diff) %*% solve(sigma) %*% diff
  denom    <- sqrt((2 * pi)^k * det(sigma))
  as.numeric(exp(exponent) / denom)
}

classify_hf_ns <- function(duration_min, auc, time_to_peak) {
  
  log_x <- c(
    log_duration     = log(duration_min + 0.1),
    log_auc          = log(auc          + 0.1),
    log_time_to_peak = log(time_to_peak + 0.1)
  )
  
  z <- (log_x - log_center) / log_scale
  
  dens <- c(
    HF = pi_k[1] * dmvnorm_manual(z, mu_k[, 1], sigma_k[, , 1]),
    NS = pi_k[2] * dmvnorm_manual(z, mu_k[, 2], sigma_k[, , 2])
  )
  
  post_prob <- dens / sum(dens)
  
  list(
    predicted_type = names(which.max(dens)),
    prob_HF        = round(post_prob["HF"], 4),
    prob_NS        = round(post_prob["NS"], 4)
  )
}

# ===== QUICK EXAMPLES =====

cat("========================================\n")
cat("   EXAMPLE CLASSIFICATIONS\n")
cat("========================================\n\n")

cat("Event: duration = 10 min, AUC = 5, time to peak = 3 min\n")
print(classify_hf_ns(10, 5, 3))

cat("\nEvent: duration = 45 min, AUC = 20, time to peak = 30 min\n")
print(classify_hf_ns(45, 20, 30))

# =============================================================================
# VALIDATION - does the formula reproduce the actual GMM's own labels?
# -----------------------------------------------------------------------------
# This isn't a test of whether the clustering itself is "correct" - it's a
# check that the formula above was pulled out correctly. It should
# reproduce the GMM's own decision boundary almost exactly, since it's the
# same math, just written out by hand instead of called through mclust.
# =============================================================================

cat("\n\n========================================\n")
cat("   VALIDATION: FORMULA vs ACTUAL GMM LABELS\n")
cat("========================================\n\n")

predicted <- mapply(
  function(d, a, t) classify_hf_ns(d, a, t)$predicted_type,
  cluster_data$duration_min,
  cluster_data$auc,
  cluster_data$time_to_peak
)

actual_type <- ifelse(gmm_k2$classification == comp_order[1], "HF", "NS")

agreement_table <- table(Formula = predicted, Actual_GMM = actual_type)

cat("Confusion matrix (formula vs actual GMM assignment):\n")
print(agreement_table)

accuracy <- sum(diag(agreement_table)) / sum(agreement_table)
cat("\nOverall agreement:", round(100 * accuracy, 2), "%\n")

if (accuracy < 0.99) {
  cat("\nAgreement isn't ~100% - double check the covariance matrices and\n")
  cat("component ordering above before trusting/sharing this formula.\n")
} else {
  cat("\nFormula reproduces the GMM's own classification almost exactly -\n")
  cat("safe to share or report as the classification rule.\n")
}

write.csv(
  data.frame(
    duration_min = cluster_data$duration_min,
    auc          = cluster_data$auc,
    time_to_peak = cluster_data$time_to_peak,
    formula_pred = predicted,
    actual_gmm   = actual_type
  ),
  paste0(viz_folder, "/HF_NS_Formula_Validation.csv"),
  row.names = FALSE
)

cat("\nSaved: HF_NS_Formula_Validation.csv\n")


###thanks for following along!### :3 