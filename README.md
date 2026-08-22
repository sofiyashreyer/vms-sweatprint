# vms-sweatprint

Analysis code supporting the manuscript investigating nocturnal electrodermal
activity (EDA) as an objective marker for distinguishing hot flashes (HF) from
night sweats (NS) during sleep, submitted to *Menopause*.

## What's included

- **`Analysis_Pipeline_Reorganized.R`** — the full analysis pipeline, in order:
  1. Data import (raw EDA → 1-minute summaries)
  2. Preprocessing (sleep-window trimming, exclusions, rolling CV, unit conversion)
  3. Event detection, including the threshold combinations tested before
     arriving at the final algorithm (all labeled with their parameters)
  4. Validation against participant-marked events
  5. Event metric extraction (duration, AUC, time to peak, recovery time)
  6. Variable selection (correlation check)
  7. The clustering analysis

- **`HF_NS_Classification_Formula.R`** — derives and validates a quadratic
  discriminant function (from the final two-component GMM) that classifies a
  single detected EDA event as HF-like or NS-like from its raw duration, AUC,
  and time-to-peak values. Includes a check confirming the formula reproduces
  the GMM's own classification.

## What's NOT included

No raw or de-identified participant data is included in this repository,
consistent with the terms of participant consent and IRB approval. File paths
in the scripts point to the original local directory structure and are shown
for transparency only — they will not resolve outside the original analysis
environment.

## Requirements

R (≥ 4.0) with the following packages:

```
tidyverse
dplyr
lubridate
zoo
mclust
cluster
moments
factoextra
ggplot2
ggdist
patchwork
corrplot
car
psych
readxl
```

Install any missing packages with:

```r
install.packages(c("tidyverse", "dplyr", "lubridate", "zoo", "mclust",
                    "cluster", "moments", "factoextra", "ggplot2", "ggdist",
                    "patchwork", "corrplot", "car", "psych", "readxl"))
```


## Contact

Questions about this repository or the analysis can be directed to
Sofiya Shreyer sofiyashreyer@umass.edu
