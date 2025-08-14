#--------------------------------------------
#
# mordor-spillover/R/01a-exposure_dm.R
#
# Geographic spillover of antimicrobial resistance 
# from mass distribution of azithromycin in MORDOR Niger
#
# Reading in and formatting treatment data from the MORDOR main (mortality) and morbidity trials.
#
#------------------------------------

rm(list = ls())

library(here)
source(here("R", "00-config.R"))

# Load data ---------------------------------------------------------------

# Load morbidity trial data 
#
# CSV files loaded are: "/final-gps/mordor-morbidity-cluster-gps.csv", "/final-7-12-amr/mordor-7-12-amr-analysis.csv", 
# and "/final-morbidity/Stata-mordor-Niger-long-dedupe-forYing.csv"
# morb_gps_raw contains coordinates for each grappe in the morbidity trial
# morb_amr_raw contains the grappe name and treatment assignment for morbidity trial villages
# morb_indiv_raw contains individual-level data for the morbidity trial, to compute number of children per grappe

morb_gps_raw <- read_csv(paste0(data_folder_path, "/final-gps/mordor-morbidity-cluster-gps.csv"))
morb_amr_raw <- read_csv(paste0(data_folder_path, "/final-7-12-amr/mordor-7-12-amr-analysis.csv"))
morb_indiv_raw <- read_csv(paste0(data_folder_path, "/final-morbidity/Stata-mordor-Niger-long-dedupe-forYing.csv"))

# Load main trial data 
#
# CSV files loaded are: "/final-mortality/MORDOR-treatment-history.csv" and "/final-mortality/tabulation_ne_24.csv"
# MORDOR-treatment-history.csv contains treatment history and coordinates for each grappe in the main trial
# tabulation_ne_24.csv contains individual-level data for the main trial

main_grappe_raw <- read_csv(paste0(data_folder_path, "/final-mortality/MORDOR-treatment-history.csv"))
main_indiv_raw <- read_csv(paste0(data_folder_path, "/final-mortality/tabulation_ne_24.csv"))

# Format main trial data ---------------------------------------------------

## Extracting grappe-level data (main_grappe_raw) for the main trial to get treatment assignment and coordinates

main_grappe <- main_grappe_raw %>% 
  mutate(arm = ifelse(substr(TREATMENT.HISTORY, 1, 4) == "AAAA", "azithro", "placebo")) %>%
  dplyr::select(grappe, arm, latitude, longitude) %>%
  distinct()

table(main_grappe$arm)  # Should be 303 azithro and 291 placebo grappes

## Extracting coordinates from main_grappe

main_gps <- main_grappe %>%
  dplyr::select(grappe, latitude, longitude) %>%
  distinct()

# Remove coordinates from main_grappe

main_grappe <- main_grappe %>%
  dplyr::select(-latitude, -longitude)

## Using individual-level data (main_indiv_raw) from the main trial to compute number of doses distributed per phase, per grappe

main_indiv_vars <- c("grappe", "masterperson", "phase", "dose.given.p")

table(main_indiv_raw$phase)
table(main_indiv_raw$dose.given.p)

main_indiv <- main_indiv_raw %>%
  rowwise() %>% mutate(grappe = strsplit(randu, ":")[[1]][1]) %>% ungroup() %>%
  dplyr::select(all_of(main_indiv_vars))

colSums(is.na(main_indiv)) ## No missing data

## Number of doses given per phase, by grappe

main_indiv_doses <- main_indiv %>% dplyr::filter(dose.given.p == 1) %>% group_by(grappe, phase) %>% 
  summarise(n_doses = sum(dose.given.p), .groups = "keep") %>% 
  ungroup() %>% distinct() %>% group_by(grappe) %>% summarise(n_doses = sum(n_doses), .groups = "keep") %>% ungroup()

## Number of children, by grappe
  
main_indiv_children <- main_indiv %>% dplyr::select(grappe, masterperson) %>%
  distinct() %>% group_by(grappe) %>% 
  summarise(n_children = n_distinct(masterperson)) %>% ungroup() 

## Number of children given at least 1 dose in 24 months, by grappe

main_indiv_treated <- main_indiv %>% dplyr::filter(dose.given.p == 1) %>% group_by(grappe) %>% 
  summarise(n_treated = n_distinct(masterperson), .groups = "keep") %>% 
  ungroup() %>% distinct()

## Merging grappe-level data and doses from individual-level data

main_grappe <- main_grappe %>% 
  left_join(main_indiv_treated, by = c("grappe")) %>% 
  left_join(main_indiv_children, by = c("grappe")) %>%
  left_join(main_indiv_doses, by = c("grappe")) %>%
  mutate(treat_bin = ifelse(arm == "placebo", 0, 1))


# Format morbidity trial data ---------------------------------------------

# Extracting individual-level data (morb_indiv_raw) from the morbidity trial to get number of children per grappe

# Restricting to phase < 5 rounds of MDA: Not including doses given in the 5th phase, 
# which would've been at 24 months

table(morb_indiv_raw$phase)
names(morb_indiv_raw)

morb_indiv <- morb_indiv_raw %>% dplyr::filter(phase < 24) %>%
  dplyr::select(masterperson, randomunit999) %>% distinct() %>%
  group_by(randomunit999) %>% summarise(n_children = n(), .groups = "keep") %>% ungroup() %>%
  mutate(grappe = sub('.*-', '', randomunit999)) %>% dplyr::select(-randomunit999)

# Using morb_amr_raw to get treatment assignment for the morbidity trial

morb_amr <- morb_amr_raw %>% dplyr::select("cluster_id", "arm") %>% distinct() %>%
  left_join(morb_indiv, by = c("cluster_id" = "grappe"))

# Merging morbidity trial data

morb_grappe <- morb_amr %>%
  mutate(treat_bin = ifelse(arm == "placebo", 0, 1)) %>% rename(grappe = cluster_id)

morb_gps <- morb_gps_raw %>%
  rename(grappe = gwu)

# Save data ---------------------------------------------------------------

saveRDS(main_grappe, here("data", "clean", "main_grappe.rds"))
saveRDS(morb_grappe, here("data", "clean", "morb_grappe.rds"))

# Save GPS data

saveRDS(main_gps, here("data", "clean", "main_gps.rds"))
saveRDS(morb_gps, here("data", "clean", "morb_gps.rds"))
