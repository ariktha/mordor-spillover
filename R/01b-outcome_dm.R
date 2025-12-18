#--------------------------------------------
#
# mordor-spillover/R/01b-outcome_dm.R
#
# Geographic spillover of antimicrobial resistance 
# from mass distribution of azithromycin in MORDOR Niger
#
# Reading in and formatting data from the MORDOR morbidity trial: 
# 24 month rectal swab DNASeq normalized abundance results for the morbidity trial
#
#------------------------------------

rm(list = ls())             

library(here)
source(here("R", "00-config.R"))
# source(here("R", "00-functions.R"))

# Load data ---------------------------------------------------------------

amr_raw <- read_csv(here("data", "untouched", "outcome", 
                         "MORDOR24_DNA_resistome.rM.class_annotated_WHG added.csv"))


# Clean data --------------------------------------------------------------

## Select and rename columns
## Check for and remove duplicates

amr_desc_cols <- c("Sample.ID", "Sample..", "Duplicate", "Phase.Id.MEP", 
              "WH.Geographic.Work.Unit..WH.Geographic.Work.Unit.Name")

amr_dat <- amr_raw %>% 
  dplyr::select(all_of(c(amr_desc_cols, ab_classes_of_interest))) %>%
  rename("sample_id" = "Sample.ID", 
         "sample" = "Sample..", 
         "dup" = "Duplicate", 
         "phase" = "Phase.Id.MEP", 
         "grappe" = "WH.Geographic.Work.Unit..WH.Geographic.Work.Unit.Name")

table(amr_dat$dup)

amr_dat <- amr_dat %>% 
  dplyr::filter(dup == 0)

## Convert to long format 

amr_dat_long <- amr_dat %>% 
  pivot_longer(cols = all_of(ab_classes_of_interest), 
               names_to = "ab_class", 
               values_to = "resistance")

amr_indiv <- amr_dat_long %>%
  dplyr::filter(ab_class == "MLS") %>%
  dplyr::select(sample_id, grappe, phase, resistance) %>%
  rename(mls_resistance = resistance)

amr_indiv_pub <- amr_indiv %>%
  dplyr::select(-sample_id)

# Add summary columns -----------------------------------------------------

##  Average resistance, SE of resistance, samples with non-zero resistance, 
##  total samples, percent samples with non-zero resistance

# se(): calculates the standard error of a numeric vector

se <- function(x) sqrt(var(x)/length(x))


amr_data_summ <- amr_dat_long %>% 
  group_by(phase, grappe, ab_class) %>%
  summarise(avg_res = mean(resistance, na.rm = FALSE),
            se_res = se(resistance),
            non_zero_res = sum(resistance != 0), 
            total_samp = n(), .groups = "keep") %>% 
  ungroup() %>% 
  mutate(perc_res = 100*non_zero_res/total_samp) %>%
  mutate(phase = factor(phase, levels = c(0, 24), labels = c("Baseline", "24 months")))

amr_dat_final <- amr_data_summ %>% 
  dplyr::select(phase, grappe, ab_class, avg_res)

# Save outcome dataset ----------------------------------------------------

saveRDS(amr_data_summ, here("data", "output", "amr_dat_full.rds"))
saveRDS(amr_dat_final, file.path(path_internal, "amr_dat.rds"))
saveRDS(amr_indiv, file.path(path_internal, "amr_indiv.rds"))

saveRDS(amr_indiv_pub, file.path(path_data_clean, "amr_indiv.rds"))

# Plots to check data -----------------------------------------------------

ggplot(amr_data_summ) + geom_boxplot(aes(y = avg_res, x = phase, group = phase)) +
  facet_wrap(~ab_class, scales = "free_y")

ggplot(amr_data_summ) + geom_boxplot(aes(y = perc_res, x = phase, group = phase)) +
  facet_wrap(~ab_class)

ggplot(amr_dat_long) + geom_boxplot(aes(y = resistance, x = grappe)) +
  facet_wrap(~phase+ab_class, scales = "free_y") +
  theme_minimal() + theme(axis.text.x = element_blank()) 


# NP data -----------------------------------------------------------------

np_raw <- read_csv(here("data", "untouched", "outcome", "MORDOR_0-12-24m_NP.csv"), show_col_types = FALSE)

np_dat <- np_raw %>%
  dplyr::select(
    Arm,
    Phase,
    WHG,
    contains("total sample"),
    contains("Erythr")) %>%
  drop_na(WHG) %>%
  rename(
    arm = Arm,
    phase = Phase,
    grappe = WHG,
    n_tested = `total sample tested#`,
    n_grew = `total sample grew#`,
    n_res = `#Resistant_Erythr`) %>%
  mutate(
    prop_res = n_res / n_grew,
    prop_res_test = n_res / n_tested) %>%
  mutate(phase = case_when(
    phase == "0" ~ "Baseline",
    phase == "12" ~ "12 months",
    phase == "24" ~ "24 months",
    TRUE ~ NA_character_
  )) %>%
  mutate(phase = factor(phase, levels = c("Baseline", "12 months", "24 months")))

saveRDS(np_dat, file.path(path_internal, "np_dat.rds"))

np_dat_pub <- np_dat %>%
  mutate(phase = as.character(phase)) %>%
  dplyr::filter(phase %in% c("Baseline", "24 months")) %>%
  mutate(phase = factor(phase, levels = c("Baseline", "24 months")))

saveRDS(np_dat, file.path(path_data_clean, "np_dat.rds"))
