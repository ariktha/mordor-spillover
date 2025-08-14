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

mls_dat <- amr_dat_final %>% 
  dplyr::filter(ab_class == "MLS") %>%
  dplyr::select(-ab_class)


# Save outcome dataset ----------------------------------------------------

saveRDS(amr_data_summ, here("data", "output", "amr_dat_full.rds"))
saveRDS(amr_dat_final, here("data", "clean", "amr_dat.rds"))
saveRDS(mls_dat, here("data", "output", "mls_dat.rds"))


# Plots to check data -----------------------------------------------------

ggplot(amr_data_summ) + geom_boxplot(aes(y = avg_res, x = phase, group = phase)) +
  facet_wrap(~ab_class, scales = "free_y")

ggplot(amr_data_summ) + geom_boxplot(aes(y = perc_res, x = phase, group = phase)) +
  facet_wrap(~ab_class)

ggplot(amr_dat_long) + geom_boxplot(aes(y = resistance, x = grappe)) +
  facet_wrap(~phase+ab_class, scales = "free_y") +
  theme_minimal() + theme(axis.text.x = element_blank()) 
