########################################
##    FIRST DESCRIPTIVE ANALYSES
########################################

rm(list = ls())
library(tidyverse)
library(readxl)
library(gtsummary)
library(gt)
library(ggnetwork)
library(data.table)
library(ggpubr)
library(pscl)
library(stringr)
source("R/dictionaries.R")
source("R/helperFunctions.R")

#################################################
## Load data ------------------------------------
#################################################
microbio = read.delim("data/AMR-COVID-Data-230615/AMR_COVID_MICROBIO_20230615_094106.txt", 
                      sep = "#", comment.char="", header = T) %>%
  mutate(PVTDAT = as.Date(PVTDAT, "%d/%m/%Y"),
         PVTRES_IC = factor(PVTRES_IC, c(1,2), c("Infection", "Colonisation")))
microbio$PVTDAT[microbio$PVTDAT == as.Date("5202-05-05")] = as.Date("2022-05-05")

stays = read.table("data/AMR-COVID-Data-230615/AMR_COVID_SEJOURS_20230615_094432.txt", 
                   sep = "#", comment.char = "", header = T, encoding = 'latin1') %>%
  mutate(
    PATDDN = as.Date(PATDDN, "%d/%m/%Y"),
    REAENT = as.Date(REAENT, "%d/%m/%Y"),
    REASOR = as.Date(REASOR, "%d/%m/%Y"),
    PATINCD = as.Date(PATINCD, "%d/%m/%Y"),
    HOSENT = as.Date(HOSENT, "%d/%m/%Y"),
    REABCAD = as.Date(REABCAD, "%d/%m/%Y"),
    REASOR1 = as.Date(REASOR1, "%d/%m/%Y"),
    REASOR2 = as.Date(REASOR2, "%d/%m/%Y"),
    REASOR3 = as.Date(REASOR3, "%d/%m/%Y"),
    REASOR4 = as.Date(REASOR4, "%d/%m/%Y"),
    REASOR5 = as.Date(REASOR5, "%d/%m/%Y"),
    REAINT1D = as.Date(REAINT1D, "%d/%m/%Y"), 
    REAINT1F = as.Date(REAINT1F, "%d/%m/%Y"),
    REAINT2D = as.Date(REAINT2D, "%d/%m/%Y"),
    REAINT2F = as.Date(REAINT2F, "%d/%m/%Y"),
    REAINT3D = as.Date(REAINT3D, "%d/%m/%Y"),
    REAINT3F = as.Date(REAINT3F, "%d/%m/%Y"),
    REAINT4D = as.Date(REAINT4D, "%d/%m/%Y"),
    REAINT4F = as.Date(REAINT4F, "%d/%m/%Y"),
    ATBDEB1 = as.Date(ATBDEB1, "%d/%m/%Y"),
    ATBFIN1 = as.Date(ATBFIN1, "%d/%m/%Y"),
    ATBDEB2 = as.Date(ATBDEB2, "%d/%m/%Y"),
    ATBFIN2 = as.Date(ATBFIN2, "%d/%m/%Y"),
    ATBDEB3 = as.Date(ATBDEB3, "%d/%m/%Y"),
    ATBFIN3 = as.Date(ATBFIN3, "%d/%m/%Y"),
    ATBDEB4 = as.Date(ATBDEB4, "%d/%m/%Y"),
    ATBFIN4 = as.Date(ATBFIN4, "%d/%m/%Y"),
    ATBDEB5 = as.Date(ATBDEB5, "%d/%m/%Y"),
    ATBFIN5 = as.Date(ATBFIN5, "%d/%m/%Y"),
    ATBDEB6 = as.Date(ATBDEB6, "%d/%m/%Y"),
    ATBFIN6 = as.Date(ATBFIN6, "%d/%m/%Y"),
    ATBDEB7 = as.Date(ATBDEB7, "%d/%m/%Y"),
    ATBFIN7 = as.Date(ATBFIN7, "%d/%m/%Y"),
    PATSEX = factor(PATSEX, c(1,2), c("male", "female")),
    stay_days = as.numeric(difftime(REASOR, REAENT, units = "days"))
    ) 

n_patient_days = stays %>%
  group_by(PATPER) %>%
  summarise(patient_days = sum(stay_days))

periods = sort(unique(stays$PATPER))
length(unique(microbio$SUBJID))

#################################################
## Data cleaning --------------------------------
#################################################
summary(stays)
# NAs in 
  # REASOR4
  # REASOR5
  # ATBFIN7

# Known infection in 2018?
stays$REABCAD[stays$SUBJID == "045-1"] = NA

# Entries and discharges are not properly ordered
sum(stays$REAENT > stays$REASOR, na.rm = T)
sum(stays$REASOR1 > stays$REASOR2, na.rm = T)
sum(stays$REASOR2 > stays$REASOR3, na.rm = T)

# Start and end of intubation are not in agreement
sum(stays$REAINT1D > stays$REAINT1F, na.rm = T)
sum(stays$REAINT1F > stays$REAINT2D, na.rm = T)
sum(stays$REAINT2D > stays$REAINT2F, na.rm = T)
sum(stays$REAINT2F > stays$REAINT3D, na.rm = T)
sum(stays$REAINT3D > stays$REAINT3F, na.rm = T)
sum(stays$REAINT3F > stays$REAINT4D, na.rm = T)
sum(stays$REAINT4D > stays$REAINT4F, na.rm = T)

# Start and end of antibiotic treatments are not in agreement
sum(stays$ATBDEB1 > stays$ATBFIN1, na.rm = T)
sum(!is.na(stays$ATBDEB1) & is.na(stays$ATBFIN1))
sum(stays$ATBDEB2 > stays$ATBFIN2, na.rm = T)
sum(!is.na(stays$ATBDEB2) & is.na(stays$ATBFIN2))
sum(stays$ATBDEB3 > stays$ATBFIN3, na.rm = T)
sum(!is.na(stays$ATBDEB3) & is.na(stays$ATBFIN3))
sum(stays$ATBDEB4 > stays$ATBFIN4, na.rm = T)
sum(!is.na(stays$ATBDEB4) & is.na(stays$ATBFIN4))
sum(stays$ATBDEB5 > stays$ATBFIN5, na.rm = T)
sum(!is.na(stays$ATBDEB5) & is.na(stays$ATBFIN5))
sum(stays$ATBDEB6 > stays$ATBFIN6, na.rm = T)
sum(!is.na(stays$ATBDEB6) & is.na(stays$ATBFIN6))
sum(stays$ATBDEB7 > stays$ATBFIN7, na.rm = T)
sum(!is.na(stays$ATBDEB7) & is.na(stays$ATBFIN7))

# Inclusion and patient agreement
table(stays$CI1, useNA = "always")
table(stays$CN1, useNA = "always")
nrow(stays)

# Bacterial sampling before or after patient stay in ICU
stays %>%
  select(SUBJID, PATPER, PATAGE, PATSEX, REAENT, REASOR) %>%
  right_join(., microbio, by = "SUBJID") %>%
  filter(PVTDAT > REASOR | PVTDAT < REAENT, !is.na(PVTDAT)) %>%
  select(SUBJID, PATPER) %>%
  distinct()

microbio %>%
  filter(SUBJID %in% c("004-1", "006-1", "007-1", "008-1", "014-1", "015-1", 
                       "018-1", "043-1", "044-1")) %>%
  group_by(SUBJID) %>%
  summarise(n = n())

# Number of species identified by sample
microbio$n_species = apply(microbio[, c("PVTKP", "PVTECO", "PVTECC", "PVTEAG", "PVTBA1", "PVTBA2", "PVTBA3")], 
      1, 
      function(x) {
  sum(x, na.rm = T)
})

table(microbio$n_species, microbio$PVTRES, 
      useNA = "always", dnn = c("Number of species", "Result"))
microbio %>%
  filter(PVTRES == 1, n_species == 0)

# Total number of patients in microbio data frame stratified by period
n_tot_patients = microbio %>% 
  select(SUBJID, PVTDAT) %>%  
  mutate(
    PATPER = case_when(
      PVTDAT <= as.Date("2021-01-01") ~ period_names[1],
      PVTDAT <= as.Date("2021-06-01") & PVTDAT >= as.Date("2021-02-01") ~ period_names[2],
      PVTDAT <= as.Date("2021-11-01") & PVTDAT >= as.Date("2021-06-01") ~ period_names[3],
      PVTDAT > as.Date("2022-01-01") ~ period_names[4]
    ),
  ) %>% 
  select(PATPER, SUBJID) %>% 
  distinct() %>% 
  group_by(PATPER) %>% 
  summarise(n = n())

#################################################
## Summary statistics on stays ------------------
#################################################
tab_list = list()
index = 1

for (i in periods) {
  tab_list[[index]] = stays %>% 
    filter(PATPER == i, !is.na(PATSEX)) %>%
    select(PATSEX, stay_days, PATAGE, REABCA, REAINT, ATB, ATB48H, PATDEC) %>%
    gtsummary::tbl_summary(
      by = PATSEX,
      list(
        PATAGE ~ "Age",
        stay_days ~ "Stay duration (in days)",
        REABCA ~ "Infection pre-ICU",
        REAINT ~ "Intubation",
        ATB ~ "Antibiotics at admission",
        ATB48H ~ "Antibiotics > 48h",
        PATDEC ~ "Death"
      )
    ) %>%
    bold_labels() %>%
    modify_header(label ~ "") %>%
    add_overall(last = TRUE) %>%
    add_p()
  
  index = index + 1
}
to_save = tbl_merge(
  tbls = tab_list,
  tab_spanner = period_names_bold[periods]
)
gtsave(as_gt(to_save), path = "tables/AMR-COVID-Data-230615/", filename = "stay_description_periodsex.png")


to_save = stays %>% 
  select(PATPER, PATSEX, stay_days, PATAGE, REABCA, REAINT, ATB, ATB48H, PATDEC) %>%
  mutate(PATPER = recode(PATPER, !!!period_names_bold)) %>%
  gtsummary::tbl_summary(
    by = PATPER,
    list(
      PATSEX ~ "Gender",
      PATAGE ~ "Age",
      stay_days ~ "Stay duration (in days)",
      REABCA ~ "Infection pre-ICU",
      REAINT ~ "Intubation",
      ATB ~ "Antibiotics at admission",
      ATB48H ~ "Antibiotics > 48h",
      PATDEC ~ "Death"
    )
  ) %>%
  bold_labels() %>%
  modify_header(label ~ "") %>%
  add_overall(last = TRUE) %>%
  add_p()
gtsave(as_gt(to_save), path = "tables/AMR-COVID-Data-230615/", filename = "stay_description_period.png")


p1 = stays %>%
  mutate(PATPER = recode(PATPER, !!!period_names)) %>%
  ggplot(., aes(x = PATSEX)) +
  facet_grid(cols = vars(PATPER)) +
  geom_bar() +
  theme_classic() +
  labs(x = "", y = "Count", title = "Gender")
  
p2 = stays %>%
  mutate(PATPER = recode(PATPER, !!!period_names)) %>%
  ggplot(., aes(x = PATAGE)) +
  facet_grid(cols = vars(PATPER)) +
  geom_histogram() + 
  theme_classic() +
  labs(x = "", title = "Patient age", y = "Count")

p3 = stays %>%
  mutate(PATPER = recode(PATPER, !!!period_names)) %>%
  ggplot(., aes(x = stay_days)) +
  facet_grid(cols = vars(PATPER)) +
  geom_histogram() + 
  theme_classic() +
  labs(x = "", title = "Stay length (in days)", y = "Count")

to_save = ggarrange(p1,p2,p3, nrow = 3, ncol = 1)
ggsave("figures/AMR-COVID-Data-230615/patient_characteristic_description.png", 
       to_save, height = 10, width = 8)

## Plot stays ------------------------------------
# Period 1
stays %>%
  mutate(
    INFCOV = factor(INFCOV, c(0,1), c("non", "oui")),
    PATDEC = factor(PATDEC, c(0,1), c("non", "oui")),
    ATB = factor(ATB, c(0,1), c("non", "oui")),
    ATB48H = factor(ATB48H, c(0,1), c("non", "oui")),
  ) %>%
  filter(REAENT < as.Date("2021-01-01")) %>%
  ggplot(.) +
  geom_rect(aes(xmin = as.Date("2020-01-01"), xmax = as.Date("2020-01-31"), ymin = -Inf, ymax = Inf), fill = "grey90") +
  geom_point(aes(x = REASOR, y = SUBJID, shape = "death",  color = PATDEC)) +
  geom_segment(aes(x = REAENT, xend = REASOR, y = SUBJID, yend = SUBJID, color = PATSEX)) +
  geom_point(aes(x = HOSENT, y = SUBJID, shape = "hospitalization", color = PATSEX)) +
  geom_point(aes(x = REABCAD, y = SUBJID, shape = "pre-ICU infection",  color = PATSEX)) +
  geom_point(aes(x = REAENT, y = SUBJID, shape = "antibiotics",  color = ATB)) +
  geom_point(aes(x = REAENT + days(2), y = SUBJID, shape = "antibiotics",  color = ATB48H)) +
  geom_rug(aes(y = SUBJID, color = INFCOV), sides = "r", linewidth = 1) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 16), 
        axis.text.x = element_text(size = 14), 
        legend.text = element_text(size = 14)) +
  scale_color_manual(name = "",
                     values = c("orange", "cornflowerblue", "#FFFFFF00", "red"),
                     labels = c("female", "male", "non", "oui"), 
                     na.value = "grey70") +
  labs(x = "", y = "Patients ID") +
  scale_shape_manual(name = "", values = c("hospitalization" = 2, 
                                           "death" = 5,
                                           "antibiotics" = 3,
                                           "pre-ICU infection" = 20))
ggsave("figures/AMR-COVID-Data-230615/stays_period1.png", height = 10, width = 8)

# Period 4
stays %>%
  mutate(
    INFCOV = factor(INFCOV, c(0,1), c("non", "oui")),
    PATDEC = factor(PATDEC, c(0,1), c("non", "oui")),
    ATB = factor(ATB, c(0,1), c("non", "oui")),
    ATB48H = factor(ATB48H, c(0,1), c("non", "oui")),
  ) %>%
  filter(REAENT > as.Date("2021-01-01")) %>%
  ggplot(.) +
  geom_rect(aes(xmin = as.Date("2022-05-01"), xmax = as.Date("2022-05-31"), ymin = -Inf, ymax = Inf), fill = "grey90") +
  geom_point(aes(x = REASOR, y = SUBJID, shape = "death",  color = PATDEC)) +
  geom_point(aes(x = HOSENT, y = SUBJID, shape = "hospitalization", color = PATSEX)) +
  geom_point(aes(x = REAENT, y = SUBJID, shape = "antibiotics",  color = ATB)) +
  geom_segment(aes(x = REAENT, xend = REASOR, y = SUBJID, yend = SUBJID, color = PATSEX)) +
  geom_point(aes(x = REABCAD, y = SUBJID, shape = "pre-ICU infection",  color = PATSEX)) +
  geom_point(aes(x = REAENT + days(2), y = SUBJID, shape = "antibiotics",  color = ATB48H)) +
  geom_rug(aes(y = SUBJID, color = INFCOV), sides = "r", linewidth = 1) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 16), 
        axis.text.x = element_text(size = 14), 
        legend.text = element_text(size = 14)) +
  scale_color_manual(name = "",
                     values = c("orange", "cornflowerblue", "#FFFFFF00", "red"),
                     labels = c("female", "male", "non", "oui"), 
                     na.value = "grey70") +
  labs(x = "", y = "Patients ID") +
  scale_shape_manual(name = "", values = c("hospitalization" = 2, 
                                           "death" = 5,
                                           "antibiotics" = 3,
                                           "pre-ICU infection" = 20))
ggsave("figures/AMR-COVID-Data-230615/stays_period4.png", height = 10, width = 8)

## Summary statistics infections -----------------
# Table of samples for period 1 stratified by sex
to_save = stays %>% 
  filter(PATPER == 1) %>%
  select(SUBJID, PATPER, PATAGE, PATSEX, REAENT, REASOR) %>%
  left_join(., microbio, by = "SUBJID") %>%
  filter(PVTDAT <= REASOR, PVTDAT >= REAENT) %>%
  group_by(SUBJID, PATPER, PATAGE, PATSEX) %>%
  summarise(
    N_samples = n(),
    N_pos_colo = sum(PVTRES_IC %in% "Colonisation" & PVTRES %in% 1),
    N_pos_inf = sum(PVTRES_IC %in% "Infection" & PVTRES %in% 1),
    N_pos_unknown = sum(is.na(PVTRES_IC) & PVTRES %in% 1),
    PVTKP = sum(PVTKP, na.rm = T),
    PVTECO = sum(PVTECO, na.rm = T),
    PVTECC = sum(PVTECC, na.rm = T),
    PVTEAG = sum(PVTEAG, na.rm = T),
    PVTBA1 = sum(PVTBA1, na.rm = T),
    .groups = "drop"
  ) %>%
  select(-c(SUBJID, PATPER)) %>%
  gtsummary::tbl_summary(
  by = PATSEX,
  type = list(all_continuous() ~ "continuous2",
              all_categorical() ~ "continuous2"),
  label = list(
    PATAGE ~ "Age",
    N_samples ~ "Total number of samples",
    N_pos_inf ~ "Positive infection samples", 
    N_pos_colo ~ "Positive colonisation samples", 
    N_pos_unknown ~ "Positive colonisation or infection samples",
    PVTKP ~ "K pneumoniae (ESBL)",
    PVTECO ~ "E coli (ESBL)",
    PVTECC ~ "E cloacae complex (ESBL)",
    PVTEAG ~ "E aerogenes (ESBL)",
    PVTBA1 ~ "Other 1 (ESBL)"
  )
) %>%
  bold_labels() %>%
  modify_header(label ~ "**Period 1**") %>%
  add_overall(last = TRUE)
gtsave(as_gt(to_save), path = "tables/AMR-COVID-Data-230615/", filename = "microbio_description.png")

# Bar plot of samples for period 1 stratified by sex and colonisation/infection status
positive_period1 = stays %>% 
  filter(PATPER == 1) %>%
  select(SUBJID, PATPER, PATAGE, PATSEX, REAENT, REASOR) %>%
  left_join(., microbio, by = "SUBJID") %>%
  filter(PVTDAT <= REASOR, PVTDAT >= REAENT) %>%
  pivot_longer(c(PVTKP, PVTECO, PVTECC, PVTEAG, PVTBA1), values_to = "result", names_to = "BACTERIA") %>%
  group_by(SUBJID, PATPER, PATAGE, PATSEX, PVTRES_IC, BACTERIA) %>%
  summarise(
    positive_sample = sum(result, na.rm = T),
    .groups = "drop"
  ) %>%
  mutate(BACTERIA = recode(BACTERIA, !!!bacteria_species)) %>%
  ggplot(., aes(x = positive_sample, fill = PVTRES_IC, col = PVTRES_IC)) +
  geom_bar() +
  theme_classic() +
  facet_grid(cols = vars(interaction(PATSEX, PVTRES_IC)), rows = vars(BACTERIA)) +
  labs(x = "Number of positive samples of ESBL-producing\nbacteria per patient during ICU stay",
       y = "Count", 
       fill = "", col = "")
ggsave("figures/AMR-COVID-Data-230615/microbio_description_sex_p1.png", positive_period1, height = 6, width = 8.5)

# Table of samples for periods 1, 3 and 4  
to_save = microbio %>%
  mutate(
    PATPER = case_when(
    PVTDAT <= as.Date("2021-01-01") ~ period_names_bold[1],
    PVTDAT <= as.Date("2021-06-01") & PVTDAT >= as.Date("2021-02-01") ~ period_names_bold[2],
    PVTDAT <= as.Date("2021-11-01") & PVTDAT >= as.Date("2021-06-01") ~ period_names_bold[3],
    PVTDAT > as.Date("2022-01-01") ~ period_names_bold[4]
  ),
  ) %>%
  group_by(SUBJID, PATPER) %>%
  summarise(
    N_samples = n(),
    N_pos_colo = sum(PVTRES_IC %in% "Colonisation" & PVTRES %in% 1),
    N_pos_inf = sum(PVTRES_IC %in% "Infection" & PVTRES %in% 1),
    N_pos_unknown = sum(is.na(PVTRES_IC) & PVTRES %in% 1),
    PVTKP = sum(PVTKP, na.rm = T),
    PVTECO = sum(PVTECO, na.rm = T),
    PVTECC = sum(PVTECC, na.rm = T),
    PVTEAG = sum(PVTEAG, na.rm = T),
    PVTBA1 = sum(PVTBA1, na.rm = T),
    .groups = "drop"
  ) %>%
  select(-SUBJID) %>%
  gtsummary::tbl_summary(
    by = PATPER,
    type = list(all_continuous() ~ "continuous2",
                all_categorical() ~ "continuous2"),
    statistic = all_continuous() ~ "{median} ({min},{max})",
    label = list(
      N_samples ~ "Total number of samples",
      N_pos_inf ~ "Positive infection samples", 
      N_pos_colo ~ "Positive colonisation samples", 
      N_pos_unknown ~ "Positive colonisation or infection samples",
      PVTKP ~ "K pneumoniae (ESBL)",
      PVTECO ~ "E coli (ESBL)",
      PVTECC ~ "E cloacae complex (ESBL)",
      PVTEAG ~ "E aerogenes (ESBL)",
      PVTBA1 ~ "Other 1 (ESBL)"
    )
  ) %>%
  bold_labels() %>%
  modify_header(label ~ "**Characteristics**") %>%
  add_overall(last = TRUE) %>%
  add_p()
gtsave(as_gt(to_save), path = "tables/AMR-COVID-Data-230615/", filename = "microbio_description_p1-3-4.png")

# Time interval between two consecutive tests
time_int = microbio %>%
  mutate(
    PATPER = case_when(
      PVTDAT <= as.Date("2021-01-01") ~ period_names[1],
      PVTDAT <= as.Date("2021-06-01") & PVTDAT >= as.Date("2021-02-01") ~ period_names[2],
      PVTDAT <= as.Date("2021-11-01") & PVTDAT >= as.Date("2021-06-01") ~ period_names[3],
      PVTDAT > as.Date("2022-01-01") ~ period_names[4]
    )
  ) %>%
  arrange(PATPER, SUBJID, PVTRES_IC, PVTDAT) %>%
  group_by(PATPER, SUBJID, PVTRES_IC) %>%
  summarise(test_interval = mean(as.numeric(difftime(PVTDAT, lag(PVTDAT), units = "day")), na.rm = T),
            .groups = "drop") %>%
  filter(!is.na(test_interval)) %>%
  ggplot(., aes(x = test_interval)) +
  geom_histogram() +
  facet_grid(cols = vars(PATPER), rows = vars(PVTRES_IC)) +
  theme_classic() +
  labs(x = "Average time interval between two microbiological tests (in days) per individual", 
       y = "Count")
ggsave("figures/AMR-COVID-Data-230615/time_interval_tests.png", time_int, height = 6, width = 6)

# Bar plot of samples for periods 1, 3 and 4  
positive_samples = microbio %>%
  mutate(
    PATPER = case_when(
      PVTDAT <= as.Date("2021-01-01") ~ period_names[1],
      PVTDAT <= as.Date("2021-06-01") & PVTDAT >= as.Date("2021-02-01") ~ period_names[2],
      PVTDAT <= as.Date("2021-11-01") & PVTDAT >= as.Date("2021-06-01") ~ period_names[3],
      PVTDAT > as.Date("2022-01-01") ~ period_names[4]
    ),
  ) %>%
  group_by(SUBJID, PATPER) %>%
  summarise(
    PVTKP = sum(PVTKP, na.rm = T),
    PVTECO = sum(PVTECO, na.rm = T),
    PVTECC = sum(PVTECC, na.rm = T),
    PVTEAG = sum(PVTEAG, na.rm = T),
    PVTBA1 = sum(PVTBA1, na.rm = T),
    .groups = "drop"
  ) %>%
  pivot_longer(c(PVTKP, PVTECO, PVTECC, PVTEAG, PVTBA1), names_to = "BACTERIA",values_to = "n_pos") %>%
  mutate(BACTERIA = recode(BACTERIA, !!!bacteria_species)) %>%
  ggplot(., aes(x = n_pos)) +
  geom_bar() +
  facet_grid(cols = vars(PATPER), rows = vars(BACTERIA)) +
  theme_classic() +
  labs(x = "Number of positive samples of ESBL-producing\nbacteria per patient during ICU stay", 
       y = "Count")
ggsave("figures/AMR-COVID-Data-230615/microbio_description_per_period.png", positive_samples, height = 6, width = 6)

#################################################
## Information on episodes ----------------------
#################################################
# Individuals with known infection pre-ICU admission
table(stays$REABCA, stays$PATPER, useNA = "always", dnn = c("Pre-ICU infection", "Period")) 

stays %>%
  filter(REABCA %in% 1) %>%
  select(SUBJID, PATSEX, PATAGE, REAENT, REASOR, PATPER, REABCA, REABCAE, 
         REABCAD, stay_days) %>%
  left_join(., microbio, by = "SUBJID") %>%
  mutate(REABCAE = ifelse(REABCAE == " ", NA, REABCAE),
         PVTRES = factor(PVTRES, labels = c("non", "oui"))) %>%
  ggplot(.) + 
  geom_point(aes(y = SUBJID, x = PVTDAT, col = PVTRES), shape = 3) +
  geom_point(aes(y = SUBJID, x = REABCAD, col = REABCAE), shape = 20) +
  geom_segment(aes(y = SUBJID, yend = SUBJID, x = REAENT, xend = REASOR, col = PATSEX)) +  scale_color_manual(labels = c("ESCHERICHIA COLI", "female", "KLEBSIELLA PNEUMONIAE", "male", "non", "oui"),
                     values = c("deeppink", "orange", "green2", "cornflowerblue", "grey80", "black")) +
  theme_classic() +
  labs(x = "Stay in ICU",  y = "Stay ID")
ggsave("figures/AMR-COVID-Data-230615/infection_pre_admission.png", height = 5, width = 8)

# Number of infections/colonisation during stay per bacterial species
n_pos_patient = microbio %>%
  mutate(
    PATPER = case_when(
      PVTDAT <= as.Date("2021-01-01") ~ period_names[1],
      PVTDAT <= as.Date("2021-06-01") & PVTDAT >= as.Date("2021-02-01") ~ period_names[2],
      PVTDAT <= as.Date("2021-11-01") & PVTDAT >= as.Date("2021-06-01") ~ period_names[3],
      PVTDAT > as.Date("2022-01-01") ~ period_names[4]
    ),
  ) %>%
  arrange(PATPER, SUBJID, PVTDAT) %>%
  group_by(PATPER, SUBJID) %>%
  summarise(
    PVTKP = any(PVTKP %in% 1),
    PVTECO = any(PVTECO %in% 1),
    PVTECC = any(PVTECC %in% 1),
    PVTEAG = any(PVTEAG %in% 1),
    PVTBA1 = any(PVTBA1 %in% 1),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = contains("PVT"), names_to = "BACTERIA", values_to = "value") %>%
  mutate(BACTERIA = recode(BACTERIA, !!!bacteria_species)) %>%
  group_by(PATPER, BACTERIA) %>%
  summarise(value = sum(value), .groups = "drop") %>%
  complete(PATPER, BACTERIA, fill = list(value = 0)) %>%
  ggplot(., aes(x = BACTERIA, y = value, fill = BACTERIA)) +
  facet_grid(cols = vars(PATPER)) +
  geom_bar(stat = "identity", width = 0.5) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), 
        legend.position = "none"
        ) +
  labs(x = "ESBL-producing bacteria", y = "Number of positive patients")
ggsave("figures/AMR-COVID-Data-230615/n_pos_patients.png", n_pos_patient, height = 6, width = 5)

# Table of the number of infections/colonisation during stay per bacterial species
to_save= microbio %>%
  mutate(
    PATPER = case_when(
      PVTDAT <= as.Date("2021-01-01") ~ period_names[1],
      PVTDAT <= as.Date("2021-06-01") & PVTDAT >= as.Date("2021-02-01") ~ period_names[2],
      PVTDAT <= as.Date("2021-11-01") & PVTDAT >= as.Date("2021-06-01") ~ period_names[3],
      PVTDAT > as.Date("2022-01-01") ~ period_names[4]
    ),
  ) %>%
  arrange(PATPER, SUBJID, PVTDAT) %>%
  group_by(PATPER, SUBJID) %>%
  summarise(
    PVTKP = any(PVTKP %in% 1),
    PVTECO = any(PVTECO %in% 1),
    PVTECC = any(PVTECC %in% 1),
    PVTEAG = any(PVTEAG %in% 1),
    PVTBA1 = any(PVTBA1 %in% 1),
    .groups = "drop"
  ) %>%
  select(-SUBJID) %>%
  gtsummary::tbl_summary(
    by = PATPER,
    label = list(PVTKP ~ "K pneumoniae",
                 PVTECO ~ "E coli",
                 PVTECC ~ "E cloacae",
                 PVTEAG ~ "E aerogenes",
                 PVTBA1 ~ "Other"
                 )
  ) %>%
  bold_labels() %>%
  modify_header(label ~ "**ESBL-producing bacteria**") %>%
  add_overall(last = TRUE) %>%
  add_p()
gtsave(as_gt(to_save), path = "tables/AMR-COVID-Data-230615/", filename = "microbio_pos_patient_p1-3-4.png")

# Type of episodes
episodes = microbio %>%
  mutate(
    PATPER = case_when(
      PVTDAT <= as.Date("2021-01-01") ~ period_names[1],
      PVTDAT <= as.Date("2021-06-01") & PVTDAT >= as.Date("2021-02-01") ~ period_names[2],
      PVTDAT <= as.Date("2021-11-01") & PVTDAT >= as.Date("2021-06-01") ~ period_names[3],
      PVTDAT > as.Date("2022-01-01") ~ period_names[4]
    ),
  ) %>%
  arrange(PATPER, SUBJID, PVTNAT, PVTDAT) %>%
  select(PATPER, SUBJID, PVTNAT, PVTDAT, PVTKP, PVTECC, PVTECO, PVTEAG, PVTBA1) %>%
  pivot_longer(c(PVTKP, PVTECC, PVTECO, PVTEAG, PVTBA1), names_to = "BACTERIA", values_to = "value") %>%
  group_by(PATPER, SUBJID, PVTNAT, BACTERIA) %>%
  nest() %>%
  mutate(episode = map(data, function(.data) getEpisodeCategory(.data))) %>%
  select(-data) %>%
  unnest(cols = episode) %>%
  ungroup() %>%
  mutate(BACTERIA = recode(BACTERIA, !!!bacteria_species))

tab_list = list()
index = 1
for (b in unique(episodes$BACTERIA)) {
  tab_list[[index]] = 
  episodes %>% 
    filter(BACTERIA == b) %>%
    select(-c(SUBJID, BACTERIA, PVTNAT)) %>%
    gtsummary::tbl_summary(
      by = PATPER,
      label = list(episode ~ b)
    ) %>%
    bold_labels() %>%
    modify_header(label ~ "**ESBL-producing bacteria**") %>%
    add_overall(last = TRUE) %>%
    add_p()
  index = index + 1
}

to_save = tbl_stack(tbls = tab_list)
gtsave(as_gt(to_save), path = "tables/AMR-COVID-Data-230615/", filename = "microbio_episode_category.png")

# Sample at the admission for all patients ?
summary(microbio %>%
  inner_join(., stays, by = "SUBJID") %>%
  mutate(sample_admission = as.numeric(difftime(PVTDAT, REAENT, units = "day"))) %>%
  group_by(SUBJID) %>%
  summarise(sample_admission = any(sample_admission %in% 0)))

to_save=microbio %>%
  inner_join(., stays, by = "SUBJID") %>%
  mutate(sample_admission = as.numeric(difftime(PVTDAT, REAENT, units = "day"))) %>%
  group_by(SUBJID) %>%
  summarise(sample_admission = min(sample_admission)) %>%
  ggplot(., aes(x = sample_admission)) +
  geom_histogram() +
  theme_classic() +
  labs(x = "Time from admission to first sample (period 1)", y = "Count")
ggsave("figures/AMR-COVID-Data-230615/time_to_first_sample.png", to_save, height = 3, width = 4)

#################################################
## Plot individual stays ------------------------
#################################################
stays %>%
  mutate(
    PATSEX = factor(PATSEX, c(1,2), c("male", "female")),
    INFCOV = factor(INFCOV, c(0,1), c("non", "oui")),
    PATDEC = factor(PATDEC, c(0,1), c("non", "oui")),
    ATB = factor(ATB, c(0,1), c("non", "oui")),
    ATB48H = factor(ATB48H, c(0,1), c("non", "oui")),
  ) %>%
  filter(REAENT < as.Date("2021-01-01")) %>%
  ggplot(.) +
  geom_point(aes(x = REASOR, y = SUBJID, shape = "death",  color = PATDEC)) +
  geom_segment(aes(x = REAENT, xend = REASOR, y = SUBJID, yend = SUBJID, color = PATSEX)) +
  geom_point(aes(x = HOSENT, y = SUBJID, shape = "hospitalization", color = PATSEX)) +
  geom_point(aes(x = REABCAD, y = SUBJID, shape = "pre-ICU infection",  color = PATSEX)) +
  geom_point(aes(x = REAENT, y = SUBJID, shape = "antibiotics",  color = ATB)) +
  geom_point(aes(x = REAENT + days(2), y = SUBJID, shape = "antibiotics",  color = ATB48H)) +
  geom_rug(aes(y = SUBJID, color = INFCOV), sides = "r", linewidth = 1) +
  # geom_polygon(data = infectionData, aes(x = -25, y = identifiant, group = factor(hhid), fill = "household"),
  #             inherit.aes = F, color = "black", size = 1.2) +
  theme_minimal() +
  theme(#axis.text.y = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.title = element_text(size = 16), 
    axis.text.x = element_text(size = 14), 
    legend.text = element_text(size = 14)) +
  scale_color_manual(name = "",
                     values = c("orange", "cornflowerblue", "#FFFFFF", "red"),
                     labels = c("female", "male", "non", "oui"), 
                     na.value = "grey70") +
  labs(x = "", y = "Patients ID") +
  scale_shape_manual(name = "", values = c("hospitalization" = 2, 
                                           "death" = 5,
                                           "antibiotics" = 3,
                                           "pre-ICU infection" = 20))

#################################################
## Contact network ------------------------------
#################################################
## Measures of disorganization ------------------
network_tab = stays %>%
  rename(REAENT1 = REAENT) %>%
  select(PATPER, SUBJID, contains("REAENT"), contains("REASEC"), contains("REACHA"), matches("REASOR[0-9]")) %>%
  mutate(
    REAENT2 = as.Date(ifelse(!is.na(REASOR2), as.character(REASOR1), NA)), 
    REAENT3 = as.Date(ifelse(!is.na(REASOR3), as.character(REASOR2), NA)), 
    REAENT4 = as.Date(ifelse(!is.na(REASOR4), as.character(REASOR3), NA)), 
    REAENT5 = as.Date(ifelse(!is.na(REASOR5), as.character(REASOR4), NA))
  ) %>%
  pivot_longer(-c(PATPER, SUBJID), 
               names_pattern = "(.*)(.)", 
               names_to = c(".value", "STAYID")) %>%
  mutate(STAYID = as.numeric(STAYID)) %>%
  filter(!is.na(REASEC), !REASEC %in% " ") %>%
  arrange(SUBJID, REAENT)

# Diagram of stays per sector and room  
network_tab %>%
  filter(PATPER == 1) %>%
  ggplot(.) +
  geom_segment(aes(x = REAENT, xend = REASOR, y = SUBJID, yend = SUBJID, color = factor(REACHA))) +
  facet_grid(cols = vars(REASEC)) +
  theme_minimal() +
  labs(x = "ICU stay", y = "Patient stay ID", col = "Room")

# No. of room sharing events per periods
room_sector_tab = network_tab %>%
  select(REASEC, REACHA, PATPER) %>%
  distinct()

room_sector_tab$n_sharing = apply(room_sector_tab, 1, function(x) {
  tab_to_test = network_tab %>%
    filter(REASEC == x[["REASEC"]], REACHA == as.numeric(x[["REACHA"]]), PATPER == x[["PATPER"]]) %>%
    select(REAENT, REASOR) %>%
    mutate(REAENT = as.character(REAENT), REASOR = as.character(REASOR))
  
  if(nrow(tab_to_test) <= 1) {
    return(0)
  } else {
    rowids = combn(1:nrow(tab_to_test), 2)
    comb = cbind(tab_to_test[rowids[1,], ],
                 tab_to_test[rowids[2,], ])
    colnames(comb) = c("REAENT1", "REASOR1", "REAENT2", "REASOR2")
    out = sum(
      (comb$REAENT1 >= comb$REAENT2 & comb$REASOR1 <= comb$REASOR2) |
        (comb$REAENT2 >= comb$REAENT1 & comb$REASOR2 <= comb$REASOR1)
      )
    return(out)
  }
})
room_sector_tab$n_sharing[room_sector_tab$REACHA == 4 & 
                            room_sector_tab$REASEC == "SELF" &
                            room_sector_tab$PATPER == 1] = 0 

room_sector_tab %>%
  group_by(PATPER) %>%
  summarise(n_sharing = sum(n_sharing), .groups = "drop") %>%
  left_join(., n_patient_days, by = "PATPER") %>%
  mutate(prev = n_sharing / patient_days * 1000) %>%
  write.table(., "tables/AMR-COVID-Data-230615/room_sharing.txt", quote = F, row.names = F, sep = "\t")

# network_tab %>%
#   filter(PATPER == 1, REASEC == "SELF", REACHA %in% c(4,6)) %>%
#   arrange(REACHA, REAENT, REASOR)
  
# No. of transfers 
transfers = lapply(stays$SUBJID, function(x) {
  network_tab %>%
    filter(SUBJID %in% x) %>%
    arrange(REAENT) %>%
    group_by(PATPER, SUBJID) %>%
    summarise(
      n_transfers = n()-1,
      n_same_sector_diff_room = sum(REACHA != lag(REACHA) & REASEC == lag(REASEC), na.rm = T),
      n_diff_sector = sum(REASEC != lag(REASEC), na.rm = T),
      .groups = "drop"
        )
}) 

do.call("rbind", transfers) %>%
  group_by(PATPER) %>%
  summarise(
    n_transfers = sum(n_transfers),
    n_same_sector_diff_room = sum(n_same_sector_diff_room),
    n_diff_sector = sum(n_diff_sector)
    ) %>%
  left_join(., n_patient_days, by = "PATPER") %>%
  mutate(
    n_transfers = paste0(n_transfers, " (", round(n_transfers/patient_days*1000,2),")"),
    n_same_sector_diff_room = paste0(n_same_sector_diff_room, " (", round(n_same_sector_diff_room/patient_days*1000,2),")"),
    n_diff_sector = paste0(n_diff_sector, " (", round(n_diff_sector/patient_days*1000,2),")")
  ) %>%
  write.table(., "tables/AMR-COVID-Data-230615/transfers.txt", quote = F, row.names = F, sep = "\t")


#################################################
## Antibiotic exposure and infection/colonisation
#################################################
tab = stays %>%
  select(PATPER, SUBJID, PATAGE, PATSEX, HOSTRA, REATRA, REAINT, stay_days, REAENT, REASOR, ATB, ATB48H) %>%
  inner_join(., microbio, by = "SUBJID") %>%
  filter(PVTDAT >= REAENT+2 & PVTDAT <= REASOR) %>%
  group_by(PATPER, SUBJID, PATAGE, PATSEX, stay_days, REAINT, HOSTRA, REATRA, ATB, ATB48H) %>%
  summarise(
    PVTKP = as.numeric(any(PVTKP %in% 1)),
    PVTECO = as.numeric(any(PVTECO %in% 1)),
    PVTECC = as.numeric(any(PVTECC %in% 1)),
    PVTEAG = as.numeric(any(PVTEAG %in% 1)),
    PVTBA1 = as.numeric(any(PVTBA1 %in% 1)),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(PVTKP, PVTECO, PVTECC, PVTEAG, PVTBA1), 
               names_to = "BACTERIA", values_to = "RES") 

mod = glm(RES ~ PATAGE + PATSEX + stay_days + HOSTRA + REATRA + REAINT + ATB + ATB48H, 
          family = binomial(link = "logit"), 
          data = tab[tab$BACTERIA %in% "PVTKP", ])
summary(mod)
anova(mod)
qqnorm(residuals(mod))

mod2 = glm(RES ~ PATAGE + PATSEX + stay_days + REAINT + ATB + ATB48H, 
          family = binomial(link = "logit"), 
          data = tab[tab$BACTERIA %in% "PVTKP", ])
summary(mod2)
qqnorm(residuals(mod2))

mod3 = glm(RES ~ PATAGE + PATSEX + stay_days + ATB + ATB48H, 
           family = binomial(link = "logit"), 
           data = tab[tab$BACTERIA %in% "PVTKP", ])
summary(mod3)
anova(mod3)
pR2(mod3)
qqnorm(residuals(mod3, type = "deviance"))

mod4 = glm(RES ~ PATAGE + PATSEX + ATB + ATB48H, 
           family = binomial(link = "logit"), 
           data = tab[tab$BACTERIA %in% "PVTKP", ])
summary(mod4)
anova(mod4)
pR2(mod4)
qqnorm(residuals(mod4, type = "deviance"))



