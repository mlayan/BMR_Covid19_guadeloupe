########################################
##    FIRST DESCRIPTIVE ANALYSES
########################################

rm(list = ls())
library(MASS)
library(tidyverse)
library(readxl)
library(gtsummary)
library(gt)
library(ggnetwork)
library(data.table)
library(ggpubr)
library(pscl)
library(stringr)
library(ResourceSelection)
source("R/dictionaries.R")
source("R/helperFunctions.R")

#################################################
## Load data ------------------------------------
#################################################
# Microbiological data
microbio = read.delim("data/AMR_COVID-Data-230803/AMRCOVID_MICROBIO_20230803_151902.txt", 
                      sep = "#", comment.char="", header = T, encoding = "latin1") %>%
  mutate(
    PVTDAT = as.Date(PVTDAT, "%d/%m/%Y"),
    PVTKP = ifelse(PVTKP %in% "OUI", 1, 0),
    PVTECO = ifelse(PVTECO %in% "OUI", 1, 0),
    PVTECC = ifelse(PVTECC %in% "OUI", 1, 0),
    PVTEAG = ifelse(PVTEAG %in% "OUI", 1, 0),
    PVTBA1 = ifelse(PVTBA1 %in% "OUI", 1, 0),
    PVTBA2 = ifelse(PVTBA2 %in% "OUI", 1, 0),
    PVTBA3 = ifelse(PVTBA3 %in% "OUI", 1, 0)
    )
summary(microbio)
sapply(microbio %>% select(-c(SUBJID, PVTID, PVTDAT)), table, useNA = "always")
table(microbio$PVTRES_IC[microbio$PVTRES == "Positif"], microbio$PVTNAT[microbio$PVTRES == "Positif"], useNA = "always")

# Stay data
stays = read.table("data/AMR_COVID-Data-230803/AMRCOVID_SEJOURS_20230803_160305_modchar.txt", 
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
    ATBDEB8 = as.Date(ATBDEB8, "%d/%m/%Y"),
    ATBFIN8 = as.Date(ATBFIN8, "%d/%m/%Y"),
    PATSTART = recode(PATPER, !!!setNames(periods_start, periods)),
    PATEND = recode(PATPER, !!!setNames(periods_end, periods)),
    stay_days = as.numeric(difftime(REASOR, REAENT, units = "days"))
    ) 

stays %>%
  group_by(PATPER) %>%
  summarise(
    n = n(), 
    n_no_dup = length(unique(SUBJID)), 
    stay_length = median(stay_days),
    n_males = sum(PATSEX %in% "Masculin"),
    n_females = sum(PATSEX %in% "Feminin"),
    p_males = sum(PATSEX %in% "Masculin") / n(),
    .groups = 'drop'
    )

n_patient_days = stays %>%
  group_by(PATPER) %>%
  summarise(patient_days = sum(stay_days))

length(unique(microbio$SUBJID))

# Quality check of the microbiological results
table(microbio$PVTRES, microbio$PVTKP, useNA = "always", dnn = c("Resultat test", "K pneumoniae"))
table(microbio$PVTRES, microbio$PVTECC, useNA = "always", dnn = c("Resultat test", "E cloacae complex"))
table(microbio$PVTRES, microbio$PVTECO, useNA = "always", dnn = c("Resultat test", "E coli"))
table(microbio$PVTRES, microbio$PVTEAG, useNA = "always", dnn = c("Resultat test", "E aerogenes"))
table(microbio$PVTRES, microbio$PVTBA1, useNA = "always", dnn = c("Resultat test", "Other 1"))
table(microbio$PVTRES, microbio$PVTBA2, useNA = "always", dnn = c("Resultat test", "Other 2"))
table(microbio$PVTRES, microbio$PVTBA3, useNA = "always", dnn = c("Resultat test", "Other 3"))

table(microbio$PVTRES, microbio$PVTRES_IC, useNA = "always")

at_least_one_positive = apply(microbio %>% select(PVTKP, PVTECC, PVTECO, PVTEAG), 1, sum)
table(microbio$PVTRES, at_least_one_positive, useNA = "always")

# Number of subjects for which rectal samples are available
microbio %>% 
  filter(PVTNAT %in% "Pvt rectal") %>%
  select(SUBJID) %>%
  distinct() %>%
  left_join(., stays %>% select(SUBJID, PATPER), by = "SUBJID") %>%
  group_by(PATPER) %>%
  summarise(n = n())

#################################################
## Data cleaning --------------------------------
#################################################
# NAs in 
  # ATBFIN7
  # ATB48HRI
  # ATB48HCY

# Known infection in 2018?
stays$REABCAD[stays$SUBJID == "045-1"] = NA

# Entries and discharges are not properly ordered
summary(stays %>% select(contains("REASOR")))
sum(stays$REAENT > stays$REASOR, na.rm = T)
sum(stays$REASOR1 > stays$REASOR2, na.rm = T)
sum(stays$REASOR2 > stays$REASOR3, na.rm = T)
sum(stays$REASOR3 > stays$REASOR4, na.rm = T)
sum(stays$REASOR5 > stays$REASOR5, na.rm = T)

# Start and end of intubation are not in agreement
summary(stays %>% select(contains("REAINT")))
sum(stays$REAINT1D > stays$REAINT1F, na.rm = T)
sum(stays$REAINT1F > stays$REAINT2D, na.rm = T)
sum(stays$REAINT2D > stays$REAINT2F, na.rm = T)
sum(stays$REAINT2F > stays$REAINT3D, na.rm = T)
sum(stays$REAINT3D > stays$REAINT3F, na.rm = T)
sum(stays$REAINT3F > stays$REAINT4D, na.rm = T)
sum(stays$REAINT4D > stays$REAINT4F, na.rm = T)

# Start and end of antibiotic treatments are not in agreement
summary(stays %>% select(contains("ATBDEB")))
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
  select(SUBJID, PATPER, PVTDAT, REAENT, REASOR) %>%
  distinct()

stays %>%
  select(SUBJID, PATPER, PATAGE, PATSEX, REAENT, REASOR) %>%
  right_join(., microbio, by = "SUBJID") %>%
  filter(SUBJID %in% c("015-1", "097-1", "107-1", "108-1", "130-1", "273-1", 
                       "281-1", "309-1", "364-1", "531-1", "541-1", "555-1", 
                       "558-1")) %>%
  group_by(SUBJID) %>%
  summarise(
    n_tot = n(), 
    n_before = sum(PVTDAT < REAENT & !is.na(PVTDAT)),
    n_after = sum(PVTDAT > REASOR & !is.na(PVTDAT))
    )

# Number of species identified by sample
microbio$n_species = apply(microbio[, c("PVTKP", "PVTECO", "PVTECC", "PVTEAG", "PVTBA1", "PVTBA2", "PVTBA3")], 
      1, 
      function(x) {
  sum(x, na.rm = T)
})

table(microbio$n_species, microbio$PVTRES, 
      useNA = "always", dnn = c("Number of species", "Result"))
microbio %>%
  filter(n_species > 1)

# Are patients with microbiological data in the stays database ?
length(unique(microbio$SUBJID))
sum(unique(microbio$SUBJID) %in% stays$SUBJID) 
nrow(stays)

# Number of individuals with no microbiological data
# the 3rd wave is not yet consolidated and this is the 
# period for which the number of missing microbiological
# test is the highest
sum(!stays$SUBJID %in% microbio$SUBJID)
stays %>%
  filter(!SUBJID %in% microbio$SUBJID) %>%
  count(PATPER)

#################################################
## Summary statistics on stays ------------------
#################################################
# Summary statistics for each period stratified by sex
for (i in c(1,3,4)) {
  to_save = stays %>% 
    filter(PATPER == periods[i], !is.na(PATSEX)) %>%
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
    modify_header(label ~ period_names_bold[[i]]) %>%
    add_overall(last = TRUE) %>%
    add_p()
  
  gtsave(as_gt(to_save), 
         path = "tables/AMR_COVID-Data-230803/", 
         filename = paste0("stay_description_periodsex_", period_names_short[i], ".png"))

}

# Comparison of the summary statistics for each period
to_save = stays %>% 
  select(PATPER, PATSEX, stay_days, PATAGE, REAINT, ATB, ATB48H, PATDEC) %>%
  filter(PATPER != periods[2]) %>%
  mutate(
    PATPER = factor(PATPER, levels = periods[c(1,3,4)], labels = period_names[c(1,3,4)]),
    REAINT = ifelse(REAINT == " ", NA, tolower(REAINT)), 
    ATB = ifelse(ATB == " ", NA, tolower(ATB)), 
    ATB48H = ifelse(ATB48H == " ", NA, tolower(ATB48H)), 
    PATDEC = ifelse(PATDEC == " ", NA, tolower(PATDEC))
  ) %>%
  gtsummary::tbl_summary(
    by = PATPER,
    list(
      PATSEX ~ "Gender",
      PATAGE ~ "Age",
      stay_days ~ "Stay length (days)",
      #REABCA ~ "Infection pre-ICU",
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
gtsave(as_gt(to_save), path = "tables/AMR_COVID-Data-230803/", 
       filename = "stay_description_period.png")


p1 = stays %>%
  filter(PATPER != periods[2]) %>%
  mutate(PATPER = factor(PATPER, levels = periods[c(1,3,4)], labels = period_names[c(1,3,4)])) %>%
  ggplot(., aes(x = PATSEX)) +
  facet_grid(cols = vars(PATPER)) +
  geom_bar() +
  theme_classic() +
  labs(x = "", y = "Count", title = "Gender")
  
p2 = stays %>%
  filter(PATPER != periods[2]) %>%
  mutate(PATPER = factor(PATPER, levels = periods[c(1,3,4)], labels = period_names[c(1,3,4)])) %>%
  ggplot(., aes(x = PATAGE)) +
  facet_grid(cols = vars(PATPER)) +
  geom_histogram() + 
  theme_classic() +
  labs(x = "", title = "Patient age", y = "Count")

p3 = stays %>%
  filter(PATPER != periods[2]) %>%
  mutate(PATPER = factor(PATPER, levels = periods[c(1,3,4)], labels = period_names[c(1,3,4)])) %>%
  ggplot(., aes(x = stay_days)) +
  facet_grid(cols = vars(PATPER)) +
  geom_histogram() + 
  theme_classic() +
  labs(x = "", title = "Stay length (in days)", y = "Count")

to_save = ggarrange(p1,p2,p3, nrow = 3, ncol = 1)
ggsave("figures/AMR_COVID-Data-230803/patient_characteristic_description.png", 
       to_save, height = 10, width = 8)

## Plot stays
for (p in c(1,3,4)) {
  to_save = stays %>%
    mutate(
      INFCOV = factor_oui_non(INFCOV), 
      PATDEC = factor_oui_non(PATDEC),
      ATB = factor_oui_non(ATB),
      ATB48H = factor_oui_non(ATB48H)
    ) %>%
    filter(PATPER == periods[p]) %>%
    arrange(REAENT) %>%
    ggplot(.) +
    geom_rect(aes(xmin = periods_start[p], xmax = periods_end[p], ymin = -Inf, ymax = Inf), fill = "grey90") +
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
                       labels = c("Feminin", "Masculin", "non", "oui"),
                       na.value = "grey70") +
    labs(x = "", y = "Patients ID") +
    scale_shape_manual(name = "", values = c("hospitalization" = 2, 
                                             "death" = 5,
                                             "antibiotics" = 3,
                                             "pre-ICU infection" = 20))
  ggsave(paste0("figures/AMR_COVID-Data-230803/stays_", period_names_short[p], ".png"), 
         to_save, height = 12, width = 8)
}

#################################################
## COVID-19 prevalence data ---------------------
#################################################
# Plot COVID-19 prevalence per period
covid_prev = list()
i= 1
for (p in periods[-2]){
  covid_prev[[i]] = stays %>%
    filter(PATPER == p) %>%
    select(SUBJID, PATPER, REAENT, PATSTART, PATEND, REASOR, INFCOV) %>%
    mutate(DATE = map2(REAENT, REASOR, seq, by = "1 day"),
           PATPER = gsub(" \\(", "\n\\(", PATPER)) %>%
    unnest(cols = DATE) %>%
    filter(DATE >= PATSTART, DATE <= PATEND) %>%
    group_by(PATPER, DATE) %>%
    summarise(p_cov = sum(INFCOV %in% "OUI") / n(), .groups = "drop") %>%
    ggplot(., aes(x = DATE, y = p_cov)) +
    geom_line() +
    facet_grid(cols = vars(PATPER)) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x= "", y = "COVID-19 patient prevalence") +
    ylim(c(0,1))
  i=i+1
}
covid_prev_save = ggarrange(plotlist = covid_prev, ncol = 3, nrow = 1)
ggsave("figures/AMR_COVID-Data-230803/covid_prevalence_period.png", 
       covid_prev_save, height = 4, width = 10)

# Save covid-19 prevalence per period
stays %>%
  filter(PATPER != periods[2]) %>%
  select(SUBJID, PATPER, REAENT, PATSTART, PATEND, REASOR, INFCOV) %>%
  mutate(DATE = map2(REAENT, REASOR, seq, by = "1 day")) %>%
  unnest(cols = DATE) %>%
  filter(DATE >= PATSTART, DATE <= PATEND) %>%
  group_by(PATPER, DATE) %>%
  summarise(p_cov = sum(INFCOV %in% "OUI") / n(), .groups = "drop") %>%
  write.table(., "data/AMR_COVID-Data-230803/covid_prevalence_period.txt", sep = "\t")

p1 = stays %>%
  filter(PATPER != periods[2]) %>%
  select(SUBJID, PATPER, REAENT, PATSTART, PATEND, REASOR, INFCOV) %>%
  mutate(DATE = map2(REAENT, REASOR, seq, by = "1 day")) %>%
  unnest(cols = DATE) %>%
  filter(DATE >= PATSTART, DATE <= PATEND) %>%
  group_by(PATPER, DATE) %>%
  summarise(p_cov = sum(INFCOV %in% "OUI") / n(), .groups = "drop") %>%
  mutate(PATPER = factor(PATPER, periods[-2])) %>%
  ggplot(., aes(x = DATE, y = p_cov)) +
  geom_line() +
  facet_grid(cols = vars(PATPER), scales = "free_x") +
  theme_bw() +
  ylim(c(0,1)) +
  labs(x = "", y = "Prevalence of COVID-19 patients") 

# Save intubation prevalence per period
n_ind_per_period = stays %>%
  filter(PATPER != periods[2]) %>%
  select(SUBJID, PATPER, REAENT, PATSTART, PATEND, REASOR, INFCOV) %>%
  mutate(DATE = map2(REAENT, REASOR, seq, by = "1 day")) %>%
  unnest(cols = DATE) %>%
  filter(DATE >= PATSTART, DATE <= PATEND) %>%
  group_by(PATPER, DATE) %>%
  summarise(n_ind = n(), .groups = "drop")
  
stays %>%
  filter(PATPER != periods[2]) %>%
  select(SUBJID, PATPER, PATSTART, PATEND, contains("REAINT")) %>%
  pivot_longer(
    cols = matches("REAINT.+"),
    cols_vary = "slowest",
    names_to = c("set", ".value"),
    names_pattern = "REAINT(.)(.)"
  ) %>%
  rename(REAINTF = `F`, REAINTD = D) %>%
  filter(!is.na(REAINTF), !is.na(REAINTD), REAINT %in% "OUI") %>%
  mutate(DATE = map2(REAINTD, REAINTF, seq, by = "1 day")) %>%
  unnest(cols = DATE) %>%
  filter(DATE >= PATSTART, DATE <= PATEND) %>%
  group_by(PATPER, DATE) %>%
  summarise(n_intub = sum(REAINT %in% "OUI"), .groups = "drop") %>%
  left_join(., n_ind_per_period, by = c("PATPER", "DATE")) %>%
  mutate(p_intub = n_intub / n_ind) %>%
  write.table(., "data/AMR_COVID-Data-230803/intubation_prevalence_period.txt", sep = "\t")

p2 = stays %>%
  filter(PATPER != periods[2]) %>%
  select(SUBJID, PATPER, PATSTART, PATEND, contains("REAINT")) %>%
  pivot_longer(
    cols = matches("REAINT.+"),
    cols_vary = "slowest",
    names_to = c("set", ".value"),
    names_pattern = "REAINT(.)(.)"
  ) %>%
  rename(REAINTF = `F`, REAINTD = D) %>%
  filter(!is.na(REAINTF), !is.na(REAINTD), REAINT %in% "OUI") %>%
  mutate(DATE = map2(REAINTD, REAINTF, seq, by = "1 day")) %>%
  unnest(cols = DATE) %>%
  filter(DATE >= PATSTART, DATE <= PATEND) %>%
  group_by(PATPER, DATE) %>%
  summarise(n_intub = sum(REAINT %in% "OUI"), .groups = "drop") %>%
  left_join(., n_ind_per_period, by = c("PATPER", "DATE")) %>%
  mutate(p_intub = n_intub / n_ind) %>%
  mutate(PATPER = factor(PATPER, periods[-2])) %>%
  ggplot(., aes(x = DATE, y = p_intub)) +
  geom_line() +
  facet_grid(cols = vars(PATPER), scales = "free_x") +
  theme_bw() +
  ylim(c(0,1)) +
  labs(x = "", y = "Prevalence of intubated patients") 


p3 = n_ind_per_period %>%
  mutate(PATPER = factor(PATPER, periods[-2])) %>%
  ggplot(., aes(x = DATE, y = n_ind)) +
  facet_grid(cols = vars(PATPER), scales = "free_x") +
  geom_line() +
  theme_bw() +
  labs(x = "", y = "Number of individuals hospitalized in ICU") +
  ylim(c(0,65))
  
to_save = ggarrange(p1, p2, p3, ncol = 1, nrow = 3)
ggsave("figures/AMR_COVID-Data-230803/intub_covid_prev.png", to_save,
       height = 10, width = 10)

# Plot for poster
intub_prev_df =  stays %>%
  filter(PATPER != periods[2]) %>%
  select(SUBJID, PATPER, PATSTART, PATEND, contains("REAINT")) %>%
  pivot_longer(
    cols = matches("REAINT.+"),
    cols_vary = "slowest",
    names_to = c("set", ".value"),
    names_pattern = "REAINT(.)(.)"
  ) %>%
  rename(REAINTF = `F`, REAINTD = D) %>%
  filter(!is.na(REAINTF), !is.na(REAINTD), REAINT %in% "OUI") %>%
  mutate(DATE = map2(REAINTD, REAINTF, seq, by = "1 day")) %>%
  unnest(cols = DATE) %>%
  filter(DATE >= PATSTART, DATE <= PATEND) %>%
  group_by(PATPER, DATE) %>%
  summarise(n_intub = sum(REAINT %in% "OUI"), .groups = "drop") %>%
  left_join(., n_ind_per_period, by = c("PATPER", "DATE")) %>%
  mutate(p_intub = n_intub / n_ind)

p4 = stays %>%
  filter(PATPER != periods[2]) %>%
  select(SUBJID, PATPER, REAENT, PATSTART, PATEND, REASOR, INFCOV) %>%
  mutate(DATE = map2(REAENT, REASOR, seq, by = "1 day")) %>%
  unnest(cols = DATE) %>%
  filter(DATE >= PATSTART, DATE <= PATEND) %>%
  group_by(PATPER, DATE) %>%
  summarise(p_cov = sum(INFCOV %in% "OUI") / n(), .groups = "drop") %>%
  left_join(., intub_prev_df, by = c("PATPER", "DATE")) %>%
  pivot_longer(c(p_cov, p_intub), names_to = "Prevalence", values_to = "value") %>%
  mutate(PATPER = case_when(
    PATPER == periods[1] ~ "Period 1",
    PATPER == periods[3] ~ "Period 2",
    PATPER == periods[4] ~ "Period 3"
  ),
    Prevalence = ifelse(Prevalence == "p_cov", "COVID-19", "Intubation")) %>%
  ggplot(., aes(x = DATE, y = value, col = Prevalence)) +
  geom_line(linewidth = 1) +
  facet_grid(cols = vars(PATPER), scales = "free_x") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 14),
    panel.spacing = unit(0.5, "cm")
  ) +
  ylim(c(0,1)) +
  scale_color_manual(values = c("#FF7F50", "#00808080")) +
  labs(x = "", y = "Daily prevalence", col = "") 

p5 = n_ind_per_period %>%
  mutate(PATPER = case_when(
    PATPER == periods[1] ~ "Period 1",
    PATPER == periods[3] ~ "Period 2",
    PATPER == periods[4] ~ "Period 3"
  )) %>%
  ggplot(., aes(x = DATE, y = n_ind)) +
  facet_grid(cols = vars(PATPER), scales = "free_x") +
  geom_line(linewidth = 1) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14),
    panel.spacing = unit(0.5, "cm")
  ) +
  labs(x = "", y = "Number of individuals\nhospitalized in ICU") +
  ylim(c(0,65))

tosave = ggarrange(p4, p5, ncol = 1, common.legend = T, legend = "right")
ggsave("figures/AMR_COVID-Data-230803/prevalences_poster.png", 
       tosave, height = 8, width = 13)

# Save data on stays for simulations
stays %>%
  select(SUBJID, PATPER, PATSTART, PATEND, REAENT, REASOR, REASEC1, REACHA1) %>%
  filter(PATPER != periods[2]) %>%
  mutate(REACHA1 = as.numeric(ifelse(REACHA1 %in% "K", 1, REACHA1))) %>%
  mutate(
    admission = if_else(REAENT < PATSTART, 0, as.numeric(REAENT-PATSTART)),
    discharge = if_else(REASOR > PATEND, as.numeric(PATEND-PATSTART), as.numeric(REASOR-PATSTART)),
    period = case_when(PATPER == periods[1] ~ 0, PATPER == periods[3] ~ 1, PATPER == periods[4] ~ 2),
    subsector = case_when(
      REASEC1 == "BLOC"                                      ~ 0,
      REASEC1 == "EXTENSION" & REACHA1 %in% 1:5              ~ 1,
      REASEC1 == "EXTENSION"& REACHA1 %in% 6:9               ~ 2,
      REASEC1 == "EXTENSION" & REACHA1 %in% 10:12            ~ 3,
      REASEC1 == "FLAMBOYANT" & REACHA1 %in% c(1,3,5,7,9)    ~ 4,
      REASEC1 == "FLAMBOYANT" & REACHA1 %in% c(2,4,6,8,10)   ~ 5,
      REASEC1 == "HIBISCUS" & REACHA1 %in% 1:3               ~ 6,
      REASEC1 == "HIBISCUS" & REACHA1 %in% 4:6               ~ 7,
      REASEC1 == "HIBISCUS" & REACHA1 %in% 7:8               ~ 8,
      REASEC1 == "MAHOGANY" & REACHA1 %in% 1:4               ~ 9,
      REASEC1 == "MAHOGANY" & REACHA1 %in% 5:8               ~ 10,
      REASEC1 == "SELF"                                      ~ 11,
      REASEC1 == "UHCD" & REACHA1 %in% 1:3                   ~ 12,
      REASEC1 == "UHCD" & REACHA1 %in% 4:8                   ~ 13
    )
  ) %>%
  arrange(period, admission, discharge) %>%
  group_by(period) %>%
  mutate(id = 0:(n()-1)) %>%
  write.table(., "data/AMR_COVID-Data-230803/occupancy_data.txt", sep = "\t")
  
# Plot on stays
to_save = 
ggsave("figures/AMR_COVID-Data-230803/n_patients_per_period.png", 
       to_save, height = 5, width = 10)

#################################################
## Intubation prevalence data -------------------
#################################################
# Plot intubation prevalence per period
intub_prev = list()
i= 1
for (p in periods[-2]){
  intub_prev[[i]] = stays %>%
    filter(PATPER == p) %>%
    select(SUBJID, PATPER, REAENT, PATSTART, PATEND, REASOR, REAINT) %>%
    mutate(DATE = map2(REAENT, REASOR, seq, by = "1 day"),
           PATPER = gsub(" \\(", "\n\\(", PATPER)) %>%
    unnest(cols = DATE) %>%
    filter(DATE >= PATSTART, DATE <= PATEND) %>%
    group_by(PATPER, DATE) %>%
    summarise(p_cov = sum(REAINT %in% "OUI") / n(), .groups = "drop") %>%
    ggplot(., aes(x = DATE, y = p_cov)) +
    geom_line() +
    facet_grid(cols = vars(PATPER)) +
    theme_minimal() +
    theme(legend.position = "none") +
    labs(x= "", y = "Intubation prevalence") +
    ylim(c(0,1))
  i=i+1
}
intub_prev_save = ggarrange(plotlist = intub_prev, ncol = 3, nrow = 1)
ggsave("figures/AMR_COVID-Data-230803/intub_prevalence_period.png", 
       intub_prev_save, height = 4, width = 10)

# Save covid-19 prevalence per period
stays %>%
  filter(PATPER != periods[2]) %>%
  select(SUBJID, PATPER, REAENT, PATSTART, PATEND, REASOR, REAINT) %>%
  mutate(DATE = map2(REAENT, REASOR, seq, by = "1 day")) %>%
  unnest(cols = DATE) %>%
  filter(DATE >= PATSTART, DATE <= PATEND) %>%
  group_by(PATPER, DATE) %>%
  summarise(p_cov = sum(REAINT %in% "OUI") / n(), .groups = "drop") %>%
  write.table(., "data/AMR_COVID-Data-230803/intub_prevalence_period.txt", sep = "\t")

#################################################
## Summary statistics on colonizations ----------
#################################################
# Number of colonization tests per patient
for (p in c(1,3,4)) {
  n_col_samples = stays %>% 
    filter(PATPER == periods[p]) %>%
    select(SUBJID, PATPER, PATAGE, PATSEX, REAENT, REASOR) %>%
    left_join(., microbio, by = "SUBJID") %>%
    filter(PVTNAT == "Pvt rectal") %>%
    pivot_longer(c(PVTKP, PVTECO, PVTECC, PVTEAG, PVTBA1), values_to = "result", names_to = "BACTERIA") %>%
    group_by(SUBJID, PATPER, PATAGE, PATSEX, BACTERIA) %>%
    summarise(n_samples = n(), .groups = "drop") %>%
    mutate(BACTERIA = recode(BACTERIA, !!!bacteria_species)) %>%
    ggplot(., aes(x = n_samples, fill = BACTERIA, col = BACTERIA)) +
    geom_bar() +
    theme_classic() +
    theme(legend.position = "none") +
    facet_grid(rows = vars(PATSEX), cols = vars(BACTERIA)) +
    labs(x = "Number of samples per patient during ICU stay",
         y = "Count", 
         fill = "", col = "")
  ggsave(paste0("figures/AMR_COVID-Data-230803/n_colonization_tests_", period_names_short[p],".png"), 
         n_col_samples, height = 4, width = 9)
}

# Time delay between two colonization tests
time_int = stays %>%
  select(SUBJID, PATPER) %>%
  left_join(., microbio, by = "SUBJID") %>%
  filter(PATPER != periods[2], PVTNAT == "Pvt rectal") %>%
  arrange(PATPER, SUBJID, PVTDAT) %>%
  group_by(PATPER, SUBJID) %>%
  summarise(
    test_interval = mean(as.numeric(difftime(PVTDAT, lag(PVTDAT), units = "day")), na.rm = T),
    .groups = "drop"
    ) %>%
  filter(!is.na(test_interval)) %>%
  mutate(PATPER = factor(PATPER, periods[c(1,3,4)], period_names[c(1,3,4)])) %>%
  ggplot(., aes(x = test_interval)) +
  geom_histogram() +
  geom_vline(aes(xintercept = median(test_interval)), col = "red") +
  facet_grid(cols = vars(PATPER)) +
  theme_classic() +
  labs(x = "Average time interval between two microbiological\ntests (in days) per individual", 
       y = "Count")
ggsave("figures/AMR_COVID-Data-230803/time_interval_tests.png", 
       time_int, height = 3, width = 7)

# Time delay between ICU admission and first colonization tests
time_first = stays %>%
  select(SUBJID, PATPER, REAENT) %>%
  left_join(., microbio, by = "SUBJID") %>%
  filter(PATPER != periods[2], PVTNAT == "Pvt rectal") %>%
  arrange(PATPER, SUBJID, PVTDAT) %>%
  group_by(PATPER, SUBJID) %>%
  summarise(
    first_delay = min(as.numeric(difftime(PVTDAT, REAENT, units = "day")), na.rm = T),
    .groups = "drop"
  ) %>%
  filter(!is.na(first_delay)) %>%
  mutate(PATPER = factor(PATPER, periods[c(1,3,4)], period_names[c(1,3,4)])) %>%
  ggplot(., aes(x = first_delay)) +
  geom_histogram() +
  geom_vline(aes(xintercept = median(first_delay)), col = "red") +
  facet_grid(cols = vars(PATPER)) +
  theme_classic() +
  labs(x = "Average time interval between ICU admission\nand first rectal sample", 
       y = "Count")
ggsave("figures/AMR_COVID-Data-230803/time_delay_first_sample.png", 
       time_first, height = 3, width = 7)


# Bar plot of the number of colonization samples by patient for 
# each species and period
for (p in c(1,3,4)) {
  n_col_samples = stays %>% 
    filter(PATPER == periods[p]) %>%
    select(SUBJID, PATPER, PATAGE, PATSEX, REAENT, REASOR) %>%
    left_join(., microbio, by = "SUBJID") %>%
    filter(PVTNAT == "Pvt rectal") %>%
    pivot_longer(c(PVTKP, PVTECO, PVTECC, PVTEAG, PVTBA1), values_to = "result", names_to = "BACTERIA") %>%
    group_by(SUBJID, PATPER, PATAGE, PATSEX, BACTERIA) %>%
    summarise(
      positive_sample = sum(result, na.rm = T),
      .groups = "drop"
    ) %>%
    mutate(BACTERIA = recode(BACTERIA, !!!bacteria_species)) %>%
    ggplot(., aes(x = positive_sample, fill = BACTERIA, col = BACTERIA)) +
    geom_bar() +
    theme_classic() +
    theme(legend.position = "none") +
    facet_grid(rows = vars(PATSEX), cols = vars(BACTERIA)) +
    labs(x = "Number of positive samples of ESBL-producing\nbacteria per patient during ICU stay",
         y = "Count", 
         fill = "", col = "")
  ggsave(paste0("figures/AMR_COVID-Data-230803/n_colonization_positives_", period_names_short[p],".png"), 
         n_col_samples, height = 4, width = 7)
}

# Table comparing time interval, average number of tests
# number of positive patients and number of positive tests 
# per patient between periods during the inclusion period
to_save = stays %>% 
  select(SUBJID, PATPER, REAENT, PATSTART, PATEND) %>%
  left_join(., microbio, by = "SUBJID") %>%
  filter(PVTNAT == "Pvt rectal", PVTDAT >= PATSTART, PVTDAT <= PATEND, PATPER != periods[2]) %>%
  mutate(PATPER = factor(PATPER, periods[c(1,3,4)])) %>%
  arrange(PATPER, SUBJID, PVTDAT) %>%
  group_by(SUBJID, PATPER) %>%
  summarise(
    n_samples = as.integer(n()),
    time_int = mean(as.numeric(difftime(PVTDAT, lag(PVTDAT), units = "day")), na.rm = T),
    first_sample_delay = min(difftime(PVTDAT, REAENT, units = "days")),
    n_pos_kp = any(PVTKP %in% 1),
    n_pos_eco = any(PVTECO %in% 1),
    n_pos_ecc = any(PVTECC %in% 1),
    n_pos_eag = any(PVTEAG %in% 1),
    n_pos_ba1 = any(PVTBA1 %in% 1),
    .groups = "drop"
  ) %>%
  select(-SUBJID) %>%
  gtsummary::tbl_summary(
    by = PATPER,
    label = list(
      n_samples ~ "Number of rectal samples",
      time_int ~ "Time interval between two tests (days)",
      first_sample_delay ~ "Delay from ICU admission to first rectal sample",
      n_pos_kp ~ "Positive patients for K pneumoniae",
      n_pos_eco ~ "Positive patients for E coli",
      n_pos_ecc ~ "Positive patients for E cloacae complex",
      n_pos_eag ~ "Positive patients for E aerogenes",
      n_pos_ba1 ~ "Positive patients for other bacteria"
    ),
    type = list(
      n_samples ~ "continuous"
    )
  ) %>%
  bold_labels() %>%
  modify_header(label ~ "") %>%
  add_overall(last = TRUE) %>%
  add_p()
gtsave(as_gt(to_save), path = "tables/AMR_COVID-Data-230803/", 
       filename = "microbio_description_period.png")

# Patients with multiple tests
multiple_pos_id = stays %>%
  select(SUBJID, PATPER, REAENT, REASOR, PATSTART, PATEND) %>%
  left_join(., microbio, by = "SUBJID") %>%
  filter(PVTNAT == "Pvt rectal", PATPER != periods[2]) %>%
  select(-c(PVTNAT, PVTLOC, PVTRES_IC, PVTRES, PVTBA1, PVTBA1P, PVTBA2, PVTBA2P, PVTBA3, PVTBA3P, X, n_species)) %>%
  pivot_longer(c(PVTKP, PVTECC, PVTECO, PVTEAG), values_to = "result", names_to = "BACTERIA") %>%
  group_by(SUBJID, PATPER, BACTERIA) %>%
  summarise(n_pos = sum(result %in% 1), .groups = "drop") %>%
  filter(n_pos > 1)

microbio %>%
  filter(SUBJID %in% multiple_pos_id$SUBJID)

# Compare incidence rates while including patients with no 
# rectal sample
n_hosp_days = stays %>%
  filter(PATPER != periods[2]) %>%
  mutate(n_jh = difftime(
    as.Date(ifelse(REASOR > PATEND, as.character(PATEND), as.character(REASOR))), 
    as.Date(ifelse(REAENT > PATSTART, as.character(REAENT), as.character(PATSTART))), 
    units = "days") + 1
    ) %>%
  group_by(PATPER) %>%
  summarise(n_jh = as.numeric(sum(n_jh)), .groups = "drop")

stays %>% 
  select(SUBJID, PATPER, PATSTART, PATEND) %>%
  filter(PATPER != periods[2]) %>%
  left_join(., microbio, by = "SUBJID") %>%
  mutate(PATPER = factor(PATPER, periods[c(1,3,4)])) %>%
  filter(PVTNAT == "Pvt rectal", PVTDAT <= PATEND) %>%
  arrange() %>%
  group_by(SUBJID, PATPER, PATSTART) %>%
  summarise(
    PVTKP = case_when(
      any(PVTKP %in% 1 & PVTDAT >= PATSTART) ~ 1,
      any(PVTKP %in% 1 & PVTDAT < PATSTART & !(lead(PVTKP) %in% 0 & lead(PVTDAT) < PATSTART)) ~ 1,
      .default = 0
    ),
    PVTECO = case_when(
      any(PVTECO %in% 1 & PVTDAT >= PATSTART) ~ 1,
      any(PVTECO %in% 1 & PVTDAT < PATSTART & !(lead(PVTECO) %in% 0 & lead(PVTDAT) < PATSTART)) ~ 1,
      .default = 0
    ), 
    PVTECC = case_when(
      any(PVTECC %in% 1 & PVTDAT >= PATSTART) ~ 1,
      any(PVTECC %in% 1 & PVTDAT < PATSTART & !(lead(PVTECC) %in% 0 & lead(PVTDAT) < PATSTART)) ~ 1,
      .default = 0
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(PVTKP, PVTECO, PVTECC), names_to = "BACTERIA", values_to = "pos") %>%
  group_by(PATPER, BACTERIA) %>%
  summarise(n_pos_patient = sum(pos), .groups = "drop") %>%
  left_join(., n_hosp_days, by = "PATPER") %>%
  mutate(
    ti = round(n_pos_patient / n_jh * 1000, 2),
    ti_upper = round(1000/n_jh * (n_pos_patient + 1.96 * sqrt(n_pos_patient)),2),
    ti_lower = round(1000/n_jh * (n_pos_patient - 1.96 * sqrt(n_pos_patient)),2),
    ) %>%
  pivot_wider(id_cols = -c(n_pos_patient), names_from = BACTERIA, values_from = ti) %>%
  write.table("tables/AMR_COVID-Data-230803/incidence_inclusions.tsv", 
              sep = "\t", row.names = F)


getAcquisitions = function(df) {
  result_out = c()
  
  if(nrow(df) == 0) {
    result_out = 0
  } else {
    
    admin = unique(df$REAENT)
    period_start = unique(df$PATSTART)
    period_end = unique(df$PATEND)
    
    df_temp = df %>%
      arrange(PVTDAT) %>%
      mutate(lag_result = lag(result))
    
    turning_pos = which(df_temp$result - df_temp$lag_result > 0)
    
    # Positive test 2 days after admission
    if(df_temp$result[1] %in% 1 & df_temp$PVTDAT[1] >= admin+2 & df_temp$PVTDAT[1] >= period_start) {
      result_out = 1
      
    # Negative test followed by positive test during inclusion period
    } else if (any( df_temp$PVTDAT[turning_pos] >= period_start )) {
      result_out = 1
      
    # All the other types of episodes
    } else {
      result_out = 0
    }
    
  }
  
  return(result_out)
}

stays %>% 
  select(SUBJID, REAENT, PATPER, PATSTART, PATEND) %>%
  filter(PATPER != periods[2]) %>%
  left_join(., microbio %>% filter(PVTNAT == "Pvt rectal") %>% select(SUBJID, PVTDAT, PVTKP, PVTECC, PVTECO), 
            by = "SUBJID") %>%
  mutate(PATPER = factor(PATPER, periods[c(1,3,4)])) %>%
  filter(PVTDAT <= PATEND) %>%
  pivot_longer(cols = c(PVTKP, PVTECO, PVTECC), names_to = "BACTERIA", values_to = "result") %>%
  group_by(SUBJID, PATPER, BACTERIA) %>%
  nest() %>%
  mutate(acquisition = map(data, function(.data) getAcquisitions(.data))) %>%
  select(-data) %>%
  unnest(cols = acquisition) %>%
  ungroup() %>%
  group_by(PATPER, BACTERIA) %>%
  summarise(n_pos_patient = sum(acquisition), .groups = "drop") %>%
  left_join(., n_hosp_days, by = "PATPER") %>%
  mutate(
    ti = round(n_pos_patient / n_jh * 1000, 2),
    ti_upper = round(1000/n_jh * (n_pos_patient + 1.96 * sqrt(n_pos_patient)),2),
    ti_lower = round(1000/n_jh * (n_pos_patient - 1.96 * sqrt(n_pos_patient)),2),
  )

# Compare incidence rates while excluding patients with no
# rectal sample
n_hosp_days = stays %>%
  filter(PATPER != periods[2], SUBJID %in% microbio$SUBJID[microbio$PVTNAT == "Pvt rectal"]) %>%
  mutate(n_jh = difftime(
    as.Date(ifelse(REASOR > PATEND, as.character(PATEND), as.character(REASOR))), 
    as.Date(ifelse(REAENT > PATSTART, as.character(REAENT), as.character(PATSTART))), 
    units = "days") + 1
  ) %>%
  group_by(PATPER) %>%
  summarise(n_jh = as.numeric(sum(n_jh)), .groups = "drop")

stays %>% 
  select(SUBJID, PATPER, PATSTART, PATEND) %>%
  filter(PATPER != periods[2]) %>%
  left_join(., microbio, by = "SUBJID") %>%
  mutate(PATPER = factor(PATPER, periods[c(1,3,4)])) %>%
  filter(PVTNAT == "Pvt rectal", PVTDAT <= PATEND) %>%
  group_by(SUBJID, PATPER, PATSTART) %>%
  summarise(
    PVTKP = case_when(
      any(PVTKP %in% 1 & PVTDAT >= PATSTART) ~ 1,
      any(PVTKP %in% 1 & PVTDAT < PATSTART & !(lead(PVTKP) %in% 0 & lead(PVTDAT) < PATSTART)) ~ 1,
      .default = 0
    ),
    PVTECO = case_when(
      any(PVTECO %in% 1 & PVTDAT >= PATSTART) ~ 1,
      any(PVTECO %in% 1 & PVTDAT < PATSTART & !(lead(PVTECO) %in% 0 & lead(PVTDAT) < PATSTART)) ~ 1,
      .default = 0
    ), 
    PVTECC = case_when(
      any(PVTECC %in% 1 & PVTDAT >= PATSTART) ~ 1,
      any(PVTECC %in% 1 & PVTDAT < PATSTART & !(lead(PVTECC) %in% 0 & lead(PVTDAT) < PATSTART)) ~ 1,
      .default = 0
    ),
    .groups = "drop"
  ) %>%
  pivot_longer(cols = c(PVTKP, PVTECO, PVTECC), names_to = "BACTERIA", values_to = "pos") %>%
  group_by(PATPER, BACTERIA) %>%
  summarise(n_pos_patient = sum(pos), .groups = "drop") %>%
  left_join(., n_hosp_days, by = "PATPER") %>%
  mutate(ti = round(n_pos_patient / n_jh * 1000,2)) %>%
  pivot_wider(id_cols = -c(n_pos_patient), names_from = BACTERIA, values_from = ti) %>%
  write.table("tables/AMR_COVID-Data-230803/incidence_exclusions.tsv", 
              sep = "\t", row.names = F)

# Comparison of individuals without rectal sample and individuals with rectal 
# sample
for (p in c(1,3,4)) {
  to_save = stays %>%
    filter(PATPER == periods[p]) %>%
    mutate(cat = ifelse(SUBJID %in% microbio$SUBJID[microbio$PVTNAT == "Pvt rectal"], "Sampled patient", "Unsampled patient")) %>%
    select(cat, PATAGE, PATSEX, stay_days, INFCOV, ATB, ATB48H, REABCA, REAINT, PATDEC) %>%
    gtsummary::tbl_summary(
      by = cat,
      label = list(
        PATAGE ~ "Age",
        PATSEX ~ "Gender",
        stay_days ~ "Stay length (days)",
        INFCOV ~ "COVID-19 infection",
        ATB ~ "Antibiotics at admission",
        ATB48H ~ "Antibiotics > 48h",
        REABCA ~ "Known infection pre-admission",
        REAINT ~ "Intubation",
        PATDEC ~ "Death"
      )
    ) %>%
    modify_header(label ~ period_names_bold[p]) %>%
    bold_labels() %>%
    add_p()
  gtsave(as_gt(to_save), path = "tables/AMR_COVID-Data-230803/", 
         filename = paste0("microbio_sampling_comparison_", period_names_short[p] ,".png"))
}

# Patients without rectal sample but tested for infection at other sampling sites
stays %>%
  filter(PATPER != periods[2]) %>%
  filter(SUBJID %in% c(microbio %>% group_by(SUBJID) %>% summarise(include = as.numeric(!"Pvt rectal" %in% PVTNAT), .groups = "drop") %>% filter(include == 1) %>% .$SUBJID)) %>%
  select(SUBJID, PATPER, PATAGE, PATSEX, REABCA, PATDEC) %>%
  left_join(., microbio, by = c("SUBJID")) %>%
  group_by(SUBJID, PATPER, PATAGE, PATSEX, REABCA, PATDEC) %>%
  summarise(PVTKP = any(PVTKP %in% 1), PVTECC = any(PVTECC %in% 1), PVTECO = any(PVTECO %in% 1),
            PVTEAG = any(PVTEAG %in% 1), col = any(PVTRES_IC %in% "Colonisation"), 
            inf = any(PVTRES_IC %in% "Infection"), .groups = "drop") %>%
  select(-SUBJID) %>%
  gtsummary::tbl_summary(
    by = PATPER,
    label = list(
      PATAGE ~ "Age",
      PATSEX ~ "Gender",
      REABCA ~ "Known infection pre-admission",
      PATDEC ~ "Death", 
      PVTKP ~ "K pneumoniae",
      PVTECC ~ "E cloacae complex",
      PVTECO ~ "E coli",
      PVTEAG ~ "E aerogenes",
      col ~ "Colonizations",
      inf ~ "Infections"
    )
  ) %>%
  bold_labels() %>%
  add_p()


# Patients without any sample
to_save = stays %>%
  filter(PATPER != periods[2]) %>%
  mutate(
    sample_result = case_when(
      SUBJID %in% c(microbio %>% filter(PVTNAT %in% "Pvt rectal", PVTRES %in% "Positif") %>% .$SUBJID) ~ "Positive rectal sample",
      SUBJID %in% c(microbio %>% filter(PVTNAT %in% "Pvt rectal", PVTRES %in% "Négatif") %>% .$SUBJID) ~ "Negative rectal sample",
      SUBJID %in% c(microbio %>% filter(!PVTNAT %in% "Pvt rectal", PVTRES_IC %in% "Colonisation" & !PVTRES_IC %in% c(" ", "Infection")) %>% .$SUBJID) ~ "Positive colonization sample",
      SUBJID %in% c(microbio %>% filter(!PVTNAT %in% "Pvt rectal", PVTRES_IC %in% "Infection" & !PVTRES_IC %in% c(" ", "Colonisation")) %>% .$SUBJID) ~ "Positive infection sample",
      SUBJID %in% c(microbio %>% filter(!PVTNAT %in% "Pvt rectal", PVTRES_IC %in% " " & !PVTRES_IC %in% c("Infection", "Colonisation")) %>% .$SUBJID) ~ "Negative non-rectal sample type",
      !SUBJID %in% microbio$SUBJID ~ "Unsampled",
      .default = NA
      ),
    sample_type = case_when(
      SUBJID %in% c(microbio %>% filter(PVTNAT %in% "Pvt rectal") %>% .$SUBJID) ~ "Rectal sample",
      SUBJID %in% c(microbio %>% filter(!PVTNAT %in% "Pvt rectal") %>% .$SUBJID) ~ "Other sample",
      !SUBJID %in% microbio$SUBJID ~ "Unsampled",
      .default = NA
    )
    ) %>%
  select(PATPER, PATSEX, PATAGE, REABCA, sample_type, sample_result) %>%
  gtsummary::tbl_summary(
    by = PATPER,
    label = list(
      PATAGE ~ "Age",
      PATSEX ~ "Gender",
      REABCA ~ "Known infection pre-admission",
      sample_type ~ "Sample type",
      sample_result ~ "Sample result"
    )
  ) %>%
  bold_labels() %>%
  add_p()
gtsave(as_gt(to_save), path = "tables/AMR_COVID-Data-230803/", 
       filename = "infection_colonisation_samples.png")


test = microbio %>%
  mutate(PVTNAT2 = ifelse(PVTNAT %in% "Pvt rectal", "Pvt rectal", "Autre pvt")) %>%
  group_by(SUBJID, PVTNAT2, PVTRES_IC) %>%
  summarise(PVTKP = any(PVTKP %in% 1), PVTECC = any(PVTECC %in% 1), PVTECO = any(PVTECO %in% 1),
            PVTEAG = any(PVTEAG %in% 1), .groups = "drop") %>%
  mutate(n_species = as.numeric(PVTKP) + as.numeric(PVTECC) + as.numeric(PVTECO) + as.numeric(PVTEAG))

table(test$PVTRES_IC, test$n_species)

#################################################
## Summary statistics on colonizations and 
## infections -----------------------------------
#################################################
# Get classification 
classification = stays %>%
  filter(PATPER != periods[2]) %>%
  mutate(admission_rea = as.numeric(difftime(REAENT, HOSENT, units = "days"))) %>%
  select(SUBJID, PATPER, PATSTART, PATEND, REAENT, PATSEX, PATAGE, stay_days, admission_rea, 
         ATB48H, ATB, REAINT, INFCOV, REABCA) %>%
  mutate(
    ATB48H = case_when(ATB48H %in% "OUI" ~ 1, ATB48H %in% "NON" ~ 0, .default = NA),
    ATB = case_when(ATB %in% "OUI" ~ 1, ATB %in% "NON" ~ 0, .default = NA),
    REAINT = case_when(REAINT %in% "OUI" ~ 1, REAINT %in% "NON" ~ 0, .default = NA),
    INFCOV = case_when(INFCOV %in% "OUI" ~ 1, INFCOV %in% "NON" ~ 0, .default = NA),
    REABCA = case_when(REABCA %in% "OUI" ~1, REABCA %in% "NON" ~ 0, .default = NA),
    PATPER = factor(PATPER)
    ) %>%
  left_join(., microbio %>% select(SUBJID, PVTDAT, PVTNAT, PVTRES_IC, PVTKP, PVTECC, PVTECO, PVTEAG), 
            by = "SUBJID") %>%
  pivot_longer(cols = c(PVTKP, PVTECC, PVTECO, PVTEAG), names_to = "BACTERIE", values_to = "Result") %>%
  group_by(SUBJID, PATPER, PATSEX, PATAGE, stay_days, admission_rea, ATB48H, ATB, REAINT, INFCOV, REABCA, BACTERIE) %>%
  nest() %>%
  mutate(results = map(data, function(.data) getResults(.data))) %>%
  select(-data) %>%
  unnest(cols = results) %>%
  ungroup()

# Summary tables of acquisitions per period
to_save = classification %>%
  group_by(PATPER, SUBJID) %>%
  summarise(
    RESULT = case_when(
      any(RESULT %in% "Positive") ~ "Positive",
      any(RESULT %in% "Acquisition") & !any(RESULT %in% "Positive") ~ "Acquisition",
      all(RESULT %in% "Negative") ~ "Negative",
      .default = NA    
    ), .groups = "drop"
  ) %>%
  mutate(
    RESULT2 = ifelse(RESULT %in% "Positive", NA, RESULT),
    PATPER = factor(PATPER, periods[-2])
  ) %>%
  select(PATPER, RESULT, RESULT2) %>%
  gtsummary::tbl_summary(
    by = PATPER,
    label = list(RESULT ~ "All results", RESULT2 ~ "Results for susceptibles")
  ) %>%
  bold_labels() %>%
  modify_header(label ~ "") %>%
  add_p()
gtsave(as_gt(to_save), path = "tables/AMR_COVID-Data-230803/", 
       filename = "acquisitions_per_period.png")


# Intubation associated with acquisition ?
for (i in c(1,3,4)) {
  to_save = classification %>%
    filter(PATPER == periods[i]) %>%
    group_by(SUBJID, REAINT) %>%
    summarise(
      RESULT = case_when(
        any(RESULT %in% "Positive") ~ "Positive",
        any(RESULT %in% "Acquisition") & !any(RESULT %in% "Positive") ~ "Acquisition",
        all(RESULT %in% "Negative") ~ "Negative",
        .default = NA    
      ), .groups = "drop"
    ) %>%
    mutate(
      RESULT2 = ifelse(RESULT %in% "Positive", NA, RESULT),
      REAINT = ifelse(REAINT == 0, "Non intubé", "Intubé")
    ) %>%
    select(REAINT, RESULT, RESULT2) %>%
    gtsummary::tbl_summary(
      by = REAINT,
      label = list(RESULT ~ "All results", RESULT2 ~ "Results for susceptibles")
    ) %>%
    bold_labels() %>%
    modify_header(label ~ period_names_bold[i]) %>%
    add_p()
  gtsave(as_gt(to_save), path = "tables/AMR_COVID-Data-230803/", 
         filename = paste0("intubation_acquisition_", period_names_short[i], ".png"))  
}



# Summary tables for non K pneumoniae acquisitions
to_save = classification %>%
  select(BACTERIE, RESULT) %>%
  mutate(BACTERIE = case_when(
    BACTERIE == "PVTKP" ~ "Klebsiella pneumoniae",
    BACTERIE == "PVTECC" ~ "Enterobacter cloacae complex",
    BACTERIE == "PVTECO" ~ "Escherichia coli",
    BACTERIE == "PVTEAG" ~ "Enterobacter aerogenes"
  )) %>%
  gtsummary::tbl_summary(
    by = BACTERIE,
    label = list(RESULT ~ "Result")
    ) %>%
  bold_labels() 
gtsave(as_gt(to_save), path = "tables/AMR_COVID-Data-230803/", 
       filename = "result_per_species.png")

# Summary table per period for K pneumoniae
to_save = classification %>%
  filter(BACTERIE == "PVTKP") %>%
  select(-c(SUBJID, BACTERIE, PVTDAT)) %>%
  gtsummary::tbl_summary(
    by = PATPER,
    label = list(
      PATAGE ~ "Age",
      PATSEX ~ "Gender",
      stay_days ~ "Stay length",
      admission_rea ~ "Delay from hospitalisation to ICU",
      ATB48H ~ "Antibiotic > 48H",
      ATB ~ "Antibiotic < 48H",
      REAINT ~ "Intubation",
      INFCOV ~ "COVID-19",
      PVTRECTAL ~ "Rectal sample",
      REABCA ~ "Infection prior ICU admission",
      RESULT ~ "Result"
    )
  ) %>%
  bold_labels() %>%
  modify_header(label ~ "**K pneumoniae**") %>%
  add_p()
gtsave(as_gt(to_save), path = "tables/AMR_COVID-Data-230803/", 
       filename = "kp_per_period.png")

# Univariate analysis
to_save = classification %>%
  filter(BACTERIE == "PVTKP") %>%
  select(-c(SUBJID, BACTERIE, PVTDAT)) %>%
  mutate(RESULT = ifelse(is.na(RESULT), "Unsampled", RESULT)) %>%
  gtsummary::tbl_summary(
    by = RESULT,
    label = list(
      PATAGE ~ "Age",
      PATSEX ~ "Gender",
      stay_days ~ "Stay length",
      admission_rea ~ "Delay from hospitalisation to ICU",
      ATB48H ~ "Antibiotic > 48H",
      ATB ~ "Antibiotic < 48H",
      REABCA ~ "Infection prior ICU admission",
      REAINT ~ "Intubation",
      INFCOV ~ "COVID-19",
      PVTRECTAL ~ "Rectal sample",
      PATPER ~ "Inclusion period"
    )
  ) %>%
  bold_labels() %>%
  modify_header(label ~ "**K pneumoniae**") %>%
  add_p()
gtsave(as_gt(to_save), path = "tables/AMR_COVID-Data-230803/", 
       filename = "kp_univariate.png")

# Univariate analysis including only acquisitions and negative patients
to_save = classification %>%
  filter(BACTERIE == "PVTKP", RESULT %in% c("Negative", "Acquisition")) %>%
  select(-c(SUBJID, BACTERIE, PVTDAT)) %>%
  mutate(RESULT = ifelse(is.na(RESULT), "Unsampled", RESULT)) %>%
  gtsummary::tbl_summary(
    by = RESULT,
    label = list(
      PATAGE ~ "Age",
      PATSEX ~ "Gender",
      stay_days ~ "Stay length",
      admission_rea ~ "Delay from hospitalisation to ICU",
      ATB48H ~ "Antibiotic > 48H",
      ATB ~ "Antibiotic < 48H",
      REABCA ~ "Infection prior ICU admission",
      REAINT ~ "Intubation",
      INFCOV ~ "COVID-19",
      PVTRECTAL ~ "Rectal sample",
      PATPER ~ "Inclusion period"
    )
  ) %>%
  bold_labels() %>%
  modify_header(label ~ "**K pneumoniae**") %>%
  add_p()
gtsave(as_gt(to_save), path = "tables/AMR_COVID-Data-230803/", 
       filename = "kp_univariate_neg_acq.png")

# Regression models
data_model = classification %>% 
  filter(BACTERIE == "PVTKP") %>%
  mutate(RESULT = ifelse(RESULT %in% "Acquisition", 1, 0),
         PATPER = factor(PATPER, periods[-2]))

nullmod <- glm(RESULT~1, family = binomial(link = "logit"), data = data_model)

full_model = glm(RESULT ~ PATPER + stay_days + admission_rea + ATB + ATB48H + INFCOV, 
                 family = binomial(link = "logit"), 
                 data = data_model
)

step.model <- stepAIC(full_model, direction = "both", trace = FALSE)

# Diagnostics 
summary(step.model)
round(exp(step.model$coefficients), 3)
round(exp(confint(step.model)), 3)
1-logLik(step.model)/logLik(nullmod)

# Final model
final_model = glm(RESULT ~ PATPER + admission_rea + ATB48H + INFCOV, 
                  family = binomial(link = "logit"), 
                  data = data_model
)
summary(final_model)
round(exp(final_model$coefficients), 3)
round(exp(confint(final_model)), 3)
anova(final_model)
qqnorm(residuals(final_model))
hl <- hoslem.test(final_model$y, fitted(final_model), g=10)
hl
cbind(hl$observed,hl$expected)
1-logLik(final_model)/logLik(nullmod)

# Where acquisitions occured
classification %>%
  filter(BACTERIE == "PVTKP", RESULT %in% "Acquisition") %>%
  left_join(., stays %>% select(SUBJID, REASEC1, REASEC2, REASEC3, REASEC4), 
            by = "SUBJID") %>%
  select(contains("REASEC"), PATPER) %>%
  mutate(
    REASEC2 = ifelse(REASEC2 == " ", NA, REASEC2),
    REASEC3 = ifelse(REASEC3 == " ", NA, REASEC3),
    REASEC4 = ifelse(REASEC4 == " ", NA, REASEC4)
  ) %>%
  gtsummary::tbl_summary(
    by = PATPER
  )

all_plots = list()
i = 1
secteurs = c("BLOC" = "#DBA92A", 
             "EXTENSION" = "#3D5A5C", "UHCD" = "darkgreen",
             "FLAMBOYANT" = "#CF5423", "HIBISCUS" = "red", "MAHOGANY" = "darkred", "SELF" = "pink")

for (p in c(1,3,4)) {
  all_plots[[i]] = classification %>%
    filter(BACTERIE == "PVTKP", RESULT %in% "Acquisition", PATPER == periods[p]) %>%
    left_join(., stays %>% select(SUBJID, REASEC1, REASEC2, REASEC3, REASEC4, 
                                  REASOR1, REASOR2, REASOR3, REASOR4, 
                                  PATSTART, PATEND, REAENT, REASOR), 
              by = "SUBJID") %>%
    mutate(
      REASEC1 = ifelse(REASEC1 %in% " ", NA, REASEC1),
      REASEC2 = ifelse(REASEC2 %in% " ", NA, REASEC2),
      REASEC3 = ifelse(REASEC3 %in% " ", NA, REASEC3),
      REASEC4 = ifelse(REASEC4 %in% " ", NA, REASEC4)
      ) %>%
    mutate(REASEC1 = factor(REASEC1, levels = unique(c(levels(REASEC1), names(secteurs))))) %>%
    tidyr::complete(REASEC1) %>%
    ggplot(.) +
    geom_vline(aes(xintercept = PATSTART)) +
    geom_vline(aes(xintercept = PATEND)) +
    geom_segment(aes(x = REAENT, xend = REASOR1, y = SUBJID, yend = SUBJID, col = REASEC1)) +
    geom_segment(aes(x = REASOR1, xend = REASOR2, y = SUBJID, yend = SUBJID, col = REASEC2)) +
    geom_segment(aes(x = REASOR2, xend = REASOR3, y = SUBJID, yend = SUBJID, col = REASEC3)) +
    geom_segment(aes(x = REASOR3, xend = REASOR4, y = SUBJID, yend = SUBJID, col = REASEC4)) +
    geom_point(aes(x = PVTDAT, y = SUBJID)) +
    scale_color_manual(values = secteurs) +
    theme_light() +
    labs(title = periods[p], x = "", col = "")
  i=i+1
}

acq_sector = ggarrange(plotlist = all_plots, ncol = 2, nrow = 2, common.legend = T)
ggsave("figures/AMR_COVID-Data-230803/acq_secteur.png", acq_sector, height = 15, width = 15)

table(stays$PATPER, stays$REASEC1)

#################################################
## Plot episodes --------------------------------
#################################################
for (p in c(1,3,4)) {
  for (b in c("K pneumoniae", "E cloacae complex", "E coli")) {
    if (b == "K pneumoniae") bac = "PVTKP"
    if (b == "E cloacae complex") bac = "PVTECC"
    if (b == "E coli") bac = "PVTECO"
    
    to_save = stays %>% 
      filter(PATPER == periods[p]) %>%
      select(SUBJID, PATSEX, PATAGE, REAENT, REASOR, INFCOV) %>%
      left_join(., microbio %>% filter(PVTNAT == "Pvt rectal"), by = "SUBJID") %>%
      mutate(
        PVTKP = factor(PVTKP, c(1,0), c("positive", "negative")),
        PVTECO = factor(PVTECO, c(1,0), c("positive", "negative")),
        PVTECC = factor(PVTECC, c(1,0), c("positive", "negative")),
        INFCOV = factor(INFCOV, c("OUI", "NON"), c("positive", "negative"))
      ) %>%
      ggplot(.) +
      geom_rect(aes(xmin = periods_start[p], xmax = periods_end[p], ymin = -Inf, ymax = Inf), fill = "grey90") +
      geom_segment(aes(x = REAENT, xend = REASOR, y = SUBJID, yend = SUBJID, color = PATSEX)) +
      geom_point(aes(x = PVTDAT, y = SUBJID, color = !!sym(bac)), shape = 5) +
      geom_rug(aes(y = SUBJID, color = INFCOV), sides = "r", linewidth = 1) +
      theme_minimal() +
      theme(panel.grid.minor = element_blank(), 
            axis.title = element_text(size = 16), 
            axis.text.x = element_text(size = 14), 
            legend.text = element_text(size = 14)) +
      scale_color_manual(name = "",
                         values = c("orange", "cornflowerblue", "grey50", "red"),
                         labels = c("Feminin", "Masculin", "negative", "positive"),
                         na.value = "grey70") +
      labs(x = "", y = "Patients ID", title = paste(periods[p], "-", b))
    
    ggsave(paste0("figures/AMR_COVID-Data-230803/microbio_", tolower(gsub("PVT", "", bac)), "_", period_names_short[p], ".png"), 
           to_save, height = 12, width = 8)
  }
}

#################################################
## Logistic regression --------------------------
#################################################
# Regression model
glm_data = stays %>%
  mutate(STAY_LENGTH = difftime(
    as.Date(ifelse(REASOR > PATEND, as.character(PATEND), as.character(REASOR))), 
    as.Date(ifelse(REAENT > PATSTART, as.character(REAENT), as.character(PATSTART))), 
    units = "days")
    ) %>%
  left_join(., microbio %>% filter(PVTNAT == "Pvt rectal"), by = "SUBJID") %>%
  filter(PATPER != periods[2]) %>%
  select(SUBJID, PATPER, PATSEX, PATAGE, ATB, ATB48H, STAY_LENGTH, PVTDAT, PVTKP, PVTECC, PVTECO, PVTEAG, PATSTART, PATEND) %>%
  pivot_longer(cols = c(PVTKP, PVTECC, PVTECO, PVTEAG), names_to = "BACTERIE", values_to = "result") %>%
  group_by(SUBJID, PATPER, PATSEX, PATAGE, ATB, ATB48H, STAY_LENGTH, BACTERIE) %>%
  nest() %>%
  mutate(test_results = map(data, function(.data) summarizeTestResults(.data))) %>%
  select(-data) %>%
  unnest(cols = test_results) %>%
  ungroup() %>%
  mutate(
    BACTERIE = gsub("PVT", "", BACTERIE),
    test_results_final = case_when(
      test_results == 0 ~ 0,
      test_results == 1 ~ 0,
      test_results == 2 ~ 1,
      test_results == 3 & first_pos_before == 1 ~ 1,
      test_results == 3 & first_pos_before == 0 ~ 0,
      test_results == 4 & pos_inc == 1 ~ 1,
      test_results == 4 & pos_inc == 0 ~ 0,
      test_results == 5 ~ 1,
    ),
    PATPER = factor(PATPER, periods[c(1,3,4)])
    ) %>%
  select(-c(test_results,first_pos_before, pos_inc)) %>%
  pivot_wider(names_from = BACTERIE, values_from = test_results_final)


full_model = glm(KP ~ PATPER + PATSEX + PATAGE + ATB + ATB48H + STAY_LENGTH, 
                 family = binomial(link = "logit"), 
                 data = glm_data %>% select(-c(ECC, ECO, EAG)))
step.model <- stepAIC(full_model, direction = "both", trace = FALSE)
summary(step.model)
round(exp(step.model$coefficients), 3)
round(exp(confint(step.model)), 3)
anova(step.model)
qqnorm(residuals(step.model))

# Mc Fadden R2
nullmod <- glm(KP~1, family = binomial(link = "logit"), data = glm_data %>% select(-c(ECC, ECO, EAG)))
1-logLik(step.model)/logLik(nullmod)

# Individuals with complicated colonization sequence
complex_infection = microbio %>%
  filter(SUBJID %in% c("492-1", "498-1", "610-1"), PVTRES == "Positif", PVTRES_IC == "Infection") %>%
  mutate(PVTKP = ifelse(PVTKP == 1, "positive", "negative"))

microbio %>% 
  filter(SUBJID %in% c("492-1", "498-1", "610-1"), PVTNAT == "Pvt rectal") %>%
  left_join(., stays %>% select(SUBJID, PATSEX, PATPER, REAENT, REASOR, PATSTART, PATEND), by = "SUBJID") %>%
  mutate(PVTKP = ifelse(PVTKP == 1, "positive", "negative")) %>%
  ggplot(.) +
  geom_segment(aes(x = REAENT, xend = REASOR, y = SUBJID, yend = SUBJID, color = PATSEX)) +
  geom_vline(aes(xintercept = PATSTART)) +
  geom_vline(aes(xintercept = PATEND)) +
  geom_point(aes(x = PVTDAT, y = SUBJID, shape = "Colonization", color = PVTKP)) +
  geom_point(data = complex_infection, aes(x = PVTDAT, y = SUBJID, color = PVTKP, shape = PVTNAT)) +
  theme_minimal() +
  theme(panel.grid.minor = element_blank(), 
        axis.title = element_text(size = 16), 
        axis.text.x = element_text(size = 14), 
        legend.text = element_text(size = 14)) +
  scale_color_manual(name = "",
                     values = c("orange", "cornflowerblue", "grey50", "red"),
                     labels = c("Feminin", "Masculin", "negative", "positive"),
                     na.value = "grey70") +
  labs(x = "", y = "Patients ID") +
  scale_shape_manual(name = "", 
                     values = c(3, 2,1,4,5), 
                     labels = c("Colonization", "LBA", "PDP", "Sang cathéter", "Vaginal"))
ggsave("figures/AMR_COVID-Data-230803/patients_complicated_sequence.png", 
       height = 4, width = 6)


#################################################
## Information on episodes ----------------------
#################################################
# Individuals with known infection pre-ICU admission
table(stays$REABCA, stays$PATPER, useNA = "always", dnn = c("Pre-ICU infection", "Period")) 

stays %>%
  filter(REABCA %in% "OUI") %>%
  select(SUBJID, PATSEX, PATAGE, REAENT, REASOR, PATPER, REABCA, REABCAE, 
         REABCAD, stay_days) %>%
  left_join(., microbio %>% filter(PVTNAT == "Pvt rectal") %>% select(SUBJID, PVTDAT, PVTKP, PVTECC, PVTECO), by = "SUBJID") %>%
  mutate(
    REABCAE = case_when(
      REABCAE %in% c("KLEBSIELLA PNEUMONIAE", "KPE") ~ "K pneumoniae",
      REABCAE == "ESCHERICHIA COLI" ~ "E coli",
      REABCAE == "ENTEROBACTER CLOACAE" ~ "E cloacae",
      .default = NA
      ),
    PVTKP = ifelse(PVTKP == 1, "Positif", "Négatif")
    ) %>%
  ggplot(.) + 
  geom_segment(aes(y = SUBJID, yend = SUBJID, x = REAENT, xend = REASOR, col = PATSEX)) +  
  geom_point(aes(y = SUBJID, x = REABCAD, col = REABCAE), shape = 20) +
  geom_point(aes(y = SUBJID, x = PVTDAT, col = PVTKP), shape = 3) +
  scale_color_manual(labels = c("E cloacae", "E coli", "Féminin", 
                                "K pneumoniae","Masculin", "Négatif", "Positif"),
                     values = c("deeppink", "darkblue", "orange", "green2", "cornflowerblue", "grey70", "red")) +
  theme_classic() +
  labs(x = "Stay in ICU",  y = "Stay ID", col = "")
ggsave("figures/AMR_COVID-Data-230803/infection_pre_admission.png", height = 5, width = 8)

# Summary of the types of sequence
episode_type = stays %>%
  mutate(STAY_LENGTH = difftime(
    as.Date(ifelse(REASOR > PATEND, as.character(PATEND), as.character(REASOR))), 
    as.Date(ifelse(REAENT > PATSTART, as.character(REAENT), as.character(PATSTART))), 
    units = "days")
  ) %>%
  left_join(., microbio %>% filter(PVTNAT == "Pvt rectal"), by = "SUBJID") %>%
  filter(PATPER != periods[2]) %>%
  select(SUBJID, PATPER, PATSEX, PATAGE, ATB, ATB48H, STAY_LENGTH, PVTDAT, PVTKP, PVTECC, PVTECO, PVTEAG, PATSTART, PATEND) %>%
  pivot_longer(cols = c(PVTKP, PVTECC, PVTECO, PVTEAG), names_to = "BACTERIE", values_to = "result") %>%
  group_by(SUBJID, PATPER, PATSEX, PATAGE, ATB, ATB48H, STAY_LENGTH, BACTERIE) %>%
  nest() %>%
  mutate(test_results = map(data, function(.data) summarizeTestResults(.data))) %>%
  select(-data) %>%
  unnest(cols = test_results) %>%
  ungroup() %>%
  select(PATPER, BACTERIE, test_results) 

tbl_episodes = list()
index = 1
for (b in c("PVTKP", "PVTECC", "PVTECO")) {
  tbl_episodes[[index]] = episode_type %>%
    filter(BACTERIE == b) %>%
    select(-BACTERIE) %>%
    mutate(test_results = case_when(
      test_results == 0 ~ "Unsampled",
      test_results == 1 ~ "Negative",
      test_results == 2 ~ "Positive",
      test_results == 3 ~ "Acquisition",
      test_results == 4 ~ "Clearance",
      test_results == 5 ~ "Complex"
    ),
    PATPER = factor(PATPER, periods[c(1,3,4)])) %>%
    gtsummary::tbl_summary(
      by = PATPER,
      label = list(test_results ~ "Episode")) %>%
    bold_labels() %>%
    modify_header(label ~ "")
  index = index +1
}

tbl_merge(tbl_episodes, 
          c("**K pneumoniae**", "**E cloacae complex**", "**E coli**"))
# Manual saving 
# File name: "episode_type.png"
# Path: "tables/AMR_COVID-Data-230803"

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
gtsave(as_gt(to_save), path = "tables/", filename = "microbio_episode_category.png")

# Sectors and sub sectors where acquisitions occur


#################################################
## Contact network ------------------------------
#################################################
## Measures of disorganization
network_tab = stays %>%
  filter(PATPER != periods[2]) %>%
  rename(REAENT1 = REAENT) %>%
  select(PATPER, SUBJID, contains("REAENT"), contains("REASEC"), contains("REACHA"), matches("REASOR[0-9]")) %>%
  mutate(
    REAENT2 = as.Date(ifelse(!is.na(REASOR2), as.character(REASOR1), NA)), 
    REAENT3 = as.Date(ifelse(!is.na(REASOR3), as.character(REASOR2), NA)), 
    REAENT4 = as.Date(ifelse(!is.na(REASOR4), as.character(REASOR3), NA)), 
    REAENT5 = as.Date(ifelse(!is.na(REASOR5), as.character(REASOR4), NA)),
    REACHA1 = as.numeric(gsub("^0", "", ifelse(REACHA1 == "K", NA, REACHA1))),
    REACHA2 = as.numeric(gsub("^0", "", ifelse(REACHA2 == "K", NA, REACHA2))),
    REACHA3 = as.numeric(gsub("^0", "", ifelse(REACHA3 == "K", NA, REACHA3))),
    REACHA4 = as.numeric(gsub("^0", "", ifelse(REACHA4 == "K", NA, REACHA4))),
    REACHA5 = as.numeric(gsub("^0", "", ifelse(REACHA5 == "K", NA, REACHA5)))
  ) %>%
  pivot_longer(-c(PATPER, SUBJID), 
               names_pattern = "(.*)(.)", 
               names_to = c(".value", "STAYID")) %>%
  mutate(STAYID = as.numeric(STAYID)) %>%
  filter(!is.na(REASEC), !REASEC %in% " ") %>%
  arrange(SUBJID, REAENT)

# Diagram of stays per sector and room  
network_tab %>%
  ggplot(.) +
  geom_segment(aes(x = REAENT, xend = REASOR, y = SUBJID, yend = SUBJID, color = factor(REACHA))) +
  facet_grid(cols = vars(REASEC), rows = vars(PATPER)) +
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
  write.table(., "tables/room_sharing.txt", quote = F, row.names = F, sep = "\t")

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
  write.table(., "tables/transfers.txt", quote = F, row.names = F, sep = "\t")


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





