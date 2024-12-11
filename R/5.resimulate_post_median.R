##############################################
## SIMULATE FAKE DATASET TO TEST CPP CODE
##############################################

rm(list = ls())
library(tidyverse)
library(Rcpp)

source("R/helperFunctions.R")
source("R/dictionaries.R")
source("R/modelParams.R")

# Compile Rcpp code
filesToRemove = list.files(
  path = "cpp/simulation/",
  pattern = "^.*\\.o$",
  full.names = T
)
if (length(filesToRemove)) sapply(filesToRemove, file.remove)
sourceCpp("cpp/simulations/Colonization_simulation.cpp")

# Load or simulate prevalence data 
real_prevalence_data = T

# Load data
covid_data = read.table("data/AMR_COVID-Data-230803/covid_prevalence_period.txt", header = T)
intubation_data = read.table("data/AMR_COVID-Data-230803/intubation_prevalence_period.txt", header = T)
occupancy_data = read.table("data/AMR_COVID-Data-230803/occupancy_data.txt", header = T)

# Information on period length and subsectors per period
nDays = occupancy_data %>% arrange(period) %>% group_by(period) %>% summarise(m = max(discharge)+1) %>% .$m
nSS = occupancy_data %>% arrange(period) %>% group_by(period) %>% summarise(m = length(unique(subsector))) %>% .$m
listSS = lapply(0:2, function(x) unique(occupancy_data$subsector[occupancy_data$period == x]))

# Number of beds per period 
care_orga = occupancy_data %>%
  select(period, subsector, id, admission, discharge) %>%
  group_by(period, subsector) %>%
  mutate(day = map2(admission, discharge, seq, by = 1)) %>%
  unnest(cols = day) %>%
  ungroup() %>%
  group_by(period, subsector, day) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(day, nesting(period, subsector), fill = list(n = 0)) %>%
  arrange(period, subsector, day)
  
care_orga_mat = list(
  matrix(care_orga$n[care_orga$period == 0], ncol = nDays[1], byrow = T),
  matrix(care_orga$n[care_orga$period == 1], ncol = nDays[2], byrow = T),
  matrix(care_orga$n[care_orga$period == 2], ncol = nDays[3], byrow = T)
)

# Convert occupancy data into list
occupancy_data = split(occupancy_data[, c("period", "subsector", "id", "admission", "discharge")], 
                       f = occupancy_data$period)
names(occupancy_data) = NULL

# Select periods of interest in prevalence data
covid_data = covid_data %>%
  arrange(PATPER, DATE) %>%
  group_by(PATPER) %>%
  mutate(day = 1:n()) %>%
  filter(PATPER %in% periods[-2])

intubation_data = intubation_data %>%
  arrange(PATPER, DATE) %>%
  group_by(PATPER) %>%
  mutate(day = 1:n()) %>%
  filter(PATPER %in% periods[-2])

# Save prevalence data into a matrix to be used by the simulator
covid_data_mat = list(
  matrix(rep(covid_data$p_cov[covid_data$PATPER == periods[1]], nSS[1]), ncol = nDays[1], byrow = T),
  matrix(rep(covid_data$p_cov[covid_data$PATPER == periods[3]], nSS[2]), ncol = nDays[2], byrow = T),
  matrix(rep(covid_data$p_cov[covid_data$PATPER == periods[4]], nSS[3]), ncol = nDays[3], byrow = T)
)

intubation_data_mat = list(
  matrix(rep(intubation_data$p_intub[intubation_data$PATPER == periods[1]], nSS[1]), ncol = nDays[1], byrow = T),
  matrix(rep(intubation_data$p_intub[intubation_data$PATPER == periods[3]], nSS[2]), ncol = nDays[2], byrow = T),
  matrix(rep(intubation_data$p_intub[intubation_data$PATPER == periods[4]], nSS[3]), ncol = nDays[3], byrow = T)
)

# Model parameters
model_params_post_median = readRDS("results/simulated_ref/model_params_post_median.RData")

# Probability of acquisition
p_acq = vector("list", length(model_params_post_median))
names(p_acq) = names(model_params_post_median)

# Simulate colonization transmission
for (m in names(model_params_post_median)) {
  
  cat(paste0("---------- ", m, " ----------\n"))
  if (!dir.exists(paste0("data/simulated/", m))) dir.create(paste0("data/simulated/", m))
  p_acq_tmp = matrix(NA, ncol = length(occupancy_data), nrow = 100)
  
  for (j in 1:100) {
    output = vector("list", length(occupancy_data))
    for (i in 1:length(occupancy_data)) {
      p_acq_test = 0
      while(p_acq_test == 0) {
        output[[i]] = simulateCol(
          occupancy = occupancy_data[[i]],
          covid_data = covid_data_mat[[i]],
          intubation_data = intubation_data_mat[[i]],
          nBeds = care_orga_mat[[i]],
          period = i-1,
          intercept = model_params_post_median[[m]]["intercept1"],
          intercept2 = model_params_post_median[[m]]["intercept2"],
          intercept3 = model_params_post_median[[m]]["intercept3"], 
          intercept4 = model_params_post_median[[m]]["intercept4"],
          pcov2 = model_params_post_median[[m]]["pcov2"],
          pcov3 = model_params_post_median[[m]]["pcov3"],
          pcov4 = model_params_post_median[[m]]["pcov4"],
          pcov = model_params_post_median[[m]]["pcov"],
          pintub = model_params_post_median[[m]]["pintub"]
        ) %>%
          mutate(period = i-1) %>%
          arrange(period, subsector, id) %>%
          select(period, id, subsector, admission, discharge, detection_date, colType)
        
        # Percentage of positive individuals
        p_acq_test = sum(output[[i]]$colType == "acq") / sum(output[[i]]$colType %in% c("acq", "neg"))
      }
      p_acq_tmp[j,i] = p_acq_test
    }
  }
  
  p_acq[[m]] = p_acq_tmp
}

# Get mean and median estimates
p_acq = mapply(function(x, y) {out = data.frame(x) %>% rename(period1 = X1, period2 = X2, period3 = X3) %>% mutate(model = y)}, 
         p_acq, 
         names(p_acq), 
         SIMPLIFY = F) 

points_to_add = do.call("rbind", p_acq) %>%
  pivot_longer(contains("period"), names_to = "prob", values_to = "value") %>%
  group_by(model, prob) %>%
  summarise(med = median(value), me = mean(value), .groups = "drop")

# Plot
##############################################
## SIMULATE FAKE DATASET TO TEST CPP CODE
##############################################

rm(list = ls())
library(tidyverse)
library(Rcpp)

source("R/helperFunctions.R")
source("R/dictionaries.R")
source("R/modelParams.R")

# Compile Rcpp code
filesToRemove = list.files(
  path = "cpp/simulation/",
  pattern = "^.*\\.o$",
  full.names = T
)
if (length(filesToRemove)) sapply(filesToRemove, file.remove)
sourceCpp("cpp/simulations/Colonization_simulation.cpp")

# Load or simulate prevalence data 
real_prevalence_data = T

# Load data
covid_data = read.table("data/AMR_COVID-Data-230803/covid_prevalence_period.txt", header = T)
intubation_data = read.table("data/AMR_COVID-Data-230803/intubation_prevalence_period.txt", header = T)
occupancy_data = read.table("data/AMR_COVID-Data-230803/occupancy_data.txt", header = T)

# Information on period length and subsectors per period
nDays = occupancy_data %>% arrange(period) %>% group_by(period) %>% summarise(m = max(discharge)+1) %>% .$m
nSS = occupancy_data %>% arrange(period) %>% group_by(period) %>% summarise(m = length(unique(subsector))) %>% .$m
listSS = lapply(0:2, function(x) unique(occupancy_data$subsector[occupancy_data$period == x]))

# Number of beds per period 
care_orga = occupancy_data %>%
  select(period, subsector, id, admission, discharge) %>%
  group_by(period, subsector) %>%
  mutate(day = map2(admission, discharge, seq, by = 1)) %>%
  unnest(cols = day) %>%
  ungroup() %>%
  group_by(period, subsector, day) %>%
  summarise(n = n(), .groups = "drop") %>%
  complete(day, nesting(period, subsector), fill = list(n = 0)) %>%
  arrange(period, subsector, day)
  
care_orga_mat = list(
  matrix(care_orga$n[care_orga$period == 0], ncol = nDays[1], byrow = T),
  matrix(care_orga$n[care_orga$period == 1], ncol = nDays[2], byrow = T),
  matrix(care_orga$n[care_orga$period == 2], ncol = nDays[3], byrow = T)
)

# Convert occupancy data into list
occupancy_data = split(occupancy_data[, c("period", "subsector", "id", "admission", "discharge")], 
                       f = occupancy_data$period)
names(occupancy_data) = NULL

# Select periods of interest in prevalence data
covid_data = covid_data %>%
  arrange(PATPER, DATE) %>%
  group_by(PATPER) %>%
  mutate(day = 1:n()) %>%
  filter(PATPER %in% periods[-2])

intubation_data = intubation_data %>%
  arrange(PATPER, DATE) %>%
  group_by(PATPER) %>%
  mutate(day = 1:n()) %>%
  filter(PATPER %in% periods[-2])

# Save prevalence data into a matrix to be used by the simulator
covid_data_mat = list(
  matrix(rep(covid_data$p_cov[covid_data$PATPER == periods[1]], nSS[1]), ncol = nDays[1], byrow = T),
  matrix(rep(covid_data$p_cov[covid_data$PATPER == periods[3]], nSS[2]), ncol = nDays[2], byrow = T),
  matrix(rep(covid_data$p_cov[covid_data$PATPER == periods[4]], nSS[3]), ncol = nDays[3], byrow = T)
)

intubation_data_mat = list(
  matrix(rep(intubation_data$p_intub[intubation_data$PATPER == periods[1]], nSS[1]), ncol = nDays[1], byrow = T),
  matrix(rep(intubation_data$p_intub[intubation_data$PATPER == periods[3]], nSS[2]), ncol = nDays[2], byrow = T),
  matrix(rep(intubation_data$p_intub[intubation_data$PATPER == periods[4]], nSS[3]), ncol = nDays[3], byrow = T)
)

# Model parameters
model_params_post_median = readRDS("results/simulated_ref/model_params_post_median.RData")

# Probability of acquisition
p_acq = vector("list", length(model_params_post_median))
names(p_acq) = names(model_params_post_median)

# Simulate colonization transmission
for (m in names(model_params_post_median)) {
  
  cat(paste0("---------- ", m, " ----------\n"))
  if (!dir.exists(paste0("data/simulated/", m))) dir.create(paste0("data/simulated/", m))
  p_acq_tmp = matrix(NA, ncol = length(occupancy_data), nrow = 100)
  
  for (j in 1:100) {
    output = vector("list", length(occupancy_data))
    for (i in 1:length(occupancy_data)) {
      p_acq_test = 0
      while(p_acq_test == 0) {
        output[[i]] = simulateCol(
          occupancy = occupancy_data[[i]],
          covid_data = covid_data_mat[[i]],
          intubation_data = intubation_data_mat[[i]],
          nBeds = care_orga_mat[[i]],
          period = i-1,
          intercept = model_params_post_median[[m]]["intercept1"],
          intercept2 = model_params_post_median[[m]]["intercept2"],
          intercept3 = model_params_post_median[[m]]["intercept3"], 
          intercept4 = model_params_post_median[[m]]["intercept4"],
          pcov2 = model_params_post_median[[m]]["pcov2"],
          pcov3 = model_params_post_median[[m]]["pcov3"],
          pcov4 = model_params_post_median[[m]]["pcov4"],
          pcov = model_params_post_median[[m]]["pcov"],
          pintub = model_params_post_median[[m]]["pintub"]
        ) %>%
          mutate(period = i-1) %>%
          arrange(period, subsector, id) %>%
          select(period, id, subsector, admission, discharge, detection_date, colType)
        
        # Percentage of positive individuals
        p_acq_test = sum(output[[i]]$colType == "acq") / sum(output[[i]]$colType %in% c("acq", "neg"))
      }
      p_acq_tmp[j,i] = p_acq_test
      print(j)
    }
  }
  
  p_acq[[m]] = p_acq_tmp
}

# Get mean and median estimates
p_acq = mapply(function(x, y) {out = data.frame(x) %>% rename(period1 = X1, period2 = X2, period3 = X3) %>% mutate(model = y)}, 
         p_acq, 
         names(p_acq), 
         SIMPLIFY = F) 

points_to_add = do.call("rbind", p_acq) %>%
  pivot_longer(contains("period"), names_to = "period", values_to = "value") %>%
  group_by(model, period) %>%
  summarise(med = median(value), me = mean(value), .groups = "drop") %>%
  mutate(model = ifelse(model == "model2", "Model 2", "Model 1"),
         period = case_when(period == "period1" ~ "Period 1", 
                            period == "period2" ~ "Period 2",
                            period == "period3" ~ "Period 3")
         )
  
# Plot 
read.table("data/simulated/p_acq.txt", header = T, sep = "\t") %>%
  ggplot(., aes(x = period, y = p_acq)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  geom_point(data = points_to_add, aes(x = period, y = med), col = "red", shape = 8) +
  geom_point(data = points_to_add, aes(x = period, y = me), col = "cornflowerblue", shape = 18) +
  facet_grid(cols = vars(model)) +
  ylim(c(0,1)) +
  theme_bw() +
    theme(
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12),
      strip.text = element_text(size = 14)
    ) +
  labs(x = "", y = "Acquisition frequency\n(no. acquisitions / no. susceptibles)")
ggsave("figures/simulated/p_acquisition_models_post.png", height = 5, width = 8)


  



