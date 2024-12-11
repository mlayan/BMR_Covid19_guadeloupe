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
# set.seed(20231130)

# Load data
covid_data = read.table("data/AMR_COVID-Data-230803/covid_prevalence_period.txt", header = T)
intubation_data = read.table("data/AMR_COVID-Data-230803/intubation_prevalence_period.txt", header = T)
occupancy_data = read.table("data/AMR_COVID-Data-230803/occupancy_data.txt", header = T)
microbio_data = read.delim("data/AMR_COVID-Data-230803/AMRCOVID_MICROBIO_20230803_151902.txt",
                           sep = "#", comment.char="", header = T, encoding = "latin1")

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

care_orga %>%
  arrange(period, subsector, day) %>%
  select(period, subsector, day, n) %>%
  write.table(., "data/nbeds_data_sim.txt", quote = F, sep = " ", col.names = F, row.names = F)

# Convert occupancy data into list
tosplit = occupancy_data[occupancy_data$SUBJID %in% microbio_data$SUBJID, 
                         c("period", "subsector", "id", "admission", "discharge")] %>%
  arrange(period, id) %>%
  group_by(period) %>%
  mutate(id = 0:(n()-1)) %>%
  ungroup()
occupancy_data = split(tosplit, f = tosplit$period)
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

# Save prevalence data for inference
covid_data_inference = data.frame()
periods_selected = periods[-2]
for (p in 0:2) {
  for (s in listSS[[p+1]]) {
    covid_data_tmp = covid_data %>%
      ungroup() %>%
      filter(PATPER == periods_selected[p+1]) %>%
      mutate(period = p, day = day-1, subsector = s) %>%
      select(period, subsector, day, p_cov)
    covid_data_inference = bind_rows(covid_data_inference, covid_data_tmp)
  }
}

covid_data_inference %>%
  arrange(period, subsector, day) %>%
  write.table(., "data/covid_data_sim.txt", quote = F, sep = " ", col.names = F, row.names = F)

intubation_data_inference = data.frame()
for (p in 0:2) {
  for (s in listSS[[p+1]]) {
    intubation_data_tmp = intubation_data %>%
      ungroup() %>%
      filter(PATPER == periods_selected[p+1]) %>%
      mutate(period = p, day = day-1, subsector = s) %>%
      select(period, subsector, day, p_intub)
    intubation_data_inference = bind_rows(intubation_data_inference, intubation_data_tmp)
  }
}

intubation_data_inference %>%
  arrange(period, subsector, day) %>%
  write.table(., "data/intubation_data_sim.txt", quote = F, sep = " ", col.names = F, row.names = F)

# Probability of acquisition
p_acq = vector("list", length(model_params))
names(p_acq) = names(model_params)

# Simulate colonization transmission
for (m in names(model_params)) {
  
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
          nBeds = care_orga_mat[[i]],
          period = i-1,
          intercept = model_params[[m]]["intercept1"],
          intercept2 = model_params[[m]]["intercept2"],
          intercept3 = model_params[[m]]["intercept3"], 
          intercept4 = model_params[[m]]["intercept4"],
          pcov2 = model_params[[m]]["pcov2"],
          pcov3 = model_params[[m]]["pcov3"],
          pcov4 = model_params[[m]]["pcov4"],
          pcov = model_params[[m]]["pcov"]
        ) %>%
          mutate(period = i-1) %>%
          arrange(period, subsector, id) %>%
          mutate(
            last_negative = case_when(
              colType %in% c("pos", "neg") ~ -1,
              colType == "acq" & (detection_date - admission) %% 7 == 0 & admission+7*((detection_date-admission)%/%7-1) < admission ~ admission,
              colType == "acq" & (detection_date - admission) %% 7 == 0 & admission+7*((detection_date-admission)%/%7-1) >= admission ~ admission+7*((detection_date-admission)%/%7-1),
              colType == "acq" &  (detection_date - admission) %% 7 != 0 ~ admission+7*((detection_date-admission)%/%7)
              ),
            first_positive = case_when(
              colType == "neg" ~ -1,
              colType == "pos" ~ detection_date,
              colType == "acq" & (detection_date - admission) %% 7 == 0 ~ detection_date,
              colType == "acq" & (detection_date - admission) %% 7 != 0 & admission+7*((detection_date-admission)%/%7+1) > discharge ~ discharge,
              colType == "acq" & (detection_date - admission) %% 7 != 0 & admission+7*((detection_date-admission)%/%7+1) <= discharge ~ admission+7*((detection_date-admission)%/%7+1),
            )
            )
        
        # Percentage of positive individuals
        # print(sum(output[[i]]$colType == "pos") / nrow(output[[i]]))
        # print(sum(output[[i]]$colType == "acq") / nrow(output[[i]]))
        p_acq_test = sum(output[[i]]$colType == "acq") / sum(output[[i]]$colType %in% c("acq", "neg"))
      }
      p_acq_tmp[j,i] = p_acq_test
    }
    
    # Save for inference
    do.call("rbind", output) %>%
      select(period, id, subsector, admission, discharge, detection_date, colType) %>%
      arrange(period, id) %>%
      write.table(., paste0("data/simulated_fixed/", m, "/episode_data_sim", j, ".txt"),
                  quote = F, sep = " ", col.names = F, row.names = F)
    
    do.call("rbind", output) %>%
      select(period, id, subsector, admission, discharge, first_positive, last_negative, colType) %>%
      arrange(period, id) %>%
      write.table(., paste0("data/simulated/", m, "/episode_data_sim", j, ".txt"),
                  quote = F, sep = " ", col.names = F, row.names = F)
    cat(paste0("Simulation ", j, "\n"))
  }
  
  p_acq[[m]] = p_acq_tmp
}

# Plot p_acq per scenario
p_acq = mapply(function(x, y) {out = data.frame(x) %>% rename(period1 = X1, period2 = X2, period3 = X3) %>% mutate(model = y)}, 
         p_acq, 
         names(p_acq), 
         SIMPLIFY = F)

do.call("rbind", p_acq) %>%
  pivot_longer(-model, names_to = "period", values_to = "p_acq") %>%
  ggplot(., aes(x = period, y = p_acq)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  facet_grid(cols = vars(model)) +
  geom_hline(yintercept = 0.06, col = "orange") +
  geom_hline(yintercept = 0.15, col = "red") +
  ylim(c(0,1)) +
  theme_bw() +
  labs(x = "", y = "Probability of acquisition\n(acquisitions/susceptibles)")
ggsave("figures/simulated/p_acquisition_models.png", height = 4, width = 10)


do.call("rbind", p_acq) %>%
  pivot_longer(-model, names_to = "period", values_to = "p_acq") %>%
  mutate(model = case_when(model == "model1" ~ "Model 1",
                           model == "model2" ~ "Model 2"),
         period = case_when(period == "period1" ~ "Period 1",
                            period == "period2" ~ "Period 2",
                            period == "period3" ~ "Period 3")
         ) %>%
  ggplot(., aes(x = period, y = p_acq)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter() +
  facet_grid(cols = vars(model)) +
  ylim(c(0,1)) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 14)
  ) +
  labs(x = "", y = "Probability of acquisition\n(acquisitions/susceptibles)")
ggsave("figures/simulated/p_acquisition_poster.png", height = 5, width = 8)


do.call("rbind", p_acq) %>%
  pivot_longer(-model, names_to = "period", values_to = "p_acq") %>%
  mutate(model = case_when(model == "model1" ~ "Model 1", 
                           model == "model2" ~ "Model 2"), 
         period = case_when(period == "period1" ~ "Period 1", 
                            period == "period2" ~ "Period 2", 
                            period == "period3" ~ "Period 3")
  ) %>%
  write.table("data/simulated/p_acq.txt", row.names = F, quote = F, sep = "\t")

