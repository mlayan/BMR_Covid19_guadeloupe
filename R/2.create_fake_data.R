##############################################
## CREATE FAKE DATASET TO TEST CPP CODE
##############################################

rm(list = ls())
library(tidyverse)

nDays = 5
nPeriod = 2
nSubjectPerPeriod = 2

# Create episode data
set.seed(101)
episode_data = data.frame(
  period = rep(0:(nPeriod-1), each = nSubjectPerPeriod*nDays),
  id = rep(rep(0:(nSubjectPerPeriod-1), each = nDays), nPeriod),
  subsector = 0,
  coltype = rep(sample(c("acq", "neg", "pos"), nSubjectPerPeriod*nPeriod, replace = T), each = nDays)
)

episodes = expand.grid(period = 0:(nPeriod-1), id = 0:(nSubjectPerPeriod-1))
episode_data$status = unlist(apply(episodes, 1, function(x) {
  out = rep(0, nDays)
  coltype = unique(episode_data$coltype[episode_data$id == x[1] & episode_data$period == x[2]])
  if (coltype == "pos") out = rep(1, nDays)
  if (coltype == "acq") out[sample(2:nDays, 1):nDays] = 1
  return(out)
}, simplify = F))

episode_data$status[c(1,2,20)] = -1

episode_data %>%
  write.table(., "data/episode_data.txt", quote = F, sep = " ", row.names = F, col.names = F)

# Create COVID-19 prevalence data
data.frame(
  period = rep(0:(nPeriod-1), each = nDays),
  subsector = 0,
  t = rep(0:(nDays-1), nPeriod),
  val = c(rep(0, nDays), 0.26 + rnorm(nDays, 0, 0.1))
) %>% 
  mutate(val = round(val, 3)) %>%
  write.table(., "data/covid_data.txt", quote = F, sep = " ", row.names = F, col.names = F)


# Create dialysis prevalence data
data.frame(
  period = rep(0:(nPeriod-1), each = nDays),
  subsector = 0,
  t = rep(0:(nDays-1), nPeriod),
  val = 0.05 + rnorm(nDays*nPeriod, 0, 0.01)
) %>% 
  mutate(val = round(val, 3)) %>%
  write.table(., "data/dialysis_data.txt", quote = F, sep = " ", row.names = F, col.names = F)

# Create intubation prevalence data
data.frame(
  period = rep(0:(nPeriod-1), each = nDays),
  subsector = 0,
  t = rep(0:(nDays-1), nPeriod),
  val = c(0.60 + rnorm(nDays, 0, 0.1), 0.70 + rnorm(nDays, 0, 0.1))
) %>% 
  mutate(val = round(val, 3)) %>%
  write.table(., "data/intubation_data.txt", quote = F, sep = " ", row.names = F, col.names = F)
