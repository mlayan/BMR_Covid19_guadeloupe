##############################################
## ANALYZE MCMC CHAINS FOR SIMULATED DATA
##############################################

rm(list = ls())
library(tidyverse)
library(coda)
library(gt)
source("R/helperFunctions.R")
source("R/modelParams.R")

##############################################
# Load mcmc chains
##############################################
estimates = data.frame()
p_accept = data.frame()
ess = data.frame()
correlations = data.frame()
post = data.frame()

for (m in names(model_params)) {
  for (f in list.files(paste0("results/simulated/", m), pattern = "mcmc_sim[0-9]+_1.txt$", full.names = T) ) {
    
    chainID = as.numeric(regmatches(f, gregexpr("(?<=_)[0-9]+(?=.txt)", f, perl = T))[[1]])
    nSim = as.numeric(regmatches(f, gregexpr("(?<=sim)[0-9]+(?=_)", f, perl = T))[[1]])
    
    tmp = try(read.table(f, header = T, sep=" "), silent = T)
    if (inherits(tmp, 'try-error')) next
    post_tmp = read.table(f, header = T, sep=" ")
    post_tmp$chain = chainID
    post_tmp$nSim = nSim
    post_tmp$model = m
    
    # Select columns that correspond to parameters in the model
    to_remove = colSums(post_tmp[, grepl("_p", colnames(post_tmp))]) == 0
    to_remove = names(which(to_remove))
    to_remove = c(to_remove, gsub("_p", "_a", to_remove), gsub("_p", "", to_remove))
    post_tmp = post_tmp[, !colnames(post_tmp) %in% to_remove]
    
    # Remove 10% burn-in
    post_tmp = post_tmp[ceiling(0.1*nrow(post_tmp)):nrow(post_tmp),]
    
    # Check if chain is stuck or not 
    if (var(post_tmp$logLik) == 0) {
      cat(paste0("Problem with file ", f, "\n"))
      next
    }
    
    # Store posterior distribution
    post = post_tmp %>% 
      select(!contains("_")) %>% 
      pivot_longer(-c(chain, iteration, nSim, model), names_to = "param", values_to = "value") %>% 
      bind_rows(., post)
    
    # Estimates of the posterior distribution
    estimates_tmp = post_tmp %>%
      select(!contains("_")) %>%
      pivot_longer(-c(chain, nSim, model), values_to = "value", names_to = "param") %>%
      group_by(chain, nSim, model, param) %>%
      summarise(m = median(value),
                q97_5 = quantile(value, 0.975),
                q2_5 = quantile(value,0.025), 
                .groups = "drop_last") %>%
      filter(!(m == 0 & q97_5 == 0 & q2_5 == 0))
    estimates = bind_rows(estimates, estimates_tmp)
    
    # Probability of acceptance
    p_accept_tmp = post_tmp %>%
      select(., matches("_[a,p]$")) %>%
      summarise(across(where(is.numeric), ~ sum(.x))) %>%
      pivot_longer(everything()) %>%
      mutate(param = gsub("_.*$", "", name), type = gsub("^.*_", "", name)) %>%
      select(param, type, value) %>%
      pivot_wider(values_from = value, names_from = type) %>%
      filter(p > 0) %>%
      mutate(prob = a/p) %>%
      select(param, prob) %>%
      mutate(chain = chainID, 
             nSim = nSim, 
             model = m
      )
    p_accept = bind_rows(p_accept, p_accept_tmp)
    
    # Chain mixing
    ess_tmp = post_tmp %>%
      select(-c(contains("_"), iteration, nSim, chain, model)) %>%
      summarise(across(where(is.numeric), ~effectiveSize(.x))) %>%
      mutate(
        chain = chainID, 
        nSim = nSim, 
        model = m
      )
    ess = bind_rows(ess, ess_tmp)
    
    # Correlations
    correlations_tmp = post_tmp %>%
      select(-c(contains("_"), chain, nSim, model, iteration)) %>%
      cor(.) %>% 
      as.data.frame(.) %>%
      rownames_to_column(var = "param") %>% 
      mutate(chain = chainID, 
             nSim = nSim, 
             model = m
      )
    
    correlations = bind_rows(correlations, correlations_tmp)
  }  
}

all_models = sort(unique(post$model))

##############################################
# Diagnostic plots
##############################################
# MCMC chains
for (mod in all_models) {
  p = post %>%
    filter(nSim <= 8, model == mod) %>%
    ggplot(., aes(x = iteration, y = value, col = factor(chain))) +
    geom_line() +
    facet_grid(cols = vars(nSim), rows = vars(param), scales = "free_y") +
    theme_classic() +
    labs(x = "Iteration", y = "Posterior density", col = "Chain", title = mod)
  print(p)
  ggsave(paste0("figures/simulated/mixing_", mod, ".png"),
         height = 20, width = 30, units = "cm")
}

# Acceptance rates
p_accept %>%
  ggplot(., aes(x = param, y = prob, col = factor(chain))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter() +
  theme_bw() +
  facet_grid(rows = vars(model)) +
  labs(x = "", y = "Acceptance probability", col = "Chain") +
  ylim(c(0,1))
ggsave("figures/simulated/p_accept.png",
       height = 20, width = 20, units = "cm")

# ESS values
ess %>%
  pivot_longer(-c(chain, nSim, model), names_to = "param", values_to = "ess") %>%
  ggplot(., aes(x = param, y = ess, col = factor(chain))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter() +
  theme_bw() +
  facet_grid(rows = vars(model)) +
  labs(x = "", y = "ESS", col = "Chain")

# Correlations across all simulations
for (mod in all_models) {
  p = pairs(
    post %>%
      filter(nSim == 1, chain == 1, model == mod, param != "logLik") %>%
      pivot_wider(names_from = param, values_from = value) %>%
      select(-c(iteration, chain, nSim, model)),
    upper.panel = panel.cor,
    diag.panel = panel.hist,
    main = mod
  )
  print(p)
  png(paste0("figures/simulated/corr_", mod, ".png"),  
      width = 20, height = 20, units = "cm", res = 300)
  pairs(
    post %>%
      filter(nSim == 1, chain == 1, model == mod, param != "logLik") %>%
      pivot_wider(names_from = param, values_from = value) %>%
      select(-c(iteration, chain, nSim, model)),
    upper.panel = panel.cor,
    diag.panel = panel.hist,
    main = mod
  )
  dev.off()
}

##############################################
# Comparison with true values
##############################################
# Median + 95% CrI compared to true value for each simulation and model 
for (mod in all_models) {
  true_values = data.frame(
    param = names(model_params[[mod]]),
    true_val = as.numeric(model_params[[mod]]),
    model = rep(mod, length(model_params[[mod]]))
  ) %>%
    filter(true_val != 0)
  
  p = estimates %>%
    filter(model == mod, !param %in% c("iteration", "logLik")) %>%
    group_by(param) %>%
    arrange(m) %>%
    mutate(new_order = 1:n()) %>%
    ungroup() %>%
    ggplot(., aes(x = new_order, y = m, ymin = q2_5, ymax = q97_5)) +
    facet_grid(rows = vars(param), scales = "free_y") +
    geom_pointrange(fatten = 2) +
    geom_hline(data = true_values, aes(yintercept = true_val)) +
    theme_bw() +
    labs(x = "Simulation", y = "Median (95% CrI)", title = mod)
  print(p)
  ggsave(paste0("figures/simulated/estimation_tue_val_", mod, ".png"),
         height = 15, width = 20, units = "cm")
}



# Comparison metrics
true_values = data.frame(
  param = as.character(unlist(lapply(model_params, names))),
  true_val = as.numeric(unlist(model_params)),
  model = rep(names(model_params), each = length(model_params[[1]]))
) %>%
  filter(true_val != 0)


median_estimate = post %>%
  filter(param != "logLik") %>%
  group_by(model, param) %>%
  summarise(m1 = round(median(value), 3), .groups = "drop")

to_save = estimates %>%
  filter(!param %in% c("iteration", "logLik")) %>%
  left_join(., true_values, by = c("param", "model")) %>%
  group_by(param, model) %>%
  summarise(
    true_value = unique(true_val),
    calibration = round(mean(true_val <= q97_5 & true_val >= q2_5)*100),
    mrb = round(mean((m-true_val)/true_val), 3),
    relativewidth = round(abs(mean((q97_5 - q2_5)/true_val)), 3),
    .groups = "drop"
  ) %>%
  left_join(., median_estimate, by = c("param", "model")) %>%
  arrange(model, param) %>%
  select(model, param, true_value, m1, mrb, calibration, relativewidth) %>%
  group_by(model) %>%
  gt(.) %>%
  tab_header(title = "Model performances without data augmentation") %>%
  cols_label(
    param = "Parameter",
    true_value = "True value",
    m1 = "Estimated median",
    mrb = "Mean relative bias",
    calibration = "Calibration (%)",
    relativewidth = "Relative 95% CrI width"
  )
to_save

gtsave(to_save, filename = "model_performances.png", 
       path = "tables/simulated/")


