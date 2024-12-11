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

for (m in names(model_params_ref2)) {
  for (d in c("yes", "no")) {
    ddir = ifelse(d=="yes", "simulated", "simulated_fixed")
    for (f in list.files(paste0("results/", ddir, "/", m), pattern = "mcmc_sim[0-9]+_1.txt$", full.names = T) ) {
      
      chainID = as.numeric(regmatches(f, gregexpr("(?<=_)[0-9]+(?=.txt)", f, perl = T))[[1]])
      nSim = as.numeric(regmatches(f, gregexpr("(?<=sim)[0-9]+(?=_)", f, perl = T))[[1]])
      
      tmp = try(read.table(f, header = T, sep=" "), silent = T)
      if (inherits(tmp, 'try-error')) next
      post_tmp = read.table(f, header = T, sep=" ")
      post_tmp$chain = chainID
      post_tmp$nSim = nSim
      post_tmp$model = m
      post_tmp$augmentation = d
      
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
        pivot_longer(-c(chain, iteration, nSim, model, augmentation), names_to = "param", values_to = "value") %>% 
        bind_rows(., post)
      
      # Estimates of the posterior distribution
      estimates_tmp = post_tmp %>%
        select(!contains("_")) %>%
        pivot_longer(-c(chain, nSim, model, augmentation), values_to = "value", names_to = "param") %>%
        group_by(chain, nSim, model, augmentation, param) %>%
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
               model = m,
               augmentation = d
        )
      p_accept = bind_rows(p_accept, p_accept_tmp)
      
      # Chain mixing
      ess_tmp = post_tmp %>%
        select(-c(contains("_"), iteration, nSim, chain, model)) %>%
        summarise(across(where(is.numeric), ~effectiveSize(.x))) %>%
        mutate(
          chain = chainID, 
          nSim = nSim, 
          model = m,
          augmentation = d
        )
      ess = bind_rows(ess, ess_tmp)
      
      # Correlations
      correlations_tmp = post_tmp %>%
        select(-c(contains("_"), chain, nSim, model, iteration, augmentation)) %>%
        cor(.) %>% 
        as.data.frame(.) %>%
        rownames_to_column(var = "param") %>% 
        mutate(chain = chainID, 
               nSim = nSim, 
               model = m,
               augmentation = d
        )
      
      correlations = bind_rows(correlations, correlations_tmp)
    }   
  }
}

all_models = sort(unique(post$model))

##############################################
# Diagnostic plots
##############################################
# MCMC chains
for (mod in all_models) {
  for (a in c("yes", "no")) {
    p = post %>%
      filter(nSim <= 8, model == mod, augmentation == a) %>%
      ggplot(., aes(x = iteration, y = value, col = factor(chain))) +
      geom_line() +
      facet_grid(cols = vars(nSim), rows = vars(param), scales = "free") +
      theme_classic() +
      labs(x = "Iteration", y = "Posterior density", col = "Chain", title = paste(mod, a))
    print(p)
    ggsave(paste0("figures/simulated/mixing_", mod, "_", a, ".png"),
           height = 20, width = 30, units = "cm")
  }
}

# Acceptance rates
p_accept %>%
  ggplot(., aes(x = param, y = prob, col = factor(chain))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter() +
  theme_bw() +
  facet_grid(rows = vars(interaction(model, augmentation))) +
  labs(x = "", y = "Acceptance probability", col = "Chain") +
  ylim(c(0,1))
ggsave("figures/simulated/p_accept.png",
       height = 20, width = 20, units = "cm")

# ESS values
ess %>%
  pivot_longer(-c(chain, nSim, model, augmentation), names_to = "param", values_to = "ess") %>%
  ggplot(., aes(x = param, y = ess, col = factor(chain))) +
  geom_boxplot(outlier.shape = NA) + 
  geom_jitter() +
  theme_bw() +
  facet_grid(rows = vars(interaction(model, augmentation))) +
  labs(x = "", y = "ESS", col = "Chain")

# Correlations across all simulations
for (mod in all_models) {
  for (a in c("yes", "no")) {
    p = pairs(
      post %>%
        filter(nSim == 1, chain == 1, model == mod, augmentation == a, param != "logLik") %>%
        pivot_wider(names_from = param, values_from = value) %>%
        select(-c(iteration, chain, nSim, model, augmentation)),
      upper.panel = panel.cor,
      diag.panel = panel.hist,
      main = paste(mod, a)
    )
    print(p)
    png(paste0("figures/simulated/corr_", mod, "_", a, ".png"),
        width = 20, height = 20, units = "cm", res = 300)
    pairs(
      post %>%
        filter(nSim == 1, chain == 1, model == mod, augmentation == a, param != "logLik") %>%
        pivot_wider(names_from = param, values_from = value) %>%
        select(-c(iteration, chain, nSim, model, augmentation)),
      upper.panel = panel.cor,
      diag.panel = panel.hist,
      main = paste(mod, a)
    )
    dev.off()
  }
}

# Comparison prior and posterior
to_save = post %>%
  filter(nSim == 1, chain == 1, param != "logLik") %>%
  ggplot(., aes(x = value)) +
  stat_function(aes(col = "prior"), fun = function(x) {dnorm(x, 0, 10)}) +
  geom_density(aes(col = "posterior")) +
  theme_light() +
  facet_grid(rows = vars(param), cols = vars(interaction(model, augmentation))) +
  scale_color_manual(values = c("posterior" = "orchid4",
                                "prior" = "orange")) +
  theme(plot.title = element_text(hjust = 0.5, size = 12),
        axis.title.y = element_text(size = 10)) +
  labs(x = "", y = "Density", col = "")
ggsave("figures/simulated/prior_posterior.png",
       to_save,
       height = 20, width = 20, units = "cm")

##############################################
# Comparison with true values
##############################################
# Median + 95% CrI compared to true value for each simulation and model 
for (mod in all_models) {
  for (a in c("yes", "no")) {
    true_values = data.frame(
      param = names(model_params_ref2[[mod]]),
      true_val = as.numeric(model_params_ref2[[mod]]),
      model = rep(mod, length(model_params_ref2[[mod]]))
    ) %>%
      filter(true_val != 0) %>%
      mutate(param = case_when(
        param == "intercept1" ~ "β1",
        param == "intercept2" & model == "model2" ~ "β1",
        param == "intercept2" & model == "model1" ~ "β2",
        param == "intercept3" ~ "β3",
        param == "pcov2" ~ "γ2",
        param == "pcov3" ~ "γ3"
      ))
    
    t = ifelse(mod == "model2", "Model 2", "Model 1")
    
    p = estimates %>%
      filter(model == mod, augmentation == a, !param %in% c("iteration", "logLik")) %>%
      group_by(model, augmentation, param) %>%
      arrange(m) %>%
      mutate(new_order = 1:n()) %>%
      ungroup() %>%
      mutate(param = case_when(
        param == "intercept1" ~ "β1",
        param == "intercept2" & model == "model2" ~ "β1",
        param == "intercept2" & model == "model1" ~ "β2",
        param == "intercept3" ~ "β3",
        param == "pcov2" ~ "γ2",
        param == "pcov3" ~ "γ3"
      )) %>%
      ggplot(., aes(x = new_order, y = m, ymin = q2_5, ymax = q97_5)) +
      facet_grid(rows = vars(param), scales = "free_y") +
      geom_pointrange(fatten = 2) +
      geom_hline(data = true_values, aes(yintercept = true_val)) +
      theme_bw() +
      theme(
        axis.title = element_text(size = 16), 
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 14),
        legend.text = element_text(size = 14),
        plot.title = element_text(size = 20, hjust = 0.5)
      ) +
      labs(x = "Simulation", y = "Median (95% CrI)", title = paste(t, a))
    print(p)
    ggsave(paste0("figures/simulated/estimation_tue_val_", mod, "_", a, ".png"),
           height = 15, width = 20, units = "cm")
  }
}


# Comparison metrics
true_values = data.frame(
  param = as.character(unlist(lapply(model_params_ref2, names))),
  true_val = as.numeric(unlist(model_params_ref2)),
  model = rep(names(model_params_ref2), each = length(model_params_ref2[[1]]))
) %>%
  filter(true_val != 0)


median_estimate = post %>%
  filter(param != "logLik") %>%
  group_by(model, param, augmentation) %>%
  summarise(m1 = round(median(value), 3), .groups = "drop")
  
to_save = estimates %>%
  filter(!param %in% c("iteration", "logLik")) %>%
  left_join(., true_values, by = c("param", "model")) %>%
  group_by(param, model, augmentation) %>%
  summarise(
    true_value = unique(true_val),
    calibration = round(mean(true_val <= q97_5 & true_val >= q2_5)*100),
    mrb = round(mean((m-true_val)/true_val)*100),
    relativewidth = round(abs(mean((q97_5 - q2_5)/true_val))*100),
    .groups = "drop"
  ) %>%
  left_join(., median_estimate, by = c("param", "model", "augmentation")) %>%
  arrange(model, param) %>%
  select(model, param, augmentation, true_value, m1, mrb, calibration, relativewidth) %>%
  mutate(
    param = case_when(
      param == "intercept2" & model == "model2" ~ "β1",
      param == "intercept1" ~ "β1",
      param == "intercept2" & model == "model1" ~ "β2",
      param == "intercept3" ~ "β3",
      param == "pcov2" ~ "γ2",
      param == "pcov3" ~ "γ3"
        ) 
    ) %>%
  mutate(model = case_when(model == "model2" ~ "Model 2", model == "model1" ~ "Model 1"))%>%
  arrange(model) %>%
  group_by(model, augmentation) %>%
  gt(.) %>%
  tab_header(title = "Model performances") %>%
  cols_label(
    param = "Parameters",
    true_value = "True value",
    m1 = "Estimated median",
    mrb = "Mean relative bias (%)",
    calibration = "Calibration (%)",
    relativewidth = "Relative 95% CrI width (%)"
  )
to_save

gtsave(to_save, filename = "model_performances.png",
       path = "tables/simulated/")

##############################################
# Save pooled posterior medians
##############################################
m1 = post %>%
  filter(model == "model1", param != "logLik") %>%
  select(iteration, nSim, param, value) %>%
  pivot_wider(names_from = param, values_from = value)

m2 = post %>%
  filter(model == "model2", param != "logLik") %>%
  select(iteration, nSim, param, value) %>%
  pivot_wider(names_from = param, values_from = value)

model_params_post_median = list(
  model2 = c(
    "intercept1" = median(m2$intercept2), 
    "intercept2" = 0, 
    "intercept3" = 0, 
    "intercept4" = 0, 
    "pcov2" = median(m2$pcov2), 
    "pcov3" = median(m2$pcov3), 
    "pcov4" = 0,
    "pcov" = 0
  ),
  model1 = c(
    "intercept1" = median(m1$intercept1 + m1$intercept2), 
    "intercept2" = median(-m1$intercept1), 
    "intercept3" = median(m1$intercept3 - m1$intercept1), 
    "intercept4" = 0, 
    "pcov2" = 0, 
    "pcov3" = 0, 
    "pcov4" = 0,
    "pcov" = 0
  )
)

saveRDS(object = model_params_post_median, 
        file = "results/simulated_ref/model_params_post_median.RData")
