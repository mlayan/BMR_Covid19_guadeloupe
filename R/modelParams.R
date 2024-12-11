# Ref pre-pandemic
model_params = list(
  "model1" = c("intercept1" = -2.7, "intercept2" = 1, "intercept3" = 0.5, "intercept4" = 0, "pcov2" = 0, "pcov3" = 0, "pcov4" = 0, "pcov" = 0),
  "model2" = c("intercept1" = -2.7, "intercept2" = 0, "intercept3" = 0, "intercept4" = 0, "pcov2" = 1.3, "pcov3" = 3, "pcov4" = 0, "pcov" = 0)
)


# Ref 4th wave
model_params_ref2 = list(
  "model1" = c("intercept1" = -1, "intercept2" = -1.7, "intercept3" = -0.5),
  "model2" = c("intercept2" = -2.7, "pcov2" = 1.3, "pcov3" = 3)
)


# Independent inference
model_params_indpd = list(
  "model3" = c("intercept1" = -1.8, "intercept2" = -1.4, "intercept3" = -1.6, "intercept4" = 0, "pcov2" = 0, "pcov3" = 0, "pcov4" = 0, "pcov" = 0.1, "pintub" = 0),
  "model2" = c("intercept1" = -1.8, "intercept2" = -1.8, "intercept3" = -1.8, "intercept4" = 0, "pcov2" = 0.6, "pcov3" = 0.2, "pcov4" = 0, "pcov" = 0, "pintub" = 0),
  "model1" = c("intercept1" = -2, "intercept2" = -1.4, "intercept3" = -1.6, "intercept4" = 0, "pcov2" = 0, "pcov3" = 0, "pcov4" = 0, "pcov" = 0, "pintub" = 0)
)
