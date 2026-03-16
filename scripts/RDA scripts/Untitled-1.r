sit <- Zenv_data$Name

varr <- Zenv_data %>% relocate(Name)


model_full <- rda(sit ~ ., data = varr)


model_test <- abund_all %>%
  MCMCglmm(
    data = .,
    fixed = sit ~ formula(varr),
    random = ~site,
    family = "poisson"
  )
summary(model_all)


bruh <- aov(formula(varr), data = varr)

summary(bruh)
