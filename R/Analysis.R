pacman::p_load("dplyr","lme4", "MASS", "brms", "MCMCglmm", "quantreg","lmerTest", "latex2exp", "tidybayes", "bayesplot", "rstanarm", "plotrix", "emmeans", "ggExtra", "marginaleffects", "gridExtra", "bayestestR","performance", "cowplot", "patchwork", "png", "magick", "ggplot2")

# import data (be sure to change Y/N to 0/1)
data <- read.csv("Data/Lizard_Db_final.csv", head = T) %>% 
  mutate(temp = as.character(temp))

# separate the two species into df
guich_df <-data %>% filter(species == "guichenoti") 
deli_df <- data %>% filter(species == "delicata")

# look at developmental time differences by species and temp
data_hatch <- data %>% 
  group_by(species, temp) %>% 
  summarise(dev_time_mean = mean(dev_time, na.rm = TRUE),
            dev_sd = sd(dev_time, na.rm = TRUE))

###########
# Functions
###########
## pMCMC Function
# Calculates the, p-value or pMCMC value for a posterior distribution
# x The vector for the posterior distribution. Note that this will test the null hypothesis that the parameter of interest is significantly different from 0. 
pmcmc <- function(x){
  2*(1 - max(table(x<0) / nrow(x)))
}



###############################################
############## Guichenoti #####################
###############################################
######
# Guichenoti Survival Analysis
######
# Guich survival model
Guich_m1_brms <- brm(mortality ~ yolk + temp+ temp:yolk  + 
                       scale(egg_mass_final) + # egg mass scaled
                       (1|clutch),  # clutch random factor
                     family = bernoulli(link = "logit"), 
                     data = guich_df, 
                     iter= 5000, warmup = 1000, 
                     thin = 4, cores = 8)


####################
# Model 1 checks: lags, residuals, r2, summary
####################
# checking lags for this model 
draws <- as.array(Guich_m1_brms)
mcmc_acf(draws,  pars = c("b_Intercept","b_yolkC", "b_temp28", "b_yolkC:temp28", "b_scaleegg_mass_final"), lags =10)
# plots
plot(Guich_m1_brms)
#R2 and summary of full model
bayes_R2(Guich_m1_brms)
summary(Guich_m1_brms)

# stan plot
stanplot(Guich_m1_brms, 
         type = "areas",
         prob = 0.95)

# overall egg mass on survival
plot(conditional_effects(Guich_m1_brms, "egg_mass_final"), ask = FALSE)

# temperature X yolk treatment interaction
plot(conditional_effects(Guich_m1_brms, "temp:yolk", points = TRUE))

# egg mass by yolk treatment on survival
plot(conditional_effects(Guich_m1_brms, "egg_mass_final:yolk"), ask = FALSE)


# extracting posteriors
post_Guich_m1 <- posterior_samples(Guich_m1_brms, pars = "^b")
variable.names(post_Guich_m1)



###############################################
############## Delicata #####################
###############################################
######
# Deli Survival Analysis
######
# Deli survival model
Deli_m1_brms <- brm(mortality ~ yolk + temp+ temp:yolk  + 
                       scale(egg_mass_final) + # egg mass scaled
                       (1|clutch),  # clutch random factor
                     family = bernoulli(link = "logit"), 
                     data = deli_df, 
                     iter= 5000, warmup = 1000, 
                     thin = 4, cores = 8)


####################
# Model 1 checks: lags, residuals, r2, summary
####################
# checking lags for this model 
draws <- as.array(Deli_m1_brms)
mcmc_acf(draws,  pars = c("b_Intercept","b_yolkC", "b_temp28", "b_yolkC:temp28", "b_scaleegg_mass_final"), lags =10)
# plots
plot(Deli_m1_brms)
#R2 and summary of full model
bayes_R2(Deli_m1_brms)
summary(Deli_m1_brms)

# stan plot
stanplot(Deli_m1_brms, 
         type = "areas",
         prob = 0.95)

# overall egg mass on survival
plot(conditional_effects(Deli_m1_brms, "egg_mass_final"), ask = FALSE)

# temperature X yolk treatment interaction
plot(conditional_effects(Deli_m1_brms, "temp:yolk", points = TRUE))

# egg mass by yolk treatment on survival
plot(conditional_effects(Deli_m1_brms, "egg_mass_final:yolk"), ask = FALSE)


# extracting posteriors
post_Deli_m1 <- posterior_samples(Deli_m1_brms, pars = "^b")
variable.names(post_Deli_m1)






yolk_a <- as.array(post_Guich_m1[,"b_Intercept"])
yolk_c <- as.array(yolk_a +  post_Guich_m1[,"b_yolkC"])

Yolk_survival <- cbind(yolk_a, yolk_c)
mcmc_areas(Yolk_survival, 
           pars = c("yolk_a", "yolk_c"),
           prob = 0.95, 
           prob_outer = 0.99, 
           point_est = "mean")+
  theme_bw() +
  theme(axis.text = element_text(size=12)) +
  theme(legend.title = element_text(colour="white", size = 16, face='bold')) +
  labs(y = TeX("Sex class"), x = "Slope Differences") 














########
# old code
#######








####################  
# extract posteriors for plotting 
####################
summary(Guich_m1_brms)
post_Guich_m1_brms <- posterior_samples(Guich_m1_brms, pars = "^b")
variable.names(post_Guich_m1_brms)

## extracting posteriors for interaction of sex and mass
XXf.mass.posterior <- as.array(post_Bas_m1[,"b_logMass"])
XXm.mass.posterior <- as.array(XXf.mass.posterior +  post_Bas_m1[,"b_sexXXm:logMass"])
XYm.mass.posterior <- as.array(XXf.mass.posterior +  post_Bas_m1[,"b_sexXYm:logMass"])





guich_model <- glmer(mortality ~temp + yolk +  temp:yolk+ (1|Egg_ID), data = guich_df, family = "binomial") 
summary(guich_model)
car::Anova(guich_model, type=3)
check_model(guich_model)
# df for figure
guich_dat <- emmeans(guich_model, specs = c("temp", "yolk"), type = "response") %>%
  summary() %>%
  as.data.frame()
# percent diff for figure
guich_23 <- round((0.3066667 - 0.1509434),digits = 2)
guich_28 <- round((0.3424657 - 0.1200000),digits = 2)


######
# Guichenoti offspring mass on live subjects - no interaction; overall yolk effect
######
guich_alive <- guich_df %>% filter(liz_status == "ALIVE") 
Guich_mass_model <- lm(log(hatch_mass)~ temp + yolk +  temp:yolk, data = guich_alive)
# Checks and plots
summary(Guich_mass_model) 
anova(Guich_mass_model)
check_model(Guich_mass_model)
Guich_mass_emm <- emmeans(Guich_mass_model, specs = c("yolk"), type = "response") %>% # back transform from log scale
  summary() %>%
  as.data.frame()
plot(emmeans(Guich_mass_model, specs = c("yolk"), type = "response"))


######
# Guichenoti offspring SVL on live subjects - no interaction; no differences
######
Guich_svl_model <- lm(log(hatch_svl)~ temp + yolk +  temp:yolk, data = guich_alive)
# Checks and plots
summary(Guich_svl_model) 
anova(Guich_svl_model)
check_model(Guich_svl_model)


######
# Guichenoti offspring BCI (OLS residual difference from svl and mass
######
guich_bci <- guich_alive %>% # filter out NA's to keep residuals from lm
  filter_at(vars(hatch_svl, hatch_mass), all_vars(!is.na(.)))
# calcualting residuals from OLS Regression of mass and SVL
guich_fit <- lm(log(hatch_svl) ~ log(hatch_mass), data = guich_bci)
check_model(guich_fit)
guich_bci$BCI <- residuals(guich_fit)

# BCI test on yolk and temp - residuals are from ols of svl and mass- nothing
Guich_bci_model <- lm(BCI~ temp + yolk +  temp:yolk, data = guich_bci)
# Checks and plots
summary(Guich_bci_model) 
anova(Guich_bci_model)
check_model(Guich_bci_model)



########
# Guich figures
########
legend_title <- "Resource Treatment"
pd = position_dodge(.8)

Guich_fig <-ggplot(guich_dat, aes(x = temp,y= prob, group = yolk, color = yolk)) +
  geom_point(shape = 19, size  = 8, position = pd) +
  geom_errorbar(aes(ymin  = prob+SE,
                    ymax  = prob-SE),
                width = 0.2,
                size  = 0.7,
                position = pd)+
  scale_color_manual(legend_title, values=c("black", "goldenrod"), labels=c('Yolk removal', 'Control'))+
  annotate(geom = "rect", xmin = 0, xmax = 1.5, 
           ymin = 0, ymax = .5, 
           fill = "#2C77BF", alpha = 0.2) +
  annotate(geom = "rect", xmin = 1.5, xmax = 3, 
           ymin = 0, ymax = .5,
           fill = "red", alpha = 0.2)+
  annotate("text", x = 1, y = .22, label = "16% Increase", fontface = "bold", size = 6)+
  annotate("text", x = 2.2, y = .22, label = "22% Increase", fontface = "bold",size = 6)+
  annotate("text", x = 1.45, y = .48, label = "Lampropholis guichenoti", fontface = "italic",size = 5)+
  theme_classic()+
  theme(legend.position = "none",
        text=element_text(size=14), 
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA))+
  xlab(bquote(Treatment~temperature~(degree*C)))+
  ylab("Probability of Mortality")+
  scale_y_continuous(breaks=seq(0, .5, .1), limits = c(0,.5), expand = expansion(mult = c(0, .005))) 

# bring in Guich image
Guich_image <- magick::image_read("Figures/Guich.png") %>% 
  magick::image_background("none")
image <- image_fill(Guich_image, 'none')
Guich_raster <- as.raster(image)

# Guich final plot
Guich_final <- ggdraw() +
  draw_plot(Guich_fig) +
  draw_image(Guich_raster, scale = .7, x = -.1, y= 0.33) 




###############################################
##############  Delicata  #####################
###############################################
######
# Delicata Survival Analysis
######
# quick look trt (temp and yolk)
ggstatsplot::ggbarstats(
  data = deli_df,
  x = liz_status,
  y = trt) +
  labs(caption = NULL) 

# Survival model
deli_model <- glmer(mortality ~temp + yolk +  temp:yolk+ (1|Egg_ID), data = deli_df, family = "binomial") 
summary(deli_model)
car::Anova(deli_model, type=3)
check_model(deli_model)
# df for figure
deli_dat <- emmeans(deli_model, specs = c("temp", "yolk"), type = "response") %>% 
  summary() %>%
  as.data.frame()
# percent diff
deli_23 <- round((0.33035714 - 0.09677419),digits = 2)
deli_28 <- round((0.26415094 - 0.10752688),digits = 2)


######
# Delienoti offspring mass on live subjects - no interaction; overall yolk effect
######
Deli_alive <- deli_df %>% filter(liz_status == "ALIVE") 
Deli_mass_model <- lm(log(hatch_mass)~ temp + yolk +  temp:yolk, data = Deli_alive)
# Checks and plots
summary(Deli_mass_model) 
anova(Deli_mass_model)
check_model(Deli_mass_model)
Deli_mass_emm <- emmeans(Deli_mass_model, specs = c("yolk"), type = "response") %>% # back transform from log scale
  summary() %>%
  as.data.frame()
plot(emmeans(Deli_mass_model, specs = c("yolk"), type = "response"))


######
# Delienoti offspring SVL on live subjects - no interaction; yolk effects
######
Deli_svl_model <- lm(log(hatch_svl)~ temp + yolk +  temp:yolk, data = Deli_alive)
# Checks and plots
summary(Deli_svl_model) 
anova(Deli_svl_model)
check_model(Deli_svl_model)
# extracting differences
Deli_svl_emm <- emmeans(Deli_svl_model, specs = c("yolk"), type = "response") %>% # back transform from log scale
  summary() %>%
  as.data.frame()
plot(emmeans(Deli_svl_model, specs = c("yolk"), type = "response"))


######
# Delienoti offspring BCI (OLS residual difference from svl and mass
######
Deli_bci <- Deli_alive %>% # filter out NA's to keep residuals from lm
  filter_at(vars(hatch_svl, hatch_mass), all_vars(!is.na(.)))
# calcualting residuals from OLS Regression of mass and SVL
Deli_fit <- lm(log(hatch_svl) ~ log(hatch_mass), data = Deli_bci)
check_model(Deli_fit)
Deli_bci$BCI <- residuals(Deli_fit)

# BCI test on yolk and temp - residuals are from ols of svl and mass
# yolk and temp effects; no interaction
Deli_bci_model <- lm(BCI~ temp + yolk +  temp:yolk, data = Deli_bci)
# Checks and plots
summary(Deli_bci_model) 
anova(Deli_bci_model)
check_model(Deli_bci_model)
# extracting differences
Deli_bci_emm <- emmeans(Deli_bci_model, specs = c("yolk", "temp"), type = "response") %>% # back transform from log scale
  summary() %>%
  as.data.frame()
plot(emmeans(Deli_bci_model, specs = c("yolk", "temp"), type = "response"))


########
# Delicata figures
########
legend_title <- "Resource Treatment"
pd = position_dodge(.8)

####### Survival on temp 
deli_fig <-ggplot(deli_dat, aes(x = temp,y= prob, group = yolk, color = yolk)) +
  geom_point(shape = 19, size  = 8, position = pd) +
  geom_errorbar(aes(ymin  = prob+SE,
                    ymax  = prob-SE),
                width = 0.2,
                size  = 0.7,
                position = pd)+
  scale_color_manual(legend_title, values=c("black", "goldenrod"), labels=c('Yolk removal', 'Control'))+
  annotate(geom = "rect", xmin = 0, xmax = 1.5, 
           ymin = 0, ymax = .5, 
           fill = "#2C77BF", alpha = 0.2) +
  annotate(geom = "rect", xmin = 1.5, xmax = 3, 
           ymin = 0, ymax = .5,
           fill = "red", alpha = 0.2)+
  annotate("text", x = 1, y = .2, label = "23% Increase", fontface = "bold", size = 6)+
  annotate("text", x = 2.2, y = .2, label = "16% Increase", fontface = "bold",size = 6)+
  annotate("text", x = 1.4, y = .48, label = "Lampropholis delicata", fontface = "italic",size = 5)+
  theme_classic()+
  theme(legend.justification=c(1,0), 
        legend.position=c(1,0),
        legend.background=element_rect(fill = alpha("white", 0.5)),
        text=element_text(size=14), 
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA))+
  xlab(bquote(Treatment~temperature~(degree*C)))+
  ylab("Probability of Mortality")+
  scale_y_continuous(breaks=seq(0, .5, .1), limits = c(0,.5), expand = expansion(mult = c(0, .005))) 

# bring in Deli image
Deli_image <- magick::image_read("Figures/Deli.png") %>% 
  magick::image_background("none")
image <- image_fill(Deli_image, 'none')
Deli_raster <- as.raster(image)

# deli final plot
Deli_final <- ggdraw() +
  draw_plot(deli_fig) +
  draw_image(Deli_raster, scale = .6, x = -.15, y= 0.33) 

# combind both plots
Final_survival <- plot_grid(Guich_final, Deli_final, labels = c('A', 'B'))






