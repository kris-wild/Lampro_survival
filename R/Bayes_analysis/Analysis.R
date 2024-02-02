pacman::p_load("dplyr","lme4", "MASS", "brms", "MCMCglmm", "quantreg","lmerTest", "latex2exp", "tidybayes", "bayesplot", "rstanarm", "plotrix", "emmeans", "ggExtra", "marginaleffects", "gridExtra", "bayestestR","performance", "cowplot", "patchwork", "png", "magick", "ggplot2")

# import data (be sure to change Y/N to 0/1)
data <- read.csv("Data/Lizard_Db_final_survival_est.csv", 
                 head = T) %>% 
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


##############################
####### Extracting posteriors
##############################
post_Deli_m1 <- posterior_samples(Deli_m1_brms, pars = "^b")
variable.names(post_Deli_m1)

### yolk ablation posteriors
yolk_a_23 <- as.array(post_Deli_m1[,"b_Intercept"]) # Ablation at 23C
yolk_a_28 <- as.array(post_Deli_m1[,"b_Intercept"] + (post_Deli_m1[,"b_temp28"])) #Ablation at 28C
##probs
prob_a_23 <- exp(yolk_a_23)/(1+exp(yolk_a_23))
# prob_a_23
prob_a_23 <- prob_a_23 %>% 
  as.data.frame() %>% 
  mutate(temp = as.character("23"),
         yolk = "A")
# prob_a_28
prob_a_28 <- exp(yolk_a_28)/(1+exp(yolk_a_28))
prob_a_28 <- prob_a_28 %>% 
  as.data.frame() %>% 
  mutate(temp = as.character("28"),
         yolk = "A")

### yolk controls posteriors
yolk_c_23 <- as.array(yolk_a_23 + (post_Deli_m1[,"b_yolkC"])) # Control at 23C
yolk_c_28 <- as.array(post_Deli_m1[,"b_Intercept"]) + (post_Deli_m1[,"b_yolkC:temp28"]) # at 28C
## probs
#prob_c_23
prob_c_23 <- exp(yolk_c_23)/(1+exp(yolk_c_23))
prob_c_23 <- prob_c_23 %>% 
  as.data.frame() %>% 
  mutate(temp = as.character("23"),
         yolk = "C")
# prob_c_28
prob_c_28 <- exp(yolk_c_28)/(1+exp(yolk_c_28)) # check this
prob_c_28 <- prob_c_28 %>% 
  as.data.frame() %>% 
  mutate(temp = as.character("28"),
         yolk = "C")

# PRIOR DF 
Deli_prior_dat <- rbind(prob_a_23, prob_a_28, prob_c_23, prob_c_28) colnames(Deli_prior_dat)[1] = "Mortality" 


### arranging summary data for figure
# raw data sample size for figure
deli_df_sum <- deli_df %>% 
  group_by(trt) %>% 
  summarise(mean = mean(egg_mass_final, na.rm=TRUE),
            n = n())

# BRMS summary for figure 
Deli_effects <- conditional_effects(Deli_m1_brms, "temp:yolk")
Deli_effects <- Deli_effects$`temp:yolk`
Deli_effects$trt <- paste(Deli_effects$temp, Deli_effects$yolk)

# percent diff (yolk - control) for figure
Deli_23 <- round((0.28760162 - 0.07938782),digits = 3)
Deli_28 <- round((0.22413991 - 0.08588692),digits = 3)

# egg mass mortality plot for inset
Deli_c_eff <- conditional_effects(Deli_m1_brms, "egg_mass_final") 
Deli_inset_plot <- plot(Deli_c_eff, plot = FALSE)[[1]] + 
  labs(x = "Egg mass (g)", y = "Probability of mortality")+
  theme_bw()+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  ) 


##############################
####### Deli figure 
##############################
legend_title <- "Incubation"
pd = position_dodge(.8)
deli_labs <- c("n = 112\n23 Yolk Removal", " n = 93\n23 Control", 
               "n = 106\n28 Yolk Removal", "n = 93\n28 Control")
### GGPLOT
Deli_fig <-ggplot(Deli_effects, aes(x = trt, y= estimate__, group = temp, color = temp)) +
  geom_errorbar(aes(ymin  = lower__,
                    ymax  = upper__),
                width = 0.2,
                size  = 0.7,
                position = pd,
                color = "black")+
  geom_point(shape = 19, size  = 9, position = pd) +
  scale_color_manual(values=c("blue", "red"), labels=c('Cold', 'Hot'))+
  annotate(geom = "rect", xmin = 0, xmax = 2.5, 
           ymin = 0, ymax = .5, 
           fill = "#2C77BF", alpha = 0.2) +
  annotate(geom = "rect", xmin = 2.5, xmax = 5, 
           ymin = 0, ymax = .5,
           fill = "red", alpha = 0.2)+
  annotate("text", x = 1.6, y = .185, label = "21% Increase", fontface = "bold", size = 5)+
  annotate("text", x = 3.75, y = .185, label = "14% Increase", fontface = "bold",size = 5)+
  annotate("text", x = 1.1, y = .48, label = "Lampropholis delicata", fontface = "italic",size = 5)+
  scale_x_discrete(labels= deli_labs)+
  theme_classic()+
  theme(legend.position= "none", 
        text=element_text(size=14), 
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA))+
  xlab(bquote(Treatment~temperature~(degree*C)))+
  ylab("Probability of Mortality")+
  scale_y_continuous(breaks=seq(0, .5, .1), limits = c(0,.5), expand = expansion(mult = c(0, .005)))

# add inset elememnt to ggplot fig
Deli_gg <- Deli_fig + inset_element(Deli_inset_plot, 0.65, 0.65, 1, 1)

# bring in Deli image
Deli_image <- magick::image_read("Figures/Deli.png") %>% 
  magick::image_background("none")
image <- image_fill(Deli_image, 'none')
Deli_raster <- as.raster(image)

# deli final plot
Deli_final <- ggdraw() +
  draw_plot(Deli_gg) +
  draw_image(Deli_raster, scale = .6, x = -.2, y= 0.35) 







###########
# Functions
###########
## pMCMC Function
# Calculates the, p-value or pMCMC value for a posterior distribution
# x The vector for the posterior distribution. Note that this will test the null hypothesis that the parameter of interest is significantly different from 0. 
pmcmc <- function(x){
  2*(1 - max(table(x<0) / nrow(x)))
}


deli_df$temp_cont <- as.numeric(deli_df$temp)



####################################
# Deli Survival Analysis Numeric Temp and egg mass interaction for 
# adult model 
######
# Deli survival model
Deli_m2_brms <- brm(mortality ~ egg_mass_final + temp_cont + # egg mass and temp continious
                      temp_cont:egg_mass_final   + # interaction
                      (1|clutch),  # clutch random factor
                    family = bernoulli(link = "logit"), 
                    data = deli_df, 
                    iter= 5000, warmup = 1000, 
                    thin = 4, cores = 8)


####################
# Model 2 checks: lags, residuals, r2, summary
####################
# checking lags for this model 
variable.names(post_Deli_m2)
draws <- as.array(Deli_m2_brms)
mcmc_acf(draws,  pars = c("b_Intercept","b_egg_mass_final", "b_temp_cont", "b_egg_mass_final:temp_cont"), lags =10)
# plots
plot(Deli_m2_brms)
#R2 and summary of full model
bayes_R2(Deli_m2_brms)
summary(Deli_m2_brms)

# stan plot
stanplot(Deli_m2_brms, 
         type = "areas",
         prob = 0.95)

# overall egg mass on survival
plot(conditional_effects(Deli_m2_brms, "egg_mass_final"), 
     ask = FALSE)

# incubation temperature X mass treatment interaction
plot(conditional_effects(Deli_m2_brms, 
                         "egg_mass_final:temp_cont", 
                         points = TRUE))





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

# egg mass by yolk treatment on survival
plot(conditional_effects(Guich_m1_brms, "egg_mass_final:yolk"), ask = FALSE) 

# temperature X yolk treatment interaction
plot(conditional_effects(Guich_m1_brms, "yolk:temp", points = TRUE))


##############################
####### Extracting posteriors
##############################
post_Guich_m1 <- posterior_samples(Guich_m1_brms, pars = "^b")
variable.names(post_Guich_m1)

### yolk ablation posteriors
yolk_a_23 <- as.array(post_Guich_m1[,"b_Intercept"]) # Ablation at 23C
yolk_a_28 <- as.array(post_Guich_m1[,"b_Intercept"] + (post_Guich_m1[,"b_temp28"])) #Ablation at 28C
##probs
prob_a_23 <- exp(yolk_a_23)/(1+exp(yolk_a_23))
# prob_a_23
prob_a_23 <- prob_a_23 %>% 
  as.data.frame() %>% 
  mutate(temp = as.character("23"),
         yolk = "A")
# prob_a_28
prob_a_28 <- exp(yolk_a_28)/(1+exp(yolk_a_28))
prob_a_28 <- prob_a_28 %>% 
  as.data.frame() %>% 
  mutate(temp = as.character("28"),
         yolk = "A")

### yolk controls posteriors
yolk_c_23 <- as.array(yolk_a_23 + (post_Guich_m1[,"b_yolkC"])) # Control at 23C
yolk_c_28 <- as.array(post_Guich_m1[,"b_Intercept"]) + (post_Guich_m1[,"b_yolkC:temp28"]) # at 28C
## probs
#prob_c_23
prob_c_23 <- exp(yolk_c_23)/(1+exp(yolk_c_23))
prob_c_23 <- prob_c_23 %>% 
  as.data.frame() %>% 
  mutate(temp = as.character("23"),
         yolk = "C")
# prob_c_28
prob_c_28 <- exp(yolk_c_28)/(1+exp(yolk_c_28)) # check this
prob_c_28 <- prob_c_28 %>% 
  as.data.frame() %>% 
  mutate(temp = as.character("28"),
         yolk = "C")
# final naming 
guich_prior_dat <- rbind(prob_a_23, prob_a_28, prob_c_23, prob_c_28) 
colnames(guich_prior_dat)[1] = "Mortality" 

# summary DF of Priors
guich_95CI <- guich_prior_dat %>% 
  group_by(yolk, temp) %>% 
  summarise(mean = mean(Mortality, na.rm = TRUE),
           sd = sd(Mortality, na.rm = TRUE),
           n = n()) %>%
  mutate(se = sd / sqrt(n))

# raw data sample size for figure
guich_df_sum <- guich_df %>% 
  group_by(trt) %>% 
  summarise(mean = mean(egg_mass_final, na.rm=TRUE),
            n = n())

# BRMS summary for figure 
Guich_effects <- conditional_effects(Guich_m1_brms, "temp:yolk")
Guich_effects <- Guich_effects$`temp:yolk`
Guich_effects$trt <- paste(Guich_effects$temp, Guich_effects$yolk)


# percent diff (yolk - control) for figure
Guich_23 <- round((0.27885490 - 0.12658309),digits = 3)
Guich_28 <- round((0.31802495 - 0.09085381),digits = 3)

# egg mass mortality plot for inset
Guich_c_eff <- conditional_effects(Guich_m1_brms, "egg_mass_final") 
Guich_inset_plot <- plot(Guich_c_eff, plot = FALSE)[[1]] + 
  labs(x = "Egg mass (g)", y = "Probability of mortality")+
  theme_bw()+
  theme(
    panel.background = element_rect(fill='transparent'),
    plot.background = element_rect(fill='transparent', color=NA),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.background = element_rect(fill='transparent'),
    legend.box.background = element_rect(fill='transparent')
  ) 
  
##############################
####### Guich figure 
##############################
legend_title <- "Incubation"
pd = position_dodge(.8)
labs <- c("n = 75\n23 Yolk Removal", " n = 53\n23 Control", 
          "n = 73\n28 Yolk Removal", "n = 50\n28 Control")
### GGPLOT
Guich_fig <-ggplot(Guich_effects, aes(x = trt, y= estimate__, group = temp, color = temp)) +
  geom_errorbar(aes(ymin  = lower__,
                    ymax  = upper__),
                width = 0.2,
                size  = 0.7,
                position = pd,
                color = "black")+
  geom_point(shape = 19, size  = 9, position = pd) +
  scale_color_manual(values=c("blue", "red"), labels=c('Cold', 'Hot'))+
  annotate(geom = "rect", xmin = 0, xmax = 2.5, 
           ymin = 0, ymax = .5, 
           fill = "#2C77BF", alpha = 0.2) +
  annotate(geom = "rect", xmin = 2.5, xmax = 5, 
           ymin = 0, ymax = .5,
           fill = "red", alpha = 0.2)+
  annotate("text", x = 1.8, y = .28, label = "15% Increase", fontface = "bold", size = 5)+
  annotate("text", x = 4, y = .22, label = "23% Increase", fontface = "bold",size = 5)+
  annotate("text", x = 1.1, y = .48, label = "Lampropholis guichenoti", fontface = "italic",size = 5)+
  scale_x_discrete(labels= labs)+
  theme_classic()+
  theme(legend.justification=c(0,1), 
        legend.position=c(0,0.12),
        legend.title=element_blank(),
        legend.background=element_rect(fill = alpha("white", 0.5)),
        text=element_text(size=14), 
        panel.background = element_rect(fill='transparent'), 
        plot.background = element_rect(fill='transparent', color=NA))+
  xlab(bquote(Treatment~temperature~(degree*C)))+
  ylab("Probability of Mortality")+
  scale_y_continuous(breaks=seq(0, .5, .1), limits = c(0,.5), expand = expansion(mult = c(0, .005))) 

# add inset elememnt to ggplot fig
Guich_gg <- Guich_fig + inset_element(Guich_inset_plot, 0.65, 0.65, 1, 1)

# bring in Deli image
Guich_image <- magick::image_read("Figures/Guich.png") %>% 
  magick::image_background("none")
image <- image_fill(Guich_image, 'none')
Guich_raster <- as.raster(image)

# Guich final plot
Guich_final <- ggdraw() +
  draw_plot(Guich_gg) +
  draw_image(Guich_raster, scale = .6, x = -.2, y= 0.38) 



# combind both plots
Final_survival <- plot_grid(Guich_final, Deli_final, labels = c('A', 'B'))























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








