pacman::p_load("dplyr","lme4", "foreign", "ggplot2", "corrplot", "gridExtra", "sjPlot", "tidyverse", "lubridate","visreg", "emmeans", "performance", "cowplot", "patchwork", "png", "magick")

# importdata (be sure to change Y/N to 0/1)
data <- read.csv("Data/Lizard_Db.csv", head = T)

#Removing NAs from dataset which is used in analysis
data_raw <- data[!is.na(data$liz_status),]
table(data_raw$liz_status, data_raw$trt) # drop missing and unkowns
data_raw <- data_raw %>% 
  filter(liz_status != "MISSING" & liz_status != "UNKNOWN") 
data_raw$temp = as.factor(data_raw$temp)

# cleaning up dates from lizard db
data_final <-data_raw %>%
  # Fix the messed up date columns
  # Here I'm separating all of the hatch dates into their components of year, month, and day.
  separate(col = 'hatch_date_.year.m.day.', into = c('y', 'm', 'd'), sep = "/") %>%
  # Here I'm using position number two as a separator to remove the '20' from the part of the date excel interpreted as the year.
  separate(col = d, into = c('millenia', 'd'), sep = 2) %>%
  # Removing that pointless '20' column.
  select(-(millenia)) %>%
  # Using the unite function to join the separated year, month, and day columns back into a single column
  unite(col = hatch_date, c('y','m','d'), sep = "/") %>%
  # R turned missing values into 'NA/NA/NA' because of the way I used unite, so this is just changing those back into missing values.
  mutate(hatch_date = ifelse(hatch_date != 'NA/NA/NA', hatch_date, NA)) %>%
  separate(col = 'date_egg_laid_.year.m.day.', into = c('y', 'm', 'd'), sep = "/") %>%
  # Repeating the same process with the egg date
  separate(col = d, into = c('millenia', 'd'), sep = 2) %>%
  select(-(millenia)) %>%
  unite(col = egg_date, c('y','m','d'), sep = "/") %>%
  mutate(egg_date = ymd(egg_date),
         death_date = dmy(death_date)) %>%
  # Using the lubridate package, I'm turning the hatch_date column I just made into a datetime object with structure 'Year/Month/Day' with ymd()
  mutate(hatch_date = ymd(hatch_date),
         dev_time = ifelse(is.na(hatch_date), NA, difftime(hatch_date, egg_date, units = "days")),
         days_alive = ifelse(is.na(hatch_date), NA, difftime(death_date, hatch_date, units = "days"))) %>% 
  mutate(hatch_mass = as.numeric(hatch_mass))

# separate the two species into df
guich_df <-data_final %>% filter(species == "guichenoti")
deli_df <- data_final %>% filter(species == "delicata")


###############################################
############## Guichenoti #####################
###############################################
######
# Guichenoti Survival Analysis
######
# Quick look trt (temp and yolk)
ggstatsplot::ggbarstats(
  data = deli_df,
  x = liz_status,
  y = trt) +
  labs(caption = NULL) 

# Guich survival model
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





