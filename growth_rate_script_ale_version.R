## Hypnea and Ulva Photosynthesis and Growth -- 2021
# Angela Richards Donà
# January 2022

# open packages for mixed model effects analysis
library(lme4)
library(lmerTest)
library(effects)
library(car)
library(MuMIn)
library (dplyr)
library(emmeans)
library(DHARMa)
library(performance)
library(patchwork)

#open weight dataset and make columns for growth rate from initial and final weights
run5_growth <- read.csv("/Users/Angela/Library/Mobile Documents/com~apple~CloudDocs/research_limu/run5_hyp_ulv_temp/final_temperature_run1-4/final_data_files/run5_growth_all_012422.csv")

#make a new column for weight change (difference final from initial)
run5_growth$growth_rate_percent <- (run5_growth$final.weight - run5_growth$Inital.weight) / run5_growth$Inital.weight * 100

#make a new column for daily growth rate from 8 day study (steady growth rate assumed rather than exponential)
run5_growth$steady_growth_daily <- run5_growth$growth_rate_percent / 8


#Change existing columns and assign temperature as a factor
run5_growth$temperature <- as.factor(run5_growth$temperature)

#make a new column that changes treatment as integers to characters (note similar column name)
run5_growth$treatment <- as.factor(as.character(run5_growth$treatment))

#make new column that changes run to character
run5_growth$run <- as.factor(run5_growth$run)

run_model <- function(the_data, histogram_title, histogram_color) {
  run5_growth_model_noint <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | run), data = the_data)
  #OR make a histogram of the data for the_data
  hist(the_data$growth_rate_percent, main = paste(histogram_title), col = histogram_color, labels = TRUE)
  anova(run5_growth_model_noint, type = c("III"), ddf = "Satterthwaite")
  
  growth_model_aov <- aov(growth_rate_percent ~ treatment + temperature, the_data)
  TukeyHSD(growth_model_aov, "treatment", ordered = FALSE)

  r.squaredGLMM(run5_growth_model_noint)
  summary(run5_growth_model_noint)
  
  run5_growth_model_noint
}

do_plots <- function(model) {
  plot(allEffects(model))
  #check the performance of the model 
  plot(resid(model) ~ fitted(model))
  qqnorm(resid(model))
  qqline(resid(model))
}

the_data <- subset(run5_growth, Species == "Hm")
model <- run_model(the_data, "Hypnea musciformis Growth Rate (%)", "maroon")
do_plots(model)
performance::check_model(model)

the_data <- subset(run5_growth, Species == "Ul")
model <- run_model(the_data, "Ulva lactuca Growth Rate (%)", "olivedrab3")
do_plots(model)
performance::check_model(model)
