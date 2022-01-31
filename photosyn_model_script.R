#script to run model for photosynthesis data -- similar to growth_rate script

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

run5_photosyn_data <- read.csv("/Users/Angela/src/work/limu/algal_growth_photosynthesis/data/run5_hyp_ulv_final_alpha_ek_rounded.csv")

get_run_number <- function(d) {
  # returns the run number based on the input date  
  run = NA
  if ("2021-09-20" <= d && d <= "2021-09-29") {
    run <- 1
  } else if ("2021-09-30" <= d && d <= "2021-10-09") {
    run <- 2
  } else if ("2021-10-11" <= d && d <= "2021-10-29") {
    run <- 3
  } else if ("2021-11-01" <= d && d <= "2021-11-12") {
    run <- 4
  } else {
    print("The date is %s not in a recognized range", d)
  }
  run
}

# make new column adding run to the dataset
run5_photosyn_data$run <- sapply(run5_photosyn_data$Date, get_run_number)
run5_photosyn_data$run <- as.factor(run5_photosyn_data$run)

#assigns temperature as a factor
run5_photosyn_data$temperature <- as.factor(run5_photosyn_data$Temperature)

#assigns treatment as characters from integers then to factors
run5_photosyn_data$treatment <- as.factor(as.character(run5_photosyn_data$Treatment))

#toggle between the species for output
hypnea <- subset(run5_photosyn_data, Species == "hm")
ulva <- subset(run5_photosyn_data, Species == "ul")

# run model with interaction between temperature and treatment for rETRmax (is best?)
run5_photosyn_model <- lmer(formula = rETRmax ~ treatment * temperature + (1 | run), data = hypnea)

#check the performance of the model for each dataset. If collinearity good, proceed with this model
#if not good, use model with no interaction between the fixed effect variables
performance::check_model(run5_photosyn_model, data = hypnea)


#ULVA 
#run model without interaction
run5_photosyn_model_noint <- lmer(formula = rETRmax ~ treatment + temperature + (1 | run), data = ulva)

#make a histogram and residual plots of the data for ulva
hist(ulva$rETRmax, main = paste("Ulva lactuca rETRmax"), col = "olivedrab3", labels = TRUE)
plot(resid(run5_photosyn_model_noint) ~ fitted(run5_photosyn_model_noint))
qqnorm(resid(run5_photosyn_model_noint))
qqline(resid(run5_photosyn_model_noint))

#check the performance of the model
performance::check_model(run5_photosyn_model_noint)

anova(run5_photosyn_model_noint, type = c("III"), ddf = "Satterthwaite")
ulva_photosyn_model_aov <- aov(rETRmax ~ treatment + temperature, data = ulva)
TukeyHSD(ulva_photosyn_model_aov, "treatment", ordered = FALSE)

r.squaredGLMM(run5_photosyn_model_noint)
summary(run5_photosyn_model_noint)








plot(allEffects(run5_photosyn_model_noint))

