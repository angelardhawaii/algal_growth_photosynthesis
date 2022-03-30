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

run5_6_photosyn_data <- read.csv("/Users/Angela/Library/Mobile Documents/com~apple~CloudDocs/research_limu/run5-6_hyp_ulv_temp/final_temperature_run1-4/final_data_files/run5_run6_all_photosyn.csv")

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
  } else if ("2022-02-11" <= d && d <= "2022-03-01") {
    run <- 0
  } else {
    print("The date is %s not in a recognized range", d)
  }
  run
}

# make new column adding run to the dataset
run5_6_photosyn_data$run <- sapply(run5_6_photosyn_data$Date, get_run_number)
run5_6_photosyn_data$run <- as.factor(run5_6_photosyn_data$run)

#assigns temperature as a factor
run5_6_photosyn_data$Temperature <- as.factor(run5_6_photosyn_data$Temperature)

#assigns treatment as characters from integers then to factors
run5_6_photosyn_data$Treatment <- as.factor(as.character(run5_6_photosyn_data$Treatment))

#toggle between the species for output. Use Day 9 for final analysis
hypnea <- subset(run5_6_photosyn_data, Species == "hm" & RLC.Day == 9)
ulva <- subset(run5_6_photosyn_data, Species == "ul" & RLC.Day == 9)

# run model with interaction between temperature and treatment for rETRmax (is best?)
#run5_6_photosyn_model <- lmer(formula = rETRmax ~ Treatment * Temperature + (1 | run), data = ulva)
#check the performance of the model for each dataset. If collinearity good, proceed with this model
#if not good, use model with no interaction between the fixed effect variables
#performance::check_model(run5_6_photosyn_model, data = ulva)


#ULVA 
#run model without interaction
run5_6_photosyn_model_noint <- lmer(formula = rETRmax ~ Treatment + Temperature + (1 | run), data = ulva)

#make a histogram and residual plots of the data for ulva
hist(ulva$rETRmax, main = paste("Ulva lactuca Ek"), col = "olivedrab3", labels = TRUE)
plot(resid(run5_6_photosyn_model_noint) ~ fitted(run5_6_photosyn_model_noint))
qqnorm(resid(run5_6_photosyn_model_noint))
qqline(resid(run5_6_photosyn_model_noint))

#check the performance of the model
performance::check_model(run5_6_photosyn_model_noint)

r.squaredGLMM(run5_6_photosyn_model_noint)
summary(run5_6_photosyn_model_noint)

#run ANOVA and pairwise comparisons
anova(run5_6_photosyn_model_noint, type = c("III"), ddf = "Satterthwaite")
ulva_photosyn_model_aov <- aov(rETRmax ~ Treatment + Temperature, data = ulva)
TukeyHSD(ulva_photosyn_model_aov, "Treatment", ordered = FALSE)
TukeyHSD(ulva_photosyn_model_aov, "Temperature", ordered = FALSE)


plot(allEffects(run5_6_photosyn_model_noint))


#HYPNEA 
#run model without interaction
run5_6_photosyn_model_noint <- lmer(formula = rETRmax ~ Treatment + Temperature + (1 | run), data = hypnea)

#make a histogram and residual plots of the data for hypnea
hist(hypnea$rETRmax, main = paste("Hypnea musciformis Ek"), col = "maroon", labels = TRUE)
plot(resid(run5_6_photosyn_model_noint) ~ fitted(run5_6_photosyn_model_noint))
qqnorm(resid(run5_6_photosyn_model_noint))
qqline(resid(run5_6_photosyn_model_noint))

#check the performance of the model
performance::check_model(run5_6_photosyn_model_noint)

r.squaredGLMM(run5_6_photosyn_model_noint)
summary(run5_6_photosyn_model_noint)

#run ANOVA and pairwise comparisons
anova(run5_6_photosyn_model_noint, type = c("III"), ddf = "Satterthwaite")
hypnea_photosyn_model_aov <- aov(rETRmax ~ Treatment + Temperature, data = hypnea)
TukeyHSD(hypnea_photosyn_model_aov, "Temperature", ordered = FALSE)



plot(allEffects(run5_6_photosyn_model_noint))





