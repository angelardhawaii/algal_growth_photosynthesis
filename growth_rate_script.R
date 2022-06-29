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
library(ggplot2)

#open weight dataset and make columns for growth rate from initial and final weights
run5.6_growth <- read.csv("/Users/Angela/src/work/limu/algal_growth_photosynthesis/data_input/run5-6_growth_all_042922.csv")

#make a new column for weight change (difference final from initial)
run5.6_growth$growth_rate_percent <- (run5.6_growth$final.weight - run5.6_growth$Initial.weight) / run5.6_growth$Initial.weight * 100

#make a new column for daily growth rate from 8 day study (steady growth rate assumed rather than exponential)
run5.6_growth$steady_growth_daily <- run5.6_growth$growth_rate_percent / 8

#make a new column that keeps only the numerical values (effectively removes the C in temperatures for consistency)
run5.6_growth$temp_clean <- as.factor(substr(run5.6_growth$temperature, 1, 2))

#assigns temperature as a factor
run5.6_growth$temperature <- as.factor(run5.6_growth$temp_clean)

#assigns treatment as characters from integers then to factors
run5.6_growth$treatment <- as.factor(as.character(run5.6_growth$treatment))

#assign run as a factor
run5.6_growth$run <- as.factor(run5.6_growth$run)

#assign plant ID as a factor
run5.6_growth$plant.ID <- as.factor(run5.6_growth$plant.ID)

#assign RLC order as a factor
run5.6_growth$RLC.order <- as.factor(run5.6_growth$RLC.order)

# assign lunar phase as factor
run5.6_growth$lunar.phase <- as.factor(run5.6_growth$lunar.phase)

#toggle between the species for output
hypnea <- subset(run5.6_growth, Species == "Hm")
ulva <- subset(run5.6_growth, Species == "Ul")

# run model with interaction between temperature and treatment 
#run5.6_growth_model <- lmer(formula = growth_rate_percent ~ treatments * temp_clean + (1 | run), data = hypnea)

#ULVA 
#run model without interaction (0 in model permits display of four levels of treatments - no intercept)
run5.6_growth_model_noint <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | run) + (1 | plant.ID) + (1 | RLC.order), data = ulva)
#make a histogram of the data for ulva
hist(ulva$growth_rate_percent, main = paste("Ulva lactuca Growth Rate (%)"), col = "olivedrab3", labels = TRUE)
ulva_growth_model_aov <- aov(growth_rate_percent ~ treatment + temperature, data = ulva)
TukeyHSD(ulva_growth_model_aov, "treatment", ordered = FALSE)

#HYPNEA 
#run model without interaction  (0 in model permits display of four levels of treatments - no intercept)
run5.6_growth_model_noint <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | run) + (1 | plant.ID) + (1 | RLC.order), data = hypnea)
#OR make a histogram of the data for hypnea
hist(hypnea$growth_rate_percent, main = paste("Hypnea musciformis Growth Rate (%)"), col = "maroon", labels = TRUE)
hypnea_growth_model_aov <- aov(growth_rate_percent ~ treatment + temperature, data = hypnea)
TukeyHSD(hypnea_growth_model_aov, "treatment", ordered = FALSE)

#check the performance of the model for each dataset: ulva and hypnea
performance::check_model(run5.6_growth_model_noint)

r.squaredGLMM(run5.6_growth_model_noint)
summary(run5.6_growth_model_noint)
ranef(run5.6_growth_model_noint)



plot(allEffects(run5.6_growth_model_noint))
anova(run5.6_growth_model_noint, type = c("III"), ddf = "Satterthwaite")
plot(resid(run5.6_growth_model_noint) ~ fitted(run5.6_growth_model_noint))
qqnorm(resid(run5.6_growth_model_noint))
qqline(resid(run5.6_growth_model_noint))
