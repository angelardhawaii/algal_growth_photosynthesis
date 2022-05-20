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
library(ggplot2)
library(ggpubr)

run5_6_photosyn_data <- read.csv("data_input/run5-6_ek_alpha.csv")

# assign run as a factor
run5_6_photosyn_data$Run <- as.factor(run5_6_photosyn_data$Run)

#assign temperature as a factor
run5_6_photosyn_data$Temperature <- as.factor(run5_6_photosyn_data$Temp...C.)

#assigns treatment as characters from integers then to factors
run5_6_photosyn_data$Treatment <- as.factor(as.character(run5_6_photosyn_data$Treatment))

#toggle between the species for output. Use Day 9 for final analysis
hypnea <- subset(run5_6_photosyn_data, Species == "hm" & RLC.Day == 9)
ulva <- subset(run5_6_photosyn_data, Species == "ul" & RLC.Day == 9)

#add growth rate from other dataset to this one and subset by species
growth_rate <- read.csv("/Users/Angela/src/work/limu/algal_growth_photosynthesis/data_input/run5-6_growth_all_042922.csv")
gr_ulva <- subset(growth_rate, Species == "Ul")
gr_hypnea <- subset(growth_rate, Species == "Hm" & final.weight != 0.1017)

ulva$growth_rate <- round((gr_ulva$final.weight - gr_ulva$Initial.weight) / gr_ulva$Initial.weight * 100, digits = 2)
hypnea$growth_rate <- round((gr_hypnea$final.weight - gr_hypnea$Initial.weight) / gr_hypnea$Initial.weight * 100, digits = 2)

# run model with interaction between temperature and treatment for rETRmax (is best?)
#run5_6_photosyn_model <- lmer(formula = rETRmax ~ Treatment * Temperature + (1 | run), data = ulva)
#check the performance of the model for each dataset. If collinearity good, proceed with this model
#if not good, use model with no interaction between the fixed effect variables
#performance::check_model(run5_6_photosyn_model, data = ulva)


#ULVA 
#run model without interaction
run5_6_photosyn_model_noint <- lmer(formula = rETRmax ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva)

#make a histogram and residual plots of the data for ulva
hist(ulva$rETRmax, main = paste("Ulva lactuca"), col = "olivedrab3", labels = TRUE)
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

#plot a regression between the photosynthetic independent variables of interest and growth rate
ulva_growth_etr_graph <- ggplot(ulva, aes(x=rETRmax, y=growth_rate)) + geom_point() + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Ulva lactuca rETRmax vs Growth Rate", x = "rETRmax (μmols electrons m-2 s-1)", 
       y = "growth rate (%)") + stat_regline_equation(label.x = 25, label.y = 165) + stat_cor()
ulva_growth_etr_graph

ulva_growth_ek_graph <- ggplot(ulva, aes(x=ek.1, y=growth_rate)) + geom_point() + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Ulva lactuca Ek vs Growth Rate", x = "Ek (μmols photons m-2 s-1)", y = "growth rate (%)") + 
  stat_regline_equation() + stat_cor(label.y = 160)
ulva_growth_ek_graph

#HYPNEA 
#run model without interaction
run5_6_photosyn_model_noint <- lmer(formula = rETRmax ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea)
run5_6_photosyn_model_noint <- lmer(formula = ek.1 ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea)

#make a histogram and residual plots of the data for hypnea
hist(hypnea$rETRmax, main = paste("Hypnea musciformis"), col = "maroon", labels = TRUE)
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

#plot a regression between the photosynthetic independent variables of interest and growth rate
hypnea_growth_etr_graph <- ggplot(hypnea, aes(x=rETRmax, y=growth_rate)) + geom_point() + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Hypnea musciformis rETRmax vs Growth Rate", x = "rETRmax (μmols electrons m-2 s-1)", 
       y = "growth rate (%)") + stat_regline_equation(label.x = 25, label.y = 165) + stat_cor(label.y = 150)
hypnea_growth_etr_graph

hypnea_growth_ek_graph <- ggplot(hypnea, aes(x=ek.1, y=growth_rate)) + geom_point() + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Hypnea musciformis Ek vs Growth Rate", x = "Ek (μmols photons m-2 s-1)", y = "growth rate (%)") + 
  stat_regline_equation(label.y = 155) + stat_cor(label.y = 145)
hypnea_growth_ek_graph

# do not need this any longer for this dataset
#get_run_number <- function(d) {
# returns the run number based on the input date  
# run = NA
#if ("2021-09-20" <= d && d <= "2021-09-29") {
# run <- 1
#} else if ("2021-09-30" <= d && d <= "2021-10-09") {
#  run <- 2
#} else if ("2021-10-11" <= d && d <= "2021-10-29") {
# run <- 3
#} else if ("2021-11-01" <= d && d <= "2021-11-12") {
# run <- 4
#} else if ("2022-02-11" <= d && d <= "2022-03-01") {
# run <- 0
#} else {
# print("The date is %s not in a recognized range", d)
#}
#run
#}
#run5_6_photosyn_data$run <- sapply(run5_6_photosyn_data$Date, get_run_number)


