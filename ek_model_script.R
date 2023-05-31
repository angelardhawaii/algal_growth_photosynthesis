#script to run model for Ek -- similar to rETRmax and growth_rate scripts

library(lme4)
library(lmerTest)
library(afex)
library(effects)
library(car)
library(MuMIn)
library (dplyr)
library(emmeans)
library(DHARMa)
library(performance)
library(patchwork)
library(rstatix)
#for plots and tables
library(ggplot2)
library(ggpubr)
library(forcats)
library(RColorBrewer)
library(tidyverse)
library(sjPlot)
library(sjmisc)
library(mmtable2)
library(gt)
library(purrr)
library(stringr)
library(tidyr)

#load this file for rETR-based use per Silsbe and Kromkamp 2012
#all_runs_photosyn_data <- read.csv("data_input/hyp_ulva_all_runs_ek_alpha.csv")

#load this file for normalized to quantum efficiency of photosynthesis per same as above
all_runs_photosyn_data <- read.csv("data_input/hyp_ulva_all_runs_ek_alpha_normalized.csv")

# assign run as a factor
all_runs_photosyn_data$Run <- as.factor(all_runs_photosyn_data$Run)

#assign temperature as a factor
all_runs_photosyn_data$Temperature <- as.factor(all_runs_photosyn_data$Temp...C.)

#assigns treatment as characters from integers then to factors
all_runs_photosyn_data$Treatment <- as.factor(as.character(all_runs_photosyn_data$Treatment))

# assign deltaNPQ as a factor
all_runs_photosyn_data$deltaNPQ <- as.factor(all_runs_photosyn_data$deltaNPQ)


#subset data to get only Ulva data for day 9 Ek-- also remove the treatment 2.5 that was problematic

ulva <- subset(all_runs_photosyn_data, Species == "ul" & RLC.Day == 9 & Treatment != 2.5)
ulva$treatment_graph[ulva$Treatment == 0] <- "1) 35ppt/0.5umol"
ulva$treatment_graph[ulva$Treatment == 1] <- "2) 35ppt/14umol" 
ulva$treatment_graph[ulva$Treatment == 2] <- "3) 28ppt/27umol" 
ulva$treatment_graph[ulva$Treatment == 3] <- "5) 18ppt/53umol" 
ulva$treatment_graph[ulva$Treatment == 4] <- "6) 11ppt/80umol"
ulva$treatment_graph[ulva$Treatment == 2.5] <- "4) 28ppt/53umol"

#Ulva Ek________________________________________________________________________________________
ulva %>% ggplot(aes(ek.1)) +
  geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()
hist(ulva$ek.1, main = paste("Ulva lactuca"), col = "olivedrab3", labels = TRUE)
#Ulva Ek model

ulva_ek_model <- lmer(formula = ek.1 ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva)
summary(ulva_ek_model)

#construct null model to perform likelihood ratio test REML must be FALSE
ulva_ek_treatment_null <- lmer(formula = ek.1 ~ Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
ulva_ek_model2 <- lmer(formula = ek.1 ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
anova(ulva_ek_treatment_null, ulva_ek_model2)
ulva_ek_temperature_null <- lmer(formula = ek.1 ~ Treatment + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
ulva_ek_model3 <- lmer(formula = ek.1 ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
anova(ulva_ek_temperature_null, ulva_ek_model3)

#plot and tables
tab_model(ulva_ek_model)
plot(allEffects(ulva_ek_model))
#make a histogram and residual plots of the data for ulva ek
hist(resid(ulva_ek_model))
plot(resid(ulva_ek_model) ~ fitted(ulva_ek_model))
qqnorm(resid(ulva_ek_model))
qqline(resid(ulva_ek_model))

#check the performance of the model
performance::check_model(ulva_ek_model)
r.squaredGLMM(ulva_ek_model)

#plot raw data
ulva %>% ggplot(aes(treatment_graph, ek.1)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 3, aes(color = Temperature), position = "jitter", show.legend = FALSE) + 
  labs(x="salinity/nitrate", y= "Day 9 Ek (μmols photons m-2 s-1)", title= "C", subtitle = "Ulva lactuca") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  ylim(-1, 250) + stat_mean() +
  scale_color_manual(values = c("#295102", "#7CB950", "#BDE269")) +
  geom_hline(yintercept=50, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#Make linear regression with growth
#add growth rate from other dataset to this one and subset by species
growth_rate <- read.csv("/Users/Angela/src/work/limu/algal_growth_photosynthesis/data_input/all_runs_growth_011723.csv")
growth_rate$Species <- as.factor(growth_rate$Species)
growth_rate$treatment <- as.factor(growth_rate$treatment)
growth_rate$lunar.phase <- as.factor(growth_rate$lunar.phase)

#make a new column for weight change (difference final from initial)
growth_rate$growth_rate_percent <- (growth_rate$final.weight - growth_rate$Initial.weight) / growth_rate$Initial.weight * 100

gr_ulva <- subset(growth_rate, Species == "Ul" & treatment != 2.5)

ulva$growth_rate <- round((gr_ulva$final.weight - gr_ulva$Initial.weight) / gr_ulva$Initial.weight * 100, digits = 2)
ulva$lunar.phase <- (gr_ulva$lunar.phase)


#plot linear regression between ek and growth
ulva_growth_ek_graph <- ggplot(ulva, aes(x=ek.1, y=growth_rate)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = Treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Ulva lactuca Ek vs Growth Rate", x = "Ek (μmols photons m-2 s-1)", 
       y = "growth rate (%)") + stat_regline_equation(label.x = 25, label.y = 165) + stat_cor()
ulva_growth_ek_graph

#summarize the means for Ek
ulva %>% group_by(Treatment) %>% summarise_at(vars(ek.1), list(mean = mean))
ulva %>% group_by(Temperature) %>% summarise_at(vars(ek.1), list(mean = mean))


#Hypnea Ek____________________________________________________________________
#toggle between the species for output. Use Day 9 for final analysis
#There was no D9 RLC for hm6-4 on 11/12/21 but had to remove hm6-4 from 10/9/21 below to match growth data
hypnea <- subset(all_runs_photosyn_data, Species == "hm" & RLC.Day == 9 & uid != "2021-10-09_hm6-4")
hypnea$treatment_graph[hypnea$Treatment == 0] <- "1) 35ppt/0.5umol"
hypnea$treatment_graph[hypnea$Treatment == 1] <- "2) 35ppt/14umol" 
hypnea$treatment_graph[hypnea$Treatment == 2] <- "3) 28ppt/27umol" 
hypnea$treatment_graph[hypnea$Treatment == 3] <- "5) 18ppt/53umol" 
hypnea$treatment_graph[hypnea$Treatment == 4] <- "6) 11ppt/80umol"
hypnea$treatment_graph[hypnea$Treatment == 2.5] <- "4) 28ppt/53umol"

#hm1-1 on 10/12/22 and hm1-2 on 4/29/22 causing issues of influential observations
#hm1-2 for both rETRmax and Ek -- leaving them in dataset because no good reason to believe not good data

#make a histogram and residual plots of the data for hypnea
hist(hypnea$ek.1, main = paste("Hypnea musciformis -- Ek"), col = "maroon2", labels = TRUE)
hypnea %>% ggplot(aes(ek.1)) +
  geom_histogram(binwidth=5, fill = "#7D0033", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

###run model for Hypnea and Ek
hyp_ek_model <- lmer(formula = ek.1 ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID), data = hypnea)
summary(hyp_ek_model)

#construct null model to perform likelihood ratio test REML must be FALSE
hyp_ek_treatment_null <- lmer(formula = ek.1 ~ Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea, REML = FALSE)
hyp_ek_model2 <- lmer(formula = ek.1 ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea, REML = FALSE)
anova(hyp_ek_treatment_null, hyp_ek_model2)

hyp_ek_temp_null <- lmer(formula = ek.1 ~ Treatment + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea, REML = FALSE)
hyp_ek_model3 <- lmer(formula = ek.1 ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea, REML = FALSE)
anova(hyp_ek_temp_null, hyp_ek_model3)

hist(resid(hyp_ek_model))
plot(resid(hyp_ek_model) ~ fitted(hyp_ek_model))
qqnorm(resid(hyp_ek_model))
qqline(resid(hyp_ek_model))

#check the performance of the model
performance::check_model(hyp_ek_model)
r.squaredGLMM(hyp_ek_model)


#make a table and plot for model
tab_model(hyp_ek_model)
plot(allEffects(hyp_ek_model))


#ggplot the data in a boxplot
hypnea %>% ggplot(aes(treatment_graph, ek.1)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 3, aes(color = Temperature), position = "jitter", show.legend = FALSE) + 
  labs(x="salinity/nitrate", y= "Ek (μmols photons m-2 s-1)", title= "D", subtitle = "Hypnea musciformis") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  ylim(-1, 250) + stat_mean() + 
  scale_color_manual(values = c("#9C0627", "#BB589F", "#F4B4E2")) +
  geom_hline(yintercept=50, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#for Hypnea remove hm6-4 on 11/12 that had no d9 RLC (final weight 0.1017)
# and hm6-4 on 10/9/21 because it was white and also looked dead 
gr_hypnea <- subset(growth_rate, Species == "Hm" & final.weight != 0.1017 & growth_rate_percent > -87.96837)
hypnea$growth_rate <- round((gr_hypnea$final.weight - gr_hypnea$Initial.weight) / gr_hypnea$Initial.weight * 100, digits = 2)

#plot a regression between the photosynthetic independent variables of interest and growth rate
hypnea_growth_ek_graph <- ggplot(hypnea, aes(x=ek.1, y=growth_rate)) + 
  geom_point(alpha = 0.5, size = 3, aes(color = Treatment), show.legend = TRUE) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Hypnea musciformis Ek vs Growth Rate", x = "Ek (μmols photons m-2 s-1)", 
       y = "growth rate (%)") + stat_regline_equation(label.x = 25, label.y = 165) + stat_cor(label.y = 150)
hypnea_growth_ek_graph

#summarize the means
hypnea %>% group_by(Treatment) %>% summarise_at(vars(ek.1), list(mean = mean))
hypnea %>% group_by(Temperature) %>% summarise_at(vars(ek.1), list(mean = mean))




#no longer needed--------------------------------------------------------------
#check for equal variance
bartlett.test(ek.1 ~ Treatment, data = ulva)
bartlett.test(ek.1 ~ Temperature, data = ulva)

#run Welch's ANOVA if not equal variance
welch_anova_treatment <- oneway.test(ek.1 ~ Treatment, data = ulva, var.equal = FALSE)
welch_anova_treatment
welch_anova_temp <- oneway.test(ek.1 ~ Temperature, data = ulva, var.equal = TRUE)
welch_anova_temp
games_howell_test(ulva, ek.1 ~ Treatment, conf.level = 0.95, detailed = TRUE)
games_howell_test(ulva, ek.1 ~ Temperature, conf.level = 0.95, detailed = TRUE)

#check for equal variance
bartlett.test(ek.1 ~ Treatment, data = hypnea)
bartlett.test(ek.1 ~ Temperature, data = hypnea)
#run Welch's ANOVA if not equal variance 
#welch_anova_treatment <- oneway.test(ek.1 ~ Treatment, data = hypnea, var.equal = FALSE)
#welch_anova_treatment
welch_anova_temp <- oneway.test(ek.1 ~ Temperature, data = hypnea, var.equal = FALSE)
welch_anova_temp
games_howell_test(hypnea, ek.1 ~ Treatment, conf.level = 0.95, detailed = TRUE)

