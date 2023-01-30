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
#this seems to be the better way to go (change also fitWebb input)
all_runs_photosyn_data <- read.csv("data_input/hyp_ulva_all_runs_ek_alpha_normalized.csv")

# assign run as a factor
all_runs_photosyn_data$Run <- as.factor(all_runs_photosyn_data$Run)

#assign temperature as a factor
all_runs_photosyn_data$Temperature <- as.factor(all_runs_photosyn_data$Temp...C.)

#assigns treatment as characters from integers then to factors
all_runs_photosyn_data$Treatment <- as.factor(as.character(all_runs_photosyn_data$Treatment))

# assign deltaNPQ as a factor
all_runs_photosyn_data$deltaNPQ <- as.factor(all_runs_photosyn_data$deltaNPQ)

 
#toggle between the species for output. Use Day 9 for final analysis
#There was no D9 RLC for hm6-4 on 11/12/21 but had to remove hm6-4 from 10/9/21 below to match growth data
hypnea <- subset(all_runs_photosyn_data, Species == "hm" & RLC.Day == 9 & uid != "2021-10-09_hm6-4")
hypnea$treatment_graph[hypnea$Treatment == 0] <- "1) 35ppt/0.5umol"
hypnea$treatment_graph[hypnea$Treatment == 1] <- "2) 35ppt/14umol" 
hypnea$treatment_graph[hypnea$Treatment == 2] <- "3) 28ppt/27umol" 
hypnea$treatment_graph[hypnea$Treatment == 3] <- "5) 18ppt/53umol" 
hypnea$treatment_graph[hypnea$Treatment == 4] <- "6) 11ppt/80umol"
hypnea$treatment_graph[hypnea$Treatment == 2.5] <- "4) 28ppt/53umol"

ulva <- subset(all_runs_photosyn_data, Species == "ul" & RLC.Day == 9)
ulva$treatment_graph[ulva$Treatment == 0] <- "1) 35ppt/0.5umol"
ulva$treatment_graph[ulva$Treatment == 1] <- "2) 35ppt/14umol" 
ulva$treatment_graph[ulva$Treatment == 2] <- "3) 28ppt/27umol" 
ulva$treatment_graph[ulva$Treatment == 3] <- "5) 18ppt/53umol" 
ulva$treatment_graph[ulva$Treatment == 4] <- "6) 11ppt/80umol"
ulva$treatment_graph[ulva$Treatment == 2.5] <- "4) 28ppt/53umol"


#add growth rate from other dataset to this one and subset by species
growth_rate <- read.csv("/Users/Angela/src/work/limu/algal_growth_photosynthesis/data_input/all_runs_growth_011723.csv")
growth_rate$Species <- as.factor(growth_rate$Species)
growth_rate$treatment <- as.factor(growth_rate$treatment)

#make a new column for weight change (difference final from initial)
growth_rate$growth_rate_percent <- (growth_rate$final.weight - growth_rate$Initial.weight) / growth_rate$Initial.weight * 100

#for Hypnea remove hm6-4 on 11/12 that had no d9 RLC (final weight 0.1017)
# and hm6-4 on 10/9/21 because it was white and also looked dead 
gr_ulva <- subset(growth_rate, Species == "Ul")
gr_hypnea <- subset(growth_rate, Species == "Hm" & final.weight != 0.1017 & growth_rate_percent > -87.96837)

ulva$growth_rate <- round((gr_ulva$final.weight - gr_ulva$Initial.weight) / gr_ulva$Initial.weight * 100, digits = 2)
hypnea$growth_rate <- round((gr_hypnea$final.weight - gr_hypnea$Initial.weight) / gr_hypnea$Initial.weight * 100, digits = 2)

# run model with interaction between temperature and treatment for rETRmax (is best?)
#run5_6_photosyn_model <- lmer(formula = rETRmax ~ Treatment * Temperature + (1 | run), data = ulva)
#check the performance of the model for each dataset. If collinearity good, proceed with this model
#if not good, use model with no interaction between the fixed effect variables
#performance::check_model(run5_6_photosyn_model, data = ulva)


#ULVA rETRmax

#run model without interaction between the treatments and temperature
all_runs_photosyn_model_noint <- lmer(formula = rETRmax ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva)


#make a histogram and residual plots of the data for ulva
hist(ulva$rETRmax, main = paste("Ulva lactuca"), col = "olivedrab3", labels = TRUE)
plot(resid(all_runs_photosyn_model_noint) ~ fitted(all_runs_photosyn_model_noint))
qqnorm(resid(all_runs_photosyn_model_noint))
qqline(resid(all_runs_photosyn_model_noint))

ulva %>% ggplot(aes(rETRmax)) +
  geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()


#check the performance of the model
performance::check_model(all_runs_photosyn_model_noint)

r.squaredGLMM(all_runs_photosyn_model_noint)
summary(all_runs_photosyn_model_noint)

#check for equal variance
bartlett.test(rETRmax ~ Treatment, data = ulva)

#run Welch's ANOVA if not equal variance
welch_anova_treatment <- oneway.test(rETRmax ~ Treatment, data = ulva, var.equal = FALSE)
welch_anova_treatment
welch_anova_temp <- oneway.test(rETRmax ~ Temperature, data = ulva, var.equal = FALSE)
welch_anova_temp
games_howell_test(ulva, rETRmax ~ Treatment, conf.level = 0.95, detailed = TRUE)

#run ANOVA and pairwise comparisons
#anova(run5_6_photosyn_model_noint, type = c("III"), ddf = "Satterthwaite")
#ulva_photosyn_model_aov <- aov(rETRmax ~ Treatment + Temperature, data = ulva)
#TukeyHSD(ulva_photosyn_model_aov, "Treatment", ordered = FALSE)
#TukeyHSD(ulva_photosyn_model_aov, "Temperature", ordered = FALSE)

#make plots and tables for the data
tab_model(all_runs_photosyn_model_noint)
plot(allEffects(all_runs_photosyn_model_noint))

ulva %>% ggplot(aes(treatment_graph, rETRmax)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.5, size = 3, aes(color = Temperature), show.legend = FALSE) + 
  labs(x="salinity/nitrate", y= "rETRmax (μmols electrons m-2 s-1)", title= "A", subtitle = "Ulva lactuca") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  ylim(-1, 170) + stat_mean() + 
  geom_hline(yintercept=0, color = "purple", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), plot.subtitle = element_text(face = "italic", vjust = -20, hjust = 0.05))

#summarize the means for rETRmax
ulva %>% group_by(Treatment) %>% summarise_at(vars(rETRmax), list(mean = mean))
ulva %>% group_by(Run) %>% summarise_at(vars(rETRmax), list(mean = mean))
#plot a regression between the photosynthetic independent variables of interest and growth rate
ulva_growth_etr_graph <- ggplot(ulva, aes(x=rETRmax, y=growth_rate)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = Treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Ulva lactuca rETRmax vs Growth Rate", x = "rETRmax (μmols electrons m-2 s-1)", 
       y = "growth rate (%)") + stat_regline_equation(label.x = 25, label.y = 165) + stat_cor()
ulva_growth_etr_graph

#plot a regression between rETRmax and Ek
ulva_growth_etr_ek_graph <- ggplot(ulva, aes(x=rETRmax, y=ek.1)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = Treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Ulva lactuca rETRmax vs Ek", x = "rETRmax (μmols electrons m-2 s-1)", 
       y = "Ek (μmols photons m-2 s-1)") + stat_regline_equation(label.x = 25, label.y = 165) + stat_cor()
ulva_growth_etr_ek_graph

#________________________________________________________________________________________

#ULva Ek model
all_runs_photosyn_model_ek <- lmer(formula = ek.1 ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva)

#make a histogram and residual plots of the data for ulva ek
hist(ulva$ek.1, main = paste("Ulva lactuca"), col = "olivedrab3", labels = TRUE)
plot(resid(all_runs_photosyn_model_ek) ~ fitted(all_runs_photosyn_model_ek))
qqnorm(resid(all_runs_photosyn_model_ek))
qqline(resid(all_runs_photosyn_model_ek))

ulva %>% ggplot(aes(ek.1)) +
  geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

#check the performance of the model
performance::check_model(all_runs_photosyn_model_ek)

r.squaredGLMM(all_runs_photosyn_model_ek)
summary(all_runs_photosyn_model_ek)

#check for equal variance
bartlett.test(ek.1 ~ Treatment, data = ulva)

#run Welch's ANOVA if not equal variance
welch_anova_treatment <- oneway.test(ek.1 ~ Treatment, data = ulva, var.equal = FALSE)
welch_anova_treatment
welch_anova_temp <- oneway.test(ek.1 ~ Temperature, data = ulva, var.equal = FALSE)
welch_anova_temp
games_howell_test(ulva, ek.1 ~ Treatment, conf.level = 0.95, detailed = TRUE)

#plot and tables
tab_model(all_runs_photosyn_model_ek)
plot(allEffects(all_runs_photosyn_model_ek))

ulva %>% ggplot(aes(treatment_graph, ek.1)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.5, size = 3, aes(color = Temperature), show.legend = FALSE) + 
  labs(x="salinity/nitrate", y= "Ek (μmols photons m-2 s-1)", title= "A", subtitle = "Ulva lactuca") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  ylim(-1, 200) + stat_mean() + 
  geom_hline(yintercept=0, color = "purple", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), plot.subtitle = element_text(face = "italic", vjust = -20, hjust = 0.05))

#plot linear regression between ek and growth
ulva_growth_ek_graph <- ggplot(ulva, aes(x=ek.1, y=growth_rate)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = Treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Ulva lactuca Ek vs Growth Rate", x = "Ek (μmols photons m-2 s-1)", 
       y = "growth rate (%)") + stat_regline_equation(label.x = 25, label.y = 165) + stat_cor()
ulva_growth_ek_graph

#summarize the means for Ek
ulva %>% group_by(Treatment) %>% summarise_at(vars(ek.1), list(mean = mean))

#_____________________________________________________________________________________________________________
#HYPNEA 
#hm1-1 on 10/12/22 and hm1-2 on 4/29/22 causing issues of influential observations
#hm1-2 for both rETRmax and Ek -- leaving them in dataset because no good reason to believe not good data

#run model without interaction
all_runs_photosyn_model_hyp <- lmer(formula = rETRmax ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea, REML = FALSE)


#make a histogram and residual plots of the data for hypnea
hist(hypnea$rETRmax, main = paste("Hypnea musciformis"), col = "maroon", labels = TRUE)
plot(resid(all_runs_photosyn_model_hyp) ~ fitted(all_runs_photosyn_model_hyp))
qqnorm(resid(all_runs_photosyn_model_hyp))
qqline(resid(all_runs_photosyn_model_hyp))

hypnea %>% ggplot(aes(rETRmax)) +
  geom_histogram(binwidth=5, fill = "#7D0033", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

#check the performance of the model
performance::check_model(all_runs_photosyn_model_hyp)

r.squaredGLMM(all_runs_photosyn_model_hyp)
summary(all_runs_photosyn_model_hyp)

#check for equal variance
bartlett.test(rETRmax ~ Treatment, data = hypnea)
#run Welch's ANOVA if not equal variance (p = 0.01426, not equal)
welch_anova_treatment <- oneway.test(rETRmax ~ Treatment, data = hypnea, var.equal = FALSE)
welch_anova_treatment
welch_anova_temp <- oneway.test(rETRmax ~ Temperature, data = hypnea, var.equal = FALSE)
welch_anova_temp
games_howell_test(hypnea, rETRmax ~ Treatment, conf.level = 0.95, detailed = TRUE)
games_howell_test(hypnea, rETRmax ~ Temperature, conf.level = 0.95, detailed = TRUE)

#run ANOVA and pairwise comparisons
#anova(all_runs_photosyn_model_noint, type = c("III"), ddf = "Satterthwaite")
#hypnea_photosyn_model_aov <- aov(rETRmax ~ Treatment + Temperature, data = hypnea)
#TukeyHSD(hypnea_photosyn_model_aov, "Temperature", ordered = FALSE)

tab_model(all_runs_photosyn_model_hyp)
plot(allEffects(all_runs_photosyn_model_hyp))

hypnea %>% ggplot(aes(treatment_graph, rETRmax)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.5, size = 3, aes(color = Temperature), show.legend = TRUE) + 
  labs(x="salinity/nitrate", y= "rETRmax (μmols electrons m-2 s-1)", title= "B", subtitle = "Hypnea musciformis") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  ylim(-1, 170) + stat_mean() + 
  geom_hline(yintercept=0, color = "purple", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), plot.subtitle = element_text(face = "italic", vjust = -20, hjust = 0.05))

#plot a regression between the photosynthetic independent variables of interest and growth rate
hypnea_growth_etr_graph <- ggplot(hypnea, aes(x=rETRmax, y=growth_rate)) + 
  geom_point(alpha = 0.5, size = 3, aes(color = Treatment), show.legend = TRUE) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Hypnea musciformis rETRmax vs Growth Rate", x = "rETRmax (μmols electrons m-2 s-1)", 
       y = "growth rate (%)") + stat_regline_equation(label.x = 25, label.y = 165) + stat_cor(label.y = 150)
hypnea_growth_etr_graph

#summarize the means
hypnea %>% group_by(Treatment) %>% summarise_at(vars(rETRmax), list(mean = mean))

#plot a regression between retrmax and ek
hypnea_etr_ek_graph <- ggplot(hypnea, aes(x=rETRmax, y=ek.1)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = Treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Hypnea musciformis rETRmax vs Ek", x = "rETRmax (μmols electrons m-2 s-1)", 
       y = "Ek (μmols photons m-2 s-1)") + stat_regline_equation(label.x = 25, label.y = 165) + stat_cor()
hypnea_etr_ek_graph
#____________________________________________________________________________________

###run model for Hypnea and Ek
all_runs_photosyn_model_ek <- lmer(formula = ek.1 ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID), data = hypnea)

#make a histogram and residual plots of the data for hypnea
hist(hypnea$ek.1, main = paste("Hypnea musciformis -- Ek"), col = "maroon2", labels = TRUE)
plot(resid(all_runs_photosyn_model_ek) ~ fitted(all_runs_photosyn_model_ek))
qqnorm(resid(all_runs_photosyn_model_ek))
qqline(resid(all_runs_photosyn_model_ek))

hypnea %>% ggplot(aes(ek.1)) +
  geom_histogram(binwidth=5, fill = "#7D0033", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()
#check the performance of the model
performance::check_model(all_runs_photosyn_model_ek)

r.squaredGLMM(all_runs_photosyn_model_ek)
summary(all_runs_photosyn_model_ek)

#check for equal variance
bartlett.test(ek.1 ~ Treatment, data = hypnea)
#run Welch's ANOVA if not equal variance 
welch_anova_treatment <- oneway.test(ek.1 ~ Treatment, data = hypnea, var.equal = FALSE)
welch_anova_treatment
welch_anova_temp <- oneway.test(ek.1 ~ Temperature, data = hypnea, var.equal = FALSE)
welch_anova_temp
games_howell_test(hypnea, ek.1 ~ Treatment, conf.level = 0.95, detailed = TRUE)

#run ANOVA and pairwise comparisons
#anova(run5_6_photosyn_model_ek, type = c("III"), ddf = "Satterthwaite")
#hypnea_photosyn_model_ek_aov <- aov(ek.1 ~ Treatment + Temperature, data = hypnea)
#TukeyHSD(hypnea_photosyn_model_ek_aov, "Treatment", ordered = FALSE)
#TukeyHSD(hypnea_photosyn_model_ek_aov, "Temperature", ordered = FALSE)

tab_model(all_runs_photosyn_model_ek)
plot(allEffects(all_runs_photosyn_model_ek))

hypnea %>% ggplot(aes(treatment_graph, ek.1)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.5, size = 3, aes(color = Temperature), show.legend = TRUE) + 
  labs(x="salinity/nitrate", y= "Ek (μmols photons m-2 s-1)", title= "B", subtitle = "Hypnea musciformis") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  ylim(-1, 240) + stat_mean() + 
  geom_hline(yintercept=0, color = "purple", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), plot.subtitle = element_text(face = "italic", vjust = -20, hjust = 0.05))

#plot a regression between the photosynthetic independent variables of interest and growth rate
hypnea_growth_ek_graph <- ggplot(hypnea, aes(x=ek.1, y=growth_rate)) + 
  geom_point(alpha = 0.5, size = 3, aes(color = Treatment), show.legend = TRUE) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Hypnea musciformis Ek vs Growth Rate", x = "Ek (μmols photons m-2 s-1)", 
       y = "growth rate (%)") + stat_regline_equation(label.x = 25, label.y = 165) + stat_cor(label.y = 150)
hypnea_growth_ek_graph
plot(allEffects(all_runs_photosyn_model_ek))

#summarize the means
hypnea %>% group_by(Treatment) %>% summarise_at(vars(ek.1), list(mean = mean))
hypnea %>% group_by(Run) %>% summarise_at(vars(ek.1), list(mean = mean))
