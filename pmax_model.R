#script to run model for photosynthesis data -- similar to growth_rate script

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


#toggle between the species for output. Use Day 9 for final analysis
ulva <- subset(all_runs_photosyn_data, Species == "ul" & Treatment != 2.5)
ulva$treatment_graph[ulva$Treatment == 0] <- "1) 35ppt/0.5umol"
ulva$treatment_graph[ulva$Treatment == 1] <- "2) 35ppt/14umol" 
ulva$treatment_graph[ulva$Treatment == 2] <- "3) 28ppt/27umol" 
ulva$treatment_graph[ulva$Treatment == 3] <- "5) 18ppt/53umol" 
ulva$treatment_graph[ulva$Treatment == 4] <- "6) 11ppt/80umol"
ulva$treatment_graph[ulva$Treatment == 2.5] <- "4) 28ppt/53umol"


#ULVA pmax________________________________________________________________
ulva %>% ggplot(aes(pmax)) +
  geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

#run model without interaction between the treatments and temperature
ulva_pmax_model <- lmer(formula = pmax ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva)

#construct null model to perform likelihood ratio test REML must be FALSE
ulva_pmax_treatment_null <- lmer(formula = pmax ~ Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
ulva_pmax_model2 <- lmer(formula = pmax ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
anova(ulva_pmax_treatment_null, ulva_pmax_model2)
ulva_pmax_temperature_null <- lmer(formula = pmax ~ Treatment + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
ulva_pmax_model3 <- lmer(formula = pmax ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva, REML = FALSE)
anova(ulva_pmax_temperature_null, ulva_pmax_model3)


#make residual plots of the data for ulva
hist(resid(ulva_pmax_model))
plot(resid(ulva_pmax_model) ~ fitted(ulva_pmax_model))
qqnorm(resid(ulva_pmax_model))
qqline(resid(ulva_pmax_model))

#check the performance of the model
performance ::check_model(ulva_pmax_model)
r.squaredGLMM(ulva_pmax_model)

#make plots and tables for the data
tab_model(ulva_pmax_model, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(ulva_pmax_model))

ulva %>% ggplot(aes(treatment_graph, pmax)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 3, aes(color = Temperature), show.legend = TRUE) + 
  labs(x="salinity/nitrate", y= "Day 9 Pmax (μmols electrons m-2 s-1)", title= "A", subtitle = "Ulva lactuca") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  ylim(-1, 150) + stat_mean() + 
  scale_color_manual(values = c("#295102", "#7CB950", "#BDE269")) +
  geom_hline(yintercept=50, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#summarize the means for pmax
ulva %>% group_by(Treatment) %>% summarise_at(vars(pmax), list(mean = mean))
ulva %>% group_by(Temperature) %>% summarise_at(vars(pmax), list(mean = mean))
ulva %>% group_by(Treatment, RLC.Day) %>% summarise_at(vars(pmax), list(mean = mean))
ulva %>% group_by(Treatment, RLC.Day) %>% summarise_at(vars(rETRmax), list(mean = mean))

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

#plot a regression between the photosynthetic independent variables of interest and growth rate
ulva_growth_etr_graph <- ggplot(ulva, aes(x=pmax, y=growth_rate)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = Treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Ulva lactuca Pmax vs Growth Rate", x = "Pmax (μmols electrons m-2 s-1)", 
       y = "growth rate (%)") + stat_regline_equation(label.x = 25, label.y = 165) + stat_cor()
ulva_growth_etr_graph


#_____________________________________________________________________________________________________________
#HYPNEA 

#There was no D9 RLC for hm6-4 on 11/12/21 but had to remove hm6-4 from 10/9/21 below to match growth data
hypnea <- subset(all_runs_photosyn_data, Species == "hm" & uid != "2021-10-09_hm6-4")
hypnea$treatment_graph[hypnea$Treatment == 0] <- "1) 35ppt/0.5umol"
hypnea$treatment_graph[hypnea$Treatment == 1] <- "2) 35ppt/14umol" 
hypnea$treatment_graph[hypnea$Treatment == 2] <- "3) 28ppt/27umol" 
hypnea$treatment_graph[hypnea$Treatment == 3] <- "5) 18ppt/53umol" 
hypnea$treatment_graph[hypnea$Treatment == 4] <- "6) 11ppt/80umol"
hypnea$treatment_graph[hypnea$Treatment == 2.5] <- "4) 28ppt/53umol"

hypnea$temperature_graph[hypnea$Temp...C. == 20] <- "20°C"
hypnea$temperature_graph[hypnea$Temp...C. == 27] <- "27°C"
hypnea$temperature_graph[hypnea$Temp...C. == 30] <- "30°C"


#for Hypnea remove hm6-4 on 11/12 that had no d9 RLC (final weight 0.1017)
# and hm6-4 on 10/9/21 because it was white and also looked dead 
gr_hypnea <- subset(growth_rate, Species == "Hm" & final.weight != 0.1017 & growth_rate_percent > -87.96837)
hypnea$growth_rate <- round((gr_hypnea$final.weight - gr_hypnea$Initial.weight) / gr_hypnea$Initial.weight * 100, digits = 2)

#hm1-1 on 10/12/22 and hm1-2 on 4/29/22 causing issues of influential observations
#hm1-2 for both pmax and Ek -- leaving them in dataset because no good reason to believe not good data
#make a histogram and residual plots of the data for hypnea
hist(hypnea$pmax, main = paste("Hypnea musciformis"), col = "maroon", labels = TRUE)
hypnea %>% ggplot(aes(pmax)) +
  geom_histogram(binwidth=5, fill = "#7D0033", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

#run model for pmax
hyp_pmax_model <- lmer(formula = pmax ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea)

hist(resid(hyp_pmax_model))
plot(resid(hyp_pmax_model) ~ fitted(hyp_pmax_model))
qqnorm(resid(hyp_pmax_model))
qqline(resid(hyp_pmax_model))

#check the performance of the model
performance::check_model(hyp_pmax_model)
r.squaredGLMM(hyp_pmax_model)
summary(hyp_pmax_model)

tab_model(hyp_pmax_model, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)

plot(allEffects(hyp_pmax_model))

#construct null model to perform likelihood ratio test REML must be FALSE
hypnea_pmax_treatment_null <- lmer(formula = pmax ~ Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea, REML = FALSE)
hypnea_pmax_model2 <- lmer(formula = pmax ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea, REML = FALSE)
anova(hypnea_pmax_treatment_null, hypnea_pmax_model2)
hypnea_pmax_temperature_null <- lmer(formula = pmax ~ Treatment + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea, REML = FALSE)
hypnea_pmax_model3 <- lmer(formula = pmax ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea, REML = FALSE)
anova(hypnea_pmax_temperature_null, hypnea_pmax_model3)

hypnea %>% ggplot(aes(treatment_graph, pmax)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 3, aes(color = Temperature), show.legend = TRUE) + 
  labs(x="Treatment", y= "Day 9 Pmax (μmols electrons m-2 s-1)", title= "B", subtitle = "Hypnea musciformis") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  ylim(-1, 150) + stat_mean() + 
  scale_color_manual(values = c("#9C0627", "#BB589F", "#F4B4E2")) +
  geom_hline(yintercept=50, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#plot temperature
hypnea %>% ggplot(aes(temperature_graph, pmax)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 3, aes(color = Treatment), show.legend = TRUE) + 
  labs(x="Temperature", y= "Day 9 Pmax (μmols electrons m-2 s-1)", title= "F", subtitle = "Hypnea musciformis") + 
  scale_x_discrete(labels = c("20°C", "27°C", "30°C")) + 
  ylim(-1, 150) + stat_mean() + 
  geom_hline(yintercept=50, color = "red", size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("#9C0627", "#7CB950", "#F4B4E2", "#3311AA", "#3D96AA", "#F5B600")) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#plot a regression between the photosynthetic independent variables of interest and growth rate
hypnea_growth_etr_graph <- ggplot(hypnea, aes(x=pmax, y=growth_rate)) + 
  geom_point(alpha = 0.5, size = 3, aes(color = Treatment), show.legend = TRUE) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Hypnea musciformis Pmax vs Growth Rate", x = "Pmax (μmols electrons m-2 s-1)", 
       y = "growth rate (%)") + stat_regline_equation(label.x = 25, label.y = 165) + stat_cor(label.y = 150)
hypnea_growth_etr_graph

#summarize the means
hypnea %>% group_by(Treatment) %>% summarise_at(vars(pmax), list(mean = mean))
hypnea %>% group_by(Treatment, RLC.Day) %>% summarise_at(vars(pmax), list(mean = mean))
hypnea %>% group_by(Treatment, RLC.Day) %>% summarise_at(vars(rETRmaxYpoint1), list(mean = mean))
hypnea %>% group_by(Temperature) %>% summarise_at(vars(pmax), list(mean = mean))

#plot a regression between pmax and ek
hypnea_etr_ek_graph <- ggplot(hypnea, aes(x=pmax, y=ek.1)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE, aes(color = Treatment)) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Hypnea musciformis Pmax vs Ek", x = "Pmax (μmols electrons m-2 s-1)", 
       y = "Ek (μmols photons m-2 s-1)") + stat_regline_equation(label.x = 25, label.y = 165) + stat_cor()
hypnea_etr_ek_graph
#____________________________________________________________________________________


# No longer in use
#--------------------------------------------------------------------------------
#check for equal variance
bartlett.test(pmax ~ Treatment, data = ulva)
bartlett.test(pmax ~ Temperature, data = ulva)
#run Welch's ANOVA if not equal variance
welch_anova_treatment <- oneway.test(pmax ~ Treatment, data = ulva, var.equal = FALSE)
welch_anova_treatment
welch_anova_temp <- oneway.test(pmax ~ Temperature, data = ulva, var.equal = FALSE)
welch_anova_temp
games_howell_test(ulva, pmax ~ Treatment, conf.level = 0.95, detailed = TRUE)

#check for equal variance
bartlett.test(pmax ~ Treatment, data = hypnea)
bartlett.test(pmax ~ Temperature, data = hypnea)


#run Welch's ANOVA if not equal variance (p = 0.01426, not equal)
welch_anova_treatment <- oneway.test(pmax ~ Treatment, data = hypnea, var.equal = FALSE)
welch_anova_treatment
welch_anova_temp <- oneway.test(pmax ~ Temperature, data = hypnea, var.equal = FALSE)
welch_anova_temp
games_howell_test(hypnea, pmax ~ Treatment, conf.level = 0.95, detailed = TRUE)
games_howell_test(hypnea, pmax ~ Temperature, conf.level = 0.95, detailed = TRUE)


#run ANOVA and pairwise comparisons
#anova(all_runs_photosyn_model_noint, type = c("III"), ddf = "Satterthwaite")
#hypnea_photosyn_model_aov <- aov(pmax ~ Treatment + Temperature, data = hypnea)
#TukeyHSD(hypnea_photosyn_model_aov, "Temperature", ordered = FALSE)

#run ANOVA and pairwise comparisons
#anova(run5_6_photosyn_model_ek, type = c("III"), ddf = "Satterthwaite")
#hypnea_photosyn_model_ek_aov <- aov(ek.1 ~ Treatment + Temperature, data = hypnea)
#TukeyHSD(hypnea_photosyn_model_ek_aov, "Treatment", ordered = FALSE)
#TukeyHSD(hypnea_photosyn_model_ek_aov, "Temperature", ordered = FALSE)