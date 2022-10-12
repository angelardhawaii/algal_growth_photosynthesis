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


run5_6_photosyn_data <- read.csv("data_input/run5-6_ek_alpha.csv")

# assign run as a factor
run5_6_photosyn_data$Run <- as.factor(run5_6_photosyn_data$Run)

#assign temperature as a factor
run5_6_photosyn_data$Temperature <- as.factor(run5_6_photosyn_data$Temp...C.)

#assigns treatment as characters from integers then to factors
run5_6_photosyn_data$Treatment <- as.factor(as.character(run5_6_photosyn_data$Treatment))

# assign deltaNPQ as a factor
run5_6_photosyn_data$deltaNPQ <- as.factor(run5_6_photosyn_data$deltaNPQ)

 
#toggle between the species for output. Use Day 9 for final analysis
#recent change: removing the very odd ek.1 value of 559.4 in hypnea dataset
hypnea <- subset(run5_6_photosyn_data, Species == "hm" & RLC.Day == 9 & ek.1 < 226)
hypnea$treatment_graph[hypnea$Treatment == 1] <- "1) 35ppt/14umol" 
hypnea$treatment_graph[hypnea$Treatment == 2] <- "2) 28ppt/27umol" 
hypnea$treatment_graph[hypnea$Treatment == 3] <- "4) 18ppt/53umol" 
hypnea$treatment_graph[hypnea$Treatment == 4] <- "5) 11ppt/80umol"
hypnea$treatment_graph[hypnea$Treatment == 5] <- "3) 28ppt/53umol"

ulva <- subset(run5_6_photosyn_data, Species == "ul" & RLC.Day == 9)
ulva$treatment_graph[ulva$Treatment == 0] <- "1) 35ppt/0.5umol"
ulva$treatment_graph[ulva$Treatment == 1] <- "2) 35ppt/14umol" 
ulva$treatment_graph[ulva$Treatment == 2] <- "3) 28ppt/27umol" 
ulva$treatment_graph[ulva$Treatment == 3] <- "4) 18ppt/53umol" 
ulva$treatment_graph[ulva$Treatment == 4] <- "5) 11ppt/80umol"
#add growth rate from other dataset to this one and subset by species
growth_rate <- read.csv("/Users/Angela/src/work/limu/algal_growth_photosynthesis/data_input/all_growth_081722.csv")
growth_rate$Species <- as.factor(growth_rate$Species)
growth_rate$treatment <- as.factor(growth_rate$treatment)
growth_rate
gr_ulva <- subset(growth_rate, Species == "Ul")

#for Hypnea remove the same replicate that was removed in the previous dataset (hm3-2 on 11/9) 
gr_hypnea <- subset(growth_rate, Species == "Hm" & final.weight != 0.1017 & final.weight != 0.3224)

ulva$growth_rate <- round((gr_ulva$final.weight - gr_ulva$Initial.weight) / gr_ulva$Initial.weight * 100, digits = 2)
hypnea$growth_rate <- round((gr_hypnea$final.weight - gr_hypnea$Initial.weight) / gr_hypnea$Initial.weight * 100, digits = 2)

# run model with interaction between temperature and treatment for rETRmax (is best?)
#run5_6_photosyn_model <- lmer(formula = rETRmax ~ Treatment * Temperature + (1 | run), data = ulva)
#check the performance of the model for each dataset. If collinearity good, proceed with this model
#if not good, use model with no interaction between the fixed effect variables
#performance::check_model(run5_6_photosyn_model, data = ulva)


#ULVA 
#run model without interaction between the treatments and temperature
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

#check for equal variance
bartlett.test(rETRmax ~ Treatment, data = ulva)
#run Welch's ANOVA if not equal variance (p = 0.01426, not equal)
welch_anova_treatment <- oneway.test(rETRmax ~ Treatment, data = ulva, var.equal = FALSE)
welch_anova_treatment
games_howell_test(ulva, rETRmax ~ Treatment, conf.level = 0.95, detailed = TRUE)
#run ANOVA and pairwise comparisons
#anova(run5_6_photosyn_model_noint, type = c("III"), ddf = "Satterthwaite")
#ulva_photosyn_model_aov <- aov(rETRmax ~ Treatment + Temperature, data = ulva)
#TukeyHSD(ulva_photosyn_model_aov, "Treatment", ordered = FALSE)
#TukeyHSD(ulva_photosyn_model_aov, "Temperature", ordered = FALSE)

tab_model(run5_6_photosyn_model_noint)
plot(allEffects(run5_6_photosyn_model_noint))

###run model for Ulva and Ek
run5_6_photosyn_model_ek <- lmer(formula = ek.1 ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = ulva)
#make a histogram and residual plots of the data for ulva
hist(ulva$ek.1, main = paste("Ulva lactuca -- Ek"), col = "olivedrab3", labels = TRUE)
plot(resid(run5_6_photosyn_model_ek) ~ fitted(run5_6_photosyn_model_ek))
qqnorm(resid(run5_6_photosyn_model_ek))
qqline(resid(run5_6_photosyn_model_ek))

#check the performance of the model
performance::check_model(run5_6_photosyn_model_ek)

r.squaredGLMM(run5_6_photosyn_model_ek)
summary(run5_6_photosyn_model_ek)

#run ANOVA and pairwise comparisons
anova(run5_6_photosyn_model_ek, type = c("III"), ddf = "Satterthwaite")
ulva_photosyn_model_ek_aov <- aov(ek.1 ~ Treatment + Temperature, data = ulva)
TukeyHSD(ulva_photosyn_model_ek_aov, "Treatment", ordered = FALSE)
TukeyHSD(ulva_photosyn_model_ek_aov, "Temperature", ordered = FALSE)


plot(allEffects(run5_6_photosyn_model_ek))

#plot a regression between the photosynthetic independent variables of interest and growth rate
ulva_growth_etr_graph <- ggplot(ulva, aes(x=rETRmax, y=growth_rate)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE) + 
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
run5_6_photosyn_model_noint <- lmer(formula = rETRmax ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea, REML = FALSE)


#make a histogram and residual plots of the data for hypnea
hist(hypnea$rETRmax, main = paste("Hypnea musciformis"), col = "maroon", labels = TRUE)
plot(resid(run5_6_photosyn_model_noint) ~ fitted(run5_6_photosyn_model_noint))
qqnorm(resid(run5_6_photosyn_model_noint))
qqline(resid(run5_6_photosyn_model_noint))

#check the performance of the model
performance::check_model(run5_6_photosyn_model_noint)

r.squaredGLMM(run5_6_photosyn_model_noint)
summary(run5_6_photosyn_model_noint)

#check for equal variance
bartlett.test(rETRmax ~ Treatment, data = hypnea)
#run Welch's ANOVA if not equal variance (p = 0.01426, not equal)
welch_anova_treatment <- oneway.test(rETRmax ~ Treatment, data = hypnea, var.equal = FALSE)
welch_anova_treatment
games_howell_test(hypnea, rETRmax ~ Treatment, conf.level = 0.95, detailed = TRUE)
#run ANOVA and pairwise comparisons
anova(run5_6_photosyn_model_noint, type = c("III"), ddf = "Satterthwaite")
hypnea_photosyn_model_aov <- aov(rETRmax ~ Treatment + Temperature, data = hypnea)
TukeyHSD(hypnea_photosyn_model_aov, "Temperature", ordered = FALSE)



plot(allEffects(run5_6_photosyn_model_noint))

###run model for Hypnea and Ek
run5_6_photosyn_model_ek <- lmer(formula = ek.1 ~ Treatment + Temperature + (1 | Run) + (1 | Plant.ID) + (1 | RLC.Order), data = hypnea)
#make a histogram and residual plots of the data for hypnea
hist(hypnea$ek.1, main = paste("Hypnea musciformis -- Ek"), col = "maroon2", labels = TRUE)
plot(resid(run5_6_photosyn_model_ek) ~ fitted(run5_6_photosyn_model_ek))
qqnorm(resid(run5_6_photosyn_model_ek))
qqline(resid(run5_6_photosyn_model_ek))

#check the performance of the model
performance::check_model(run5_6_photosyn_model_ek)

r.squaredGLMM(run5_6_photosyn_model_ek)
summary(run5_6_photosyn_model_ek)

#run ANOVA and pairwise comparisons
anova(run5_6_photosyn_model_ek, type = c("III"), ddf = "Satterthwaite")
hypnea_photosyn_model_ek_aov <- aov(ek.1 ~ Treatment + Temperature, data = hypnea)
TukeyHSD(hypnea_photosyn_model_ek_aov, "Treatment", ordered = FALSE)
TukeyHSD(hypnea_photosyn_model_ek_aov, "Temperature", ordered = FALSE)


plot(allEffects(run5_6_photosyn_model_ek))

#plot a regression between the photosynthetic independent variables of interest and growth rate
hypnea_growth_etr_graph <- ggplot(hypnea, aes(x=rETRmax, y=growth_rate)) + 
  geom_point(alpha = 0.5, size = 3, show.legend = TRUE) + 
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Hypnea musciformis rETRmax vs Growth Rate", x = "rETRmax (μmols electrons m-2 s-1)", 
       y = "growth rate (%)") + stat_regline_equation(label.x = 25, label.y = 165) + stat_cor(label.y = 150)
hypnea_growth_etr_graph

hypnea_growth_ek_graph <- ggplot(hypnea, aes(x=ek.1, y=growth_rate)) + 
  geom_point(alpha = 0.5, size = 3, aes(color = Treatment), show.legend = TRUE) +
  geom_smooth(method = "lm", col = "black") + theme_bw() + 
  labs(title = "Hypnea musciformis Ek vs Growth Rate", x = "Ek (μmols photons m-2 s-1)", y = "growth rate (%)") + 
  stat_regline_equation(label.y = 155) + stat_cor(label.y = 145)
hypnea_growth_ek_graph


#PLOTS
#rETRmax
ulva %>% ggplot(aes(rETRmax)) +
  geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

ulva %>% ggplot(aes(treatment_graph, rETRmax)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.5, size = 3, aes(color = Temperature), show.legend = TRUE) + 
  labs(x="salinty/nitrate", y= "rETR (umol e- m-2s-1)", title= "rETRmax", subtitle = "Ulva lactuca") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  stat_mean() + 
  geom_hline(yintercept=0, color = "purple", size = 0.5, alpha = 0.5) +
  theme_bw()


hypnea %>% ggplot(aes(treatment_graph, rETRmax)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.5, size = 3, aes(color = Temperature), show.legend = TRUE) + 
  labs(x="salinty/nitrate", y= "rETRmax (umol e- m-2s-1)", title= "rETRmax", subtitle = "Hypnea musciformis") + 
  scale_x_discrete(labels = c("35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  stat_mean() + 
  geom_hline(yintercept=0, color = "purple", size = 0.5, alpha = 0.5) +
  theme_bw()

#PLOTS
#Ek
ulva %>% ggplot(aes(treatment_graph, ek.1)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.5, size = 3, aes(color = Temperature), show.legend = TRUE) + 
  labs(x="salinty/nitrate", y= "Ek (umol photons m-2s-1)", title= "Ek", subtitle = "Ulva lactuca") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  stat_mean() + 
  geom_hline(yintercept=0, color = "purple", size = 0.5, alpha = 0.5) +
  theme_bw()


hypnea %>% ggplot(aes(treatment_graph, ek.1)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.5, size = 3, aes(color = Temperature), show.legend = TRUE) + 
  labs(x="salinty/nitrate", y= "Ek (umol photons m-2s-1)", title= "Ek", subtitle = "Hypnea musciformis") + 
  scale_x_discrete(labels = c("35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  stat_mean() + 
  geom_hline(yintercept=0, color = "purple", size = 0.5, alpha = 0.5) +
  theme_bw()
