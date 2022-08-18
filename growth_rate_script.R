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

#open weight dataset and make columns for growth rate from initial and final weights
all_growth <- read.csv("/Users/Angela/src/work/limu/algal_growth_photosynthesis/data_input/all_runs_growth_081722.csv")

#make a new column for weight change (difference final from initial)
all_growth$growth_rate_percent <- (all_growth$final.weight - all_growth$Initial.weight) / all_growth$Initial.weight * 100

#make a new column for daily growth rate from 8 day study (steady growth rate assumed rather than exponential)
all_growth$steady_growth_daily <- all_growth$growth_rate_percent / 8

#make a new column that keeps only the numerical values (effectively removes the C in temperatures for consistency)
all_growth$temp_clean <- as.factor(substr(all_growth$temperature, 1, 2))

#assigns temperature as a factor
all_growth$temperature <- as.factor(all_growth$temp_clean)

#assigns treatment as characters from integers then to factors
all_growth$treatment <- as.factor(as.character(all_growth$treatment))

#assign run as a factor
all_growth$run <- as.factor(all_growth$run)

#assign plant ID as a factor
all_growth$plant.ID <- as.factor(all_growth$plant.ID)

#assign RLC order as a factor
all_growth$RLC.order <- as.factor(all_growth$RLC.order)

# assign lunar phase as factor
all_growth$lunar.phase <- as.factor(all_growth$lunar.phase)

#toggle between the species for output
hypnea <- subset(all_growth, Species == "Hm")
ulva <- subset(all_growth, Species == "Ul")

ulva <- subset(all_growth, Species == "Ul" & treatment != 0)
ulva$treatment_graph[ulva$treatment == 0] <- "1) 35ppt/0.5umol"
ulva$treatment_graph[ulva$treatment == 1] <- "2) 35ppt/14umol" 
ulva$treatment_graph[ulva$treatment == 2] <- "3) 28ppt/27umol" 
ulva$treatment_graph[ulva$treatment == 3] <- "4) 18ppt/53umol" 
ulva$treatment_graph[ulva$treatment == 4] <- "5) 11ppt/80umol"

hypnea <- subset(all_growth, Species == "Hm" & treatment != 0)
hypnea$treatment_graph[hypnea$treatment == 0] <- "1) 35ppt/0.5umol"
hypnea$treatment_graph[hypnea$treatment == 1] <- "2) 35ppt/14umol" 
hypnea$treatment_graph[hypnea$treatment == 2] <- "3) 28ppt/27umol" 
hypnea$treatment_graph[hypnea$treatment == 3] <- "5) 18ppt/53umol" 
hypnea$treatment_graph[hypnea$treatment == 4] <- "6) 11ppt/80umol"
hypnea$treatment_graph[hypnea$treatment == 5] <- "4) 28ppt/53umol"


#ULVA 
#run model without interaction (0 in model permits display of four levels of treatments - no intercept)
all_growth_model_noint <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | run) + (1 | plant.ID), data = ulva)

#make a histogram of the data for ulva
hist(ulva$growth_rate_percent, main = paste("Ulva lactuca Growth Rate (%)"), col = "olivedrab3", labels = TRUE)
#or
ulva %>% ggplot(aes(growth_rate_percent)) +
  geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()
plot(resid(all_growth_model_noint) ~ fitted(all_growth_model_noint))
qqnorm(resid(all_growth_model_noint))
qqline(resid(all_growth_model_noint))

#check the performance of the model for dataset: ulva
performance::check_model(all_growth_model_noint)

r.squaredGLMM(all_growth_model_noint)
summary(all_growth_model_noint)
#view random effects levels
ranef(all_growth_model_noint)

#check for equal variance
bartlett.test(growth_rate_percent ~ treatment, data = ulva)
#bartlett.test(growth_rate_percent ~ temperature, data = ulva)
#run Welch's ANOVA if not equal variance (p = 0.01426, not equal)
welch_anova_treatment <- oneway.test(growth_rate_percent ~ treatment, data = ulva, var.equal = FALSE)
welch_anova_treatment
games_howell_test(ulva, growth_rate_percent ~ treatment, conf.level = 0.95, detailed = TRUE)

#temperature has no effect

tab_model(all_growth_model_noint)
plot(allEffects(all_growth_model_noint))

ulva %>% ggplot(aes(treatment_graph, growth_rate_percent)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.5, size = 3, aes(color = temperature), show.legend = TRUE) + 
  labs(x="salinty/nitrate", y= "8-Day Growth Rate (%)", title= "Growth Rate (%)", subtitle = "Ulva lactuca") + 
  scale_x_discrete(labels = c("35ppt/14umolN", "28ppt/27umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  ylim(-100, 200) + stat_mean() + 
  geom_hline(yintercept=0, color = "purple", size = 0.5, alpha = 0.5) +
  theme_bw()

#HYPNEA 
#run model without interaction  (0 in model permits display of four levels of treatments - no intercept)
all_growth_model_noint <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | run) + (1 | plant.ID), data = hypnea)
#OR make a histogram of the data for hypnea
hist(hypnea$growth_rate_percent, main = paste("Hypnea musciformis Growth Rate (%)"), col = "maroon", labels = TRUE)

hypnea %>% ggplot(aes(growth_rate_percent)) +
  geom_histogram(binwidth=5, fill = "#990066", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()
plot(resid(all_growth_model_noint) ~ fitted(all_growth_model_noint))
qqnorm(resid(all_growth_model_noint))
qqline(resid(all_growth_model_noint))

#check the performance of the model for each dataset: ulva and hypnea
performance::check_model(all_growth_model_noint)

r.squaredGLMM(all_growth_model_noint)
summary(all_growth_model_noint)
ranef(all_growth_model_noint)

#check for equal variance
bartlett.test(growth_rate_percent ~ treatment, data = hypnea)
#bartlett.test(growth_rate_percent ~ temperature, data = hypnea)
#run Welch's ANOVA if not equal variance (p = 0.01426, not equal)
welch_anova_treatment <- oneway.test(growth_rate_percent ~ treatment, data = hypnea, var.equal = FALSE)
welch_anova_treatment
games_howell_test(hypnea, growth_rate_percent ~ treatment, conf.level = 0.95, detailed = TRUE)

#temperature has no effect

tab_model(all_growth_model_noint)
plot(allEffects(all_growth_model_noint))

hypnea %>% ggplot(aes(treatment_graph, growth_rate_percent)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.5, size = 3, aes(color = temperature), show.legend = TRUE) + 
  labs(x="salinty/nitrate", y= "8-Day Growth Rate (%)", title= "Growth Rate (%)", subtitle = "Hypnea musciformis") + 
  scale_x_discrete(labels = c("35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  ylim(-100, 200) + stat_mean() + 
  geom_hline(yintercept=0, color = "purple", size = 0.5, alpha = 0.5) +
  theme_bw()

plot(allEffects(all_growth_model_noint))
anova(all_growth_model_noint, type = c("III"), ddf = "Satterthwaite")
plot(resid(all_growth_model_noint) ~ fitted(all_growth_model_noint))
qqnorm(resid(all_growth_model_noint))
qqline(resid(all_growth_model_noint))
#--------------no longer in use------------------
#ulva_growth_model_aov <- aov(growth_rate_percent ~ treatment + temperature, data = ulva)
#TukeyHSD(ulva_growth_model_aov, "treatment", ordered = FALSE)
#hypnea_growth_model_aov <- aov(growth_rate_percent ~ treatment + temperature, data = hypnea)
#TukeyHSD(hypnea_growth_model_aov, "treatment", ordered = FALSE)
