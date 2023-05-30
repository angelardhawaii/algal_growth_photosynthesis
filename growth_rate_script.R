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
#library(mmtable2)
library(gt)
library(purrr)
library(stringr)
library(tidyr)
library(piecewiseSEM)

#open weight dataset and make columns for growth rate from initial and final weights
all_growth <- read.csv("/Users/Angela/src/work/limu/algal_growth_photosynthesis/data_input/all_runs_growth_011723.csv")

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


#ULVA
#subset data by species
ulva <- subset(all_growth, Species == "Ul" & treatment != 2.5)

#create subsets for the plots
ulva$treatment_graph[ulva$treatment == 0] <- "1) 35ppt/0.5umol"
ulva$treatment_graph[ulva$treatment == 1] <- "2) 35ppt/14umol" 
ulva$treatment_graph[ulva$treatment == 2] <- "3) 28ppt/27umol" 
ulva$treatment_graph[ulva$treatment == 3] <- "5) 18ppt/53umol" 
ulva$treatment_graph[ulva$treatment == 4] <- "6) 11ppt/80umol"
#ulva$treatment_graph[ulva$treatment == 2.5] <- "4) 28ppt/53umol"


#make a histogram of the data for ulva
hist(ulva$growth_rate_percent, main = paste("Ulva lactuca Growth Rate (%)"), col = "olivedrab3", labels = TRUE)
#or
ulva %>% ggplot(aes(growth_rate_percent)) +
  geom_histogram(binwidth=5, fill = "#5BB300", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

#run model without interaction
growth_model_ulva <- lmer(formula = growth_rate_percent ~ treatment + temperature +
                                 (1 | plant.ID) + (1 | lunar.phase), data = ulva, REML = TRUE)

hist(resid(growth_model_ulva))
plot(resid(growth_model_ulva) ~ fitted(growth_model_ulva))
qqnorm(resid(growth_model_ulva))
qqline(resid(growth_model_ulva))

#check the performance of the model for dataset: ulva
performance::check_model(growth_model_ulva)
rsquared(growth_model_ulva)
summary(growth_model_ulva)
#view random effects levels
ranef(growth_model_ulva)
tab_model(growth_model_ulva, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(growth_model_ulva))

#construct null model to perform likelihood ratio test REML must be FALSE
ulva_growth_treatment_null <- lmer(formula = growth_rate_percent ~ temperature + (1 | lunar.phase) + (1 | plant.ID), data = ulva, REML = FALSE)
ulva_growth_model2 <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 |lunar.phase) + (1 | plant.ID), data = ulva, REML = FALSE)
anova(ulva_growth_treatment_null, ulva_growth_model2)
ulva_growth_temperature_null <- lmer(formula = growth_rate_percent ~ treatment + (1 | lunar.phase) + (1 | plant.ID), data = ulva, REML = FALSE)
ulva_growth_model3 <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | lunar.phase) + (1 | plant.ID), data = ulva, REML = FALSE)
anova(ulva_growth_temperature_null, ulva_growth_model3)

#plots
ulva %>% ggplot(aes(treatment_graph, growth_rate_percent)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 3, aes(color = temperature), show.legend = TRUE) + 
  labs(x="salinity/nitrate", y= "8-Day Growth (%)", title= "A", subtitle = "Ulva lactuca") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  ylim(-75, 200) + stat_mean() + 
  scale_color_manual(values = c("#295102", "#7CB950", "#BDE269")) +
  geom_hline(yintercept=4.46, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))

#summarize the means
ulva %>% group_by(treatment) %>% summarise_at(vars(growth_rate_percent), list(mean = mean))
ulva %>% group_by(temperature) %>% summarise_at(vars(growth_rate_percent), list(mean = mean))

#HYPNEA______________________________________________________________________
#removing hm6-4 on 11/12/21 because it was dead and had no D9 RLC (final weight 0.1017)
# and hm6-4 on 10/29/21 because it was white and also looked dead
hypnea <- subset(all_growth, Species == "Hm" & growth_rate_percent > -87.96837 & final.weight != 0.1017)
#for plots
hypnea$treatment_graph[hypnea$treatment == 0] <- "1) 35ppt/0.5umol"
hypnea$treatment_graph[hypnea$treatment == 1] <- "2) 35ppt/14umol" 
hypnea$treatment_graph[hypnea$treatment == 2] <- "3) 28ppt/27umol" 
hypnea$treatment_graph[hypnea$treatment == 3] <- "5) 18ppt/53umol" 
hypnea$treatment_graph[hypnea$treatment == 4] <- "6) 11ppt/80umol"
hypnea$treatment_graph[hypnea$treatment == 2.5] <- "4) 28ppt/53umol"

hypnea$temperature_graph[hypnea$temperature == 20] <- "20°C"
hypnea$temperature_graph[hypnea$temperature == 27] <- "27°C"
hypnea$temperature_graph[hypnea$temperature == 30] <- "30°C"


#make a histogram of the data for hypnea
hist(hypnea$growth_rate_percent, main = paste("Hypnea musciformis Growth Rate (%)"), col = "maroon", labels = TRUE)

hypnea %>% ggplot(aes(growth_rate_percent)) +
  geom_histogram(binwidth=5, fill = "#990066", color = "black", size = 0.25, alpha = 0.85) +
  theme_bw()

#run model without interaction
growth_model_hypnea <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | plant.ID) +
                                 (1 | run) + (1 | RLC.order), data = hypnea, REML = TRUE)


#plot residuals for hypnea growth model
plot(resid(growth_model_hypnea) ~ fitted(growth_model_hypnea))
qqnorm(resid(growth_model_hypnea))
qqline(resid(growth_model_hypnea))

#check the performance of the model
performance::check_model(growth_model_hypnea)
r.squaredGLMM(growth_model_hypnea)
summary(growth_model_hypnea)
ranef(growth_model_hypnea)
tab_model(growth_model_hypnea, show.intercept = TRUE, show.se = TRUE, show.stat = TRUE, show.df = TRUE, show.zeroinf = TRUE)
plot(allEffects(growth_model_hypnea))

#construct null model to perform likelihood ratio test REML must be FALSE
hypnea_growth_treatment_null <- lmer(formula = growth_rate_percent ~ temperature + (1 | plant.ID) + (1 | run) + (1 | RLC.order), data = hypnea, REML = FALSE)
hypnea_growth_model2 <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | plant.ID) + (1 | run) + (1 | RLC.order), data = hypnea, REML = FALSE)
anova(hypnea_growth_treatment_null, hypnea_growth_model2)
hypnea_growth_temperature_null <- lmer(formula = growth_rate_percent ~ treatment + (1 | plant.ID) + (1 | run) + (1 | RLC.order), data = hypnea, REML = FALSE)
hypnea_growth_model3 <- lmer(formula = growth_rate_percent ~ treatment + temperature + (1 | plant.ID) + (1 | run) + (1 | RLC.order), data = hypnea, REML = FALSE)
anova(hypnea_growth_temperature_null, hypnea_growth_model3)

hypnea %>% ggplot(aes(treatment_graph, growth_rate_percent)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 3, aes(color = temperature), show.legend = TRUE) + 
  labs(x="salinity/nitrate", y= "8-Day Growth Rate (%)", title= "B", subtitle = "Hypnea musciformis") + 
  scale_x_discrete(labels = c("35ppt/0.5umolN", "35ppt/14umolN", "28ppt/27umolN", "28ppt/53umolN", "18ppt/53umolN", "11ppt/80umolN")) + 
  ylim(-75, 200) + stat_mean() + 
  scale_color_manual(values = c("#9C0627", "#BB589F", "#F4B4E2")) +
  geom_hline(yintercept=0, color = "red", size = 0.5, alpha = 0.5) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))
#plot temperature
hypnea %>% ggplot(aes(temperature_graph, growth_rate_percent)) + 
  geom_boxplot(size=0.5) + 
  geom_point(alpha = 0.75, size = 5, aes(color = treatment), show.legend = TRUE) + 
  labs(x="Temperature", y= "9-day Growth (%)", title= "F", subtitle = "Hypnea musciformis") + 
  scale_x_discrete(labels = c("20°C", "27°C", "30°C")) + 
  ylim(-1, 150) + stat_mean() + 
  geom_hline(yintercept=50, color = "red", size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c("#9C0627", "#7CB950", "#F4B4E2", "#3311AA", "#3D96AA", "#F5B600")) +
  theme_bw() +
  theme(legend.position = c(0.90,0.90), plot.title = element_text(face = "bold", vjust = -15, hjust = 0.05), 
        plot.subtitle = element_text(face = "italic", size = 14, vjust = -20, hjust = 0.05))


#summarize the means
hypnea %>% group_by(treatment) %>% summarise_at(vars(growth_rate_percent), list(mean = mean))
hypnea %>% group_by(temperature) %>% summarise_at(vars(growth_rate_percent), list(mean = mean))




# dont use
#anova(all_growth_model_noint, type = c("III"), ddf = "Satterthwaite")


#--------------no longer in use------------------
#ulva_growth_model_aov <- aov(growth_rate_percent ~ treatment + temperature, data = ulva)
#TukeyHSD(ulva_growth_model_aov, "treatment", ordered = FALSE)
#hypnea_growth_model_aov <- aov(growth_rate_percent ~ treatment + temperature, data = hypnea)
#TukeyHSD(hypnea_growth_model_aov, "treatment", ordered = FALSE)



#temperature has no significant effect

#check for equal variance
bartlett.test(growth_rate_percent ~ treatment, data = ulva)
bartlett.test(growth_rate_percent ~ temperature, data = ulva)
#run Welch's ANOVA if not equal variance
welch_anova_treatment <- oneway.test(growth_rate_percent ~ treatment, data = ulva, var.equal = FALSE)
welch_anova_treatment
welch_anova_temp <- oneway.test(growth_rate_percent ~ temperature, data = ulva, var.equal = FALSE)
welch_anova_temp
games_howell_test(ulva, growth_rate_percent ~ treatment, conf.level = 0.95, detailed = FALSE)

#check for equal variance
bartlett.test(growth_rate_percent ~ treatment, data = hypnea)
bartlett.test(growth_rate_percent ~ temperature, data = hypnea)
#run Welch's ANOVA if not equal variance (p = 0.01426, not equal)
welch_anova_treatment <- oneway.test(growth_rate_percent ~ treatment, data = hypnea, var.equal = FALSE)
welch_anova_treatment
welch_anova_temp <- oneway.test(growth_rate_percent ~ temperature, data = hypnea, var.equal = FALSE)
welch_anova_temp
games_howell_test(hypnea, growth_rate_percent ~ treatment, conf.level = 0.95, detailed = TRUE)

#temperature effect not significant