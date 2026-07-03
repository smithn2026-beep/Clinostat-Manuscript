### CLINOSTAT ROOT CURVATURE ###----
#Clinostat Disorientation Root Curvature for Manuscript
install.packages("ggplot2")
install.packages("dplyr")
install.packages("readr")
install.packages("car")
install.packages("tidyr")
install.packages("nlme")
install.packages("emmeans")
install.packages("ggsignif")

library(ggplot2)
library(dplyr)
library(readr)
library(car)
library(tidyr)
library(nlme)
library(emmeans)
library(ggsignif)

Clinostat <- read.csv("C:/Users/smith/OneDrive - Ohio University/Clinostat Manuscript Materials/Clinostat Manuscript Materials/R code + figures/Clinostat_root curvature.csv")
#outlier test
justC <- Clinostat %>%
  filter(Treatment == "Control")
just1 <- Clinostat %>%
  filter(Treatment == "1 RPM/1 RPM")
justG <- Clinostat %>%
  filter(Treatment == "1.5 RPM/3.83 RPM")

ggplot(justC, aes(x=time, y = difference)) +geom_line()
ggplot(just1, aes(x=time, y = difference)) +geom_line()
ggplot(justG, aes(x=time, y = difference, color = as.factor(identity))) +geom_line()

#make graph

#find mean and standard error 
dataset <- Clinostat %>%
  group_by(Treatment, time) %>%
  summarise(mean = mean(difference, na.rm = TRUE),  
            se = sd(difference, na.rm = TRUE) / sqrt(n()))

rc.lm <- lm(difference ~ time + Treatment, data = Clinostat)
treatment.lm <- lm(difference ~ Treatment, data = Clinostat)
time.lm <- lm(difference ~ time, data = Clinostat)

shapiro.test(residuals(rc.lm))
#p<0.001
hist(residuals(treatment.lm))
#normal


leveneTest(difference ~ Treatment, data = Clinostat)
#Pr < 0.05
plot(residuals(time.lm))

plot(residuals(treatment.lm))

boxplot(difference ~ time, data = Clinostat)

boxplot(difference ~ Treatment, data = Clinostat)

remANOVA <- lme(difference ~ time + Treatment, random = ~1|identity, data = Clinostat)

lmeresid <- residuals(remANOVA)
boxplot(lmeresid~Clinostat$Treatment, main = "Genotype")

boxplot(lmeresid~Clinostat$identity, main = "Replicate")

plot(lmeresid~Clinostat$time, main = "Time")

shapiro.test(residuals(remANOVA))
#p<0.001
hist(residuals(remANOVA))
#normal

plot(residuals(remANOVA))


Clinostat$time = as.factor(Clinostat$time)

remANOVA <- lme(difference ~ time * Treatment, random = ~1|identity, 
                data = Clinostat)
Anova(remANOVA)

lsmeans(remANOVA, pairwise ~ Treatment|time)
#this output tells you significance 


#root curvature graph
ggplot(data = dataset, aes(x=time, y=mean, color = Treatment, shape=Treatment)) +  
  geom_point(size=5, stroke=1.5) +  
  geom_errorbar(data = dataset, aes(x=time, ymin=mean-se, ymax=mean+se), width=0.2) +  
  xlab("Time (hours)") +  
  ylab("Root Curvature (degrees)") +  
  geom_signif(y_position = c(-31, -38, -49, -60, -65, -67, -68, -67), xmin = c(2, 3, 4, 5, 6, 7, 8, 9), xmax = c(2, 3, 4, 5, 6, 7, 8, 9), annotation = c("*","***","****","****","****", "****", "****", "****"),
              tip_length = 0, textsize = 6)+
  theme_classic(base_size=13) +  
  geom_line(linewidth = 1.2, aes(linetype=Treatment)) +   
  geom_hline(yintercept = 0, color = "black")  +
  scale_color_manual(values=c('grey40','grey60', 'black'))+  
  theme(axis.text=element_text(size=15, color = 'black'),   
        axis.title=element_text(size=20, color = 'black'))+   
  scale_linetype_manual(values=c("longdash", "longdash", "solid")) +
  scale_x_continuous(limits = c(0,9.5), expand = c(0,0), breaks = c(0,1,2,3,4,5,6,7,8,9))+   
  scale_y_continuous(limits = c(-70,20), expand = c(0,3), breaks = c(-80,-60,-40,-20,0,20,40,60))   

### ROOT HAIR + MASS ###-----

### PIN1:GFP ###-----
library(ggplot2)
library(dplyr)

pin1_aspb <- read.csv("C:/Users/smith/OneDrive - Ohio University/Clinostat Manuscript Materials/Clinostat Manuscript Materials/R code + figures/GFP/PIN1_manuscript_data.csv")


pin1_aspb$Bin <- as.factor(pin1_aspb$Bin)

pin1_aspb$Variables <- paste(pin1_aspb$Treatment, pin1_aspb$group, sep="_")

bigmodel <- lm(AvgGrayValue ~ Bin *Variables, data =pin1_aspb)
summary(bigmodel)    

# normality----
shapiro.test(residuals(bigmodel))

qqnorm(residuals(bigmodel))
qqline(residuals(bigmodel))

#var
boxplot(residuals(bigmodel) ~ pin1_aspb$group)
boxplot(residuals(bigmodel) ~ pin1_aspb$Bin)
boxplot(residuals(bigmodel) ~ pin1_aspb$group+pin1_aspb$Bin)

#ANOVA
anova_bigModel <- aov(AvgGrayValue ~ Bin *Variables, data =pin1_aspb)

summary(anova_bigModel)

#post-hoc
TukeyHSD(aov(AvgGrayValue ~ Bin *Variables, data =pin1_aspb))


#Graph
library(dplyr)

pin1_graph <- pin1_aspb %>%
  filter(Variables != "ground_bacteria") %>%
  mutate(Variables = recode(Variables, "clinostat_bacteria" = "Clinostat bacteria",
                            "clinostat_control" = "Clinostat control",
                            "ground_control" = "Stationary control",
                            "ground_bacteria" = "Stationary bacteria"))

library(ggpattern)

ggplot(pin1_graph, aes(x = as.factor(Bin), y = AvgGrayValue, fill = Variables, pattern = Variables))+
  geom_boxplot_pattern(pattern_fill = "#808080")+ 
  labs(x = "Bin", y = "Average (grey scale value)") +
  scale_fill_manual(values = c("Clinostat bacteria" = "white", 
                               "Clinostat control" = "#696969",
                               "Stationary control" = "lightgrey" )) +
  scale_pattern_manual(values = c("Clinostat bacteria" = "none",
                                  "Clinostat control" = "none",
                                  "Stationary control" = "stripe"))+
  theme_classic()


### DR5:GFP ###-----

library(ggplot2)
library(dplyr)
library(car)

dr5_finalfigs <- read.csv("C:/Users/smith/OneDrive - Ohio University/Clinostat Manuscript Materials/Clinostat Manuscript Materials/R code + figures/GFP/DR5_manuscript_data.csv")


dr5_finalfigs$ID <- as.factor(dr5_finalfigs$ID)
dr5_finalfigs$Group <- as.factor(dr5_finalfigs$Group)

df_norm_dr5 <- dr5_finalfigs %>%
  group_by(ID) %>%
  mutate(RelativePosition = (Distance - min(Distance)) / (max(Distance) - min(Distance))) %>%
  ungroup()

dr5_finalfigs_binned <- df_norm_dr5 %>%
  mutate(Bin = case_when(
    RelativePosition <= 0.25 ~ "1",
    RelativePosition <= 0.5 ~ "2",
    RelativePosition <= 0.75 ~ "3",
    TRUE ~ "4"))

image_bin_averages_dr5 <- dr5_finalfigs_binned %>%
  group_by(ID, Group, Bin, Treatment) %>%
  summarise(AvgGrayValue = mean(Intensity, na.rm = TRUE), .groups = "drop")

dr5_finalfigs_binned$combined_treatment <- paste(dr5_finalfigs_binned$Treatment, dr5_finalfigs_binned$Group, sep="_")

image_bin_averages_dr5$Bin <- as.factor(image_bin_averages_dr5$Bin)


image_bin_averages_dr5$Variables <- paste(image_bin_averages_dr5$Treatment, image_bin_averages_dr5$Group, sep="_")


DR5_bigmodel <- lm(AvgGrayValue ~ Bin *Variables, data =image_bin_averages_dr5)
summary(DR5_bigmodel)    


# normality
shapiro.test(residuals(DR5_bigmodel))

qqnorm(residuals(DR5_bigmodel))
qqline(residuals(DR5_bigmodel))

#var
boxplot(residuals(DR5_bigmodel) ~ image_bin_averages_dr5$Group)
boxplot(residuals(DR5_bigmodel) ~ image_bin_averages_dr5$Bin)
boxplot(residuals(DR5_bigmodel) ~ image_bin_averages_dr5$Group+image_bin_averages_dr5$Bin)


#ANOVA
anova_DR5_bigModel <- aov(AvgGrayValue ~ Bin *Variables, data =image_bin_averages_dr5)

summary(anova_DR5_bigModel)

#post-hoc
TukeyHSD(aov(AvgGrayValue ~ Bin *Variables, data =image_bin_averages_dr5))
plot(TukeyHSD(aov(AvgGrayValue ~ Bin *Variables, data =image_bin_averages_dr5)))

summary(anova_DR5_bigModel)


#Graph 1
library(dplyr)
no_clino_con <- dr5_finalfigs_binned %>%
  filter(combined_treatment != "clinostat_control") %>%
  mutate(combined_treatment = dplyr::recode(combined_treatment, "clinostat_bacteria" = "Clinostat bacteria",
                                            "stationary_control" = "Stationary control",
                                            "stationary_bacteria" = "Stationary bacteria"))

ggplot(no_clino_con, aes(x = as.factor(Bin), y = Intensity, fill = combined_treatment, pattern = combined_treatment))+
  geom_boxplot_pattern(pattern_fill = "#808080", pattern_size = .2, pattern_density = .3  , outlier.shape = NA)+ 
  labs(x = "Bin", y = "Average (grey scale value)") +
  scale_fill_manual(values = c("Clinostat bacteria" = "white", 
                               "Stationary bacteria" = "white",
                               "Stationary control" = "#c7c7c7" )) +
  scale_pattern_manual(values = c("Clinostat bacteria" = "none",
                                  "Stationary bacteria" = "circle",
                                  "Stationary control" = "stripe"))+
  theme_classic()

#Graph 2
no_clino_bac <- dr5_finalfigs_binned %>%
  filter(combined_treatment != "clinostat_bacteria") %>%
  mutate(combined_treatment = dplyr::recode(combined_treatment, "clinostat_control" = "Clinostat control",
                                            "stationary_control" = "Stationary control",
                                            "stationary_bacteria" = "Stationary bacteria"))

ggplot(no_clino_bac, aes(x = as.factor(Bin), y = Intensity, fill = combined_treatment, pattern = combined_treatment))+
  geom_boxplot_pattern(pattern_fill = "#808080", pattern_size = .2, pattern_density = .3  , outlier.shape = NA)+ 
  labs(x = "Bin", y = "Average (grey scale value)") +
  scale_fill_manual(values = c("Clinostat control" = "#696969", 
                               "Stationary bacteria" = "white",
                               "Stationary control" = "#c7c7c7" )) +
  scale_pattern_manual(values = c("Clinostat control" = "none",
                                  "Stationary bacteria" = "circle",
                                  "Stationary control" = "stripe"))+
  theme_classic()
