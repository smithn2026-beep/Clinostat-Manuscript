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
RH_counts <- read.csv("C:/Users/smith/OneDrive - Ohio University/Clinostat Manuscript Materials/GitHub repository/Clinostat Manuscript/Mass + Root hair data/UPDATED_RH_COUNTS_CSV.csv")

library(ggplot2)
library(dplyr)
library(ggsignif)
RHgroup <- RH_counts%>%
  group_by(Length, Treatment) %>%
  summarize(mean=mean(Number,na.rm = TRUE),
            se=sd(Number,na.rm = TRUE)/sqrt(n()))

# RH graph------

ggplot(RH_counts, aes(x = Length, y = Number, fill = Treatment)) +
  geom_boxplot(position = position_dodge(0.7), width = 0.7) +
  geom_jitter(aes(shape = Treatment),position = position_jitterdodge(jitter.width = 0.15,
                                                                     dodge.width = 0.7),size = 2,alpha = 0.8) +
  geom_errorbar(data = RHgroup,inherit.aes = FALSE,
                aes(x = Length,y = mean,ymin = mean - se,ymax = mean + se,group = Treatment),
                position = position_dodge(0.7),width = 0.2) +
  scale_fill_manual(values = c("#696969", "lightgrey")) +
  xlab("Segment Length") + ylab("Number of Root Hairs") +
  theme_classic(base_size = 13) +
  geom_signif(y_position = c(300,310,285,170,160), xmin = c(.8,1.8,2.8,3.8,4.8),
              xmax = c(1.2,2.2,3.2,4.2,5.2), annotation = c("****","**","**","**","***"),
              tip_length = 0, textsize = 4)+
  theme_classic(base_size=13)
# two sample t test ------

RH_counts$Treatment <- as.character(RH_counts$Treatment)

# subset data by length------
zero_to_one_cm <- subset(RH_counts, Length == "0-1 cm")
one_to_two <- subset(RH_counts, Length == "1-2 cm")
two_to_three_cm <- subset(RH_counts, Length == "2-3 cm")
three_to_four_cm <- subset(RH_counts, Length == "3-4 cm")
four_to_five_cm <- subset(RH_counts, Length == "4-5 cm")

# subset the subsets by treatment-------
con01 <- subset(zero_to_one_cm, Treatment == "Control")
inoc01 <- subset(zero_to_one_cm, Treatment == "Inoculated")


con12 <- subset(one_to_two, Treatment == "Control")
inoc12 <- subset(one_to_two, Treatment == "Inoculated")


con23 <- subset(two_to_three_cm, Treatment == "Control")
inoc23 <- subset(two_to_three_cm, Treatment == "Inoculated")


con34 <- subset(three_to_four_cm, Treatment == "Control")
inoc34 <- subset(three_to_four_cm, Treatment == "Inoculated")


con45 <- subset(four_to_five_cm, Treatment == "Control")
inoc45 <- subset(four_to_five_cm, Treatment == "Inoculated")





# test for normality  0-1  # NORMAL-------
shapiro.test(con01$Number)
shapiro.test(inoc01$Number)

hist(con01$Number)
hist(inoc01$Number)


# test for equal variance 0-1
var.test(con01$Number, inoc01$Number) # FTR, F = 0.14253, num df = 6, denom df = 6, p-value = 0.03192

# Welch two sample t test 0-1 : t = -17.724, df = 7.6763, p-value = 1.668e-07


t.test(Number ~ Treatment, data=zero_to_one_cm)
t.test(con01$Number, inoc01$Number, var.equal = FALSE)

# doing both T tests amounts to the same data and is just done two different

# t = -15.974, df = 4.8288, p-value = 2.296e-05







# prelim test + t test for 1-2------
shapiro.test(con12$Number) #W = 0.86337, p-value = 0.2725
shapiro.test(inoc12$Number) #W = 0.80357, p-value = 0.1088

hist(con12$Number)
hist(inoc12$Number)


# test for equal variance 1-2
var.test(con12$Number, inoc12$Number) # F = 0.0016044, num df = 3, denom df = 3, p-value = 0.0002176

# Welch's t test for unequal variance: t = -9.343, df = 3.0096, p-value = 0.002563
t.test(Number ~ Treatment, data=one_to_two)
t.test(con12$Number, inoc12$Number, var.equal = FALSE)


# prelim test + t test 2-3-----
shapiro.test(con23$Number) #W = 0.94466, p-value = 0.683
shapiro.test(inoc23$Number) #W = 0.97342, p-value = 0.8625

hist(con23$Number)
hist(inoc23$Number)


# test for equal variance 1-2
var.test(con23$Number, inoc23$Number) # F = 0.0017937, num df = 3, denom df = 3, p-value = 0.0002571

# Welch's t test for unequal variance: t = -6.789, df = 3.0108, p-value = 0.006459
t.test(Number ~ Treatment, data=two_to_three_cm)
t.test(con23$Number, inoc23$Number, var.equal = FALSE)







# prelim test + t test 3-4------
shapiro.test(con34$Number) #W = 0.93788, p-value = 0.6414
shapiro.test(inoc34$Number) #W = 0.7864, p-value = 0.07998

hist(con34$Number)
hist(inoc34$Number)


# test for equal variance 3-4
var.test(con34$Number, inoc34$Number) # F = 0.055582, num df = 3, denom df = 3, p-value = 0.04037

# Welch's t test for unequal variance: t = -5.2315, df = 6, p-value = 0.001954
t.test(Number~Treatment, data=three_to_four_cm, var.equal= TRUE)
t.test(con34$Number, inoc34$Number, var.equal= TRUE)




# prelim test + t test 4-5------
shapiro.test(con45$Number) #W = 0.8494, p-value = 0.2242
shapiro.test(inoc45$Number) #W = 0.9151, p-value = 0.5098

hist(con45$Number)
hist(inoc45$Number)


# test for equal variance 4-5
var.test(con45$Number, inoc45$Number) #F = 0.011926, num df = 3, denom df = 3, p-value = 0.004329

# Wt test for equal variance: t = -8.9381, df = 6, p-value = 0.0001094
t.test(Number~Treatment, data=four_to_five_cm, var.equal= TRUE)
t.test(con45$Number, inoc45$Number, var.equal= TRUE)









# Shoot mass t tests and graphs-----
Mass_SSEP <- read.csv("C:/Users/smith/OneDrive - Ohio University/Clinostat Manuscript Materials/GitHub repository/Clinostat Manuscript/Mass + Root hair data/Mass Data.csv")
shoots <- subset(Mass_SSEP, Organ == "Shoot")

Mass_SSEP$Treatment <- as.factor(Mass_SSEP$Treatment)
Mass_SSEP$Organ <- as.factor(Mass_SSEP$Organ)

con_shoot <- subset(shoots, Treatment == "Control")

inoc_shoot <- subset(shoots, Treatment == "Inoculated")


# test for shoots -------
shapiro.test(con_shoot$Mass) # norm
shapiro.test(inoc_shoot$Mass) # norm

hist(con_shoot$Mass)
hist(inoc_shoot$Mass)


# test for equal variance - equal var
var.test(con_shoot$Ind_Mass, inoc_shoot$Ind_Mass) # F = 2.5329, num df = 15, denom df = 15, p-value = 0.08184

#  two sample t test 

t.test(con_shoot$Ind_Mass, inoc_shoot$Ind_Mass, var.equal= TRUE) #t = -6.0003, df = 30, p-value = 1.393e-06

#graph
ggplot(data = shoots, aes(x = Treatment, y = Mass, fill = Treatment)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, color = "black") +
  geom_jitter(width = 0.1, size = 3) +
  scale_fill_manual(values = c("Control" = "#696969", "Inoculated" = "lightgrey")) +
  labs(x = "Treatment", y = "Average Mass (g)") +
  geom_signif(y_position = c(.082), xmin = c(1.9), xmax = c(1.1), annotation = c("****"),tip_length = 0, textsize = 4)+
  theme(legend.position = "none")+
  theme_classic(base_size=13)


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
