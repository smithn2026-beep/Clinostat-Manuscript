
library(ggplot2)
library(dplyr)
library(car)

dr5_finalfigs <- read.csv("C:/Users/smith/OneDrive - Ohio University/literally everything/Thesis/GFP experiment/DR5_clino moving/DR5_reorienting_all_R_csv.csv")


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




image_bin_averages_dr5$Bin <- as.factor(image_bin_averages_dr5$Bin)


image_bin_averages_dr5$Variables <- paste(image_bin_averages_dr5$Treatment, image_bin_averages_dr5$Group, sep="_")



ggplot(image_bin_averages_dr5, aes(x = as.factor(Bin), y = AvgGrayValue,color = Variables,))+
  geom_boxplot()+
  xlab("Bin")


DR5_bigmodel <- lm(AvgGrayValue ~ Bin *Variables, data =image_bin_averages_dr5)
summary(DR5_bigmodel)    


# normality----
shapiro.test(residuals(DR5_bigmodel))

qqnorm(residuals(DR5_bigmodel))
qqline(residuals(DR5_bigmodel))

#var
boxplot(residuals(DR5_bigmodel) ~ image_bin_averages_dr5$Group)
boxplot(residuals(DR5_bigmodel) ~ image_bin_averages_dr5$Bin)
boxplot(residuals(DR5_bigmodel) ~ image_bin_averages_dr5$Group+image_bin_averages_dr5$Bin)



anova_DR5_bigModel <- aov(AvgGrayValue ~ Bin *Variables, data =image_bin_averages_dr5)

summary(anova_DR5_bigModel)

TukeyHSD(aov(AvgGrayValue ~ Bin *Variables, data =image_bin_averages_dr5))


library(multcomp)
anova_DR5_bigModel <- aov(AvgGrayValue ~ Bin *Variables, data =image_bin_averages_dr5)
DR5_Anova_TUK_bigModel <- glht(anova_DR5_bigModel, linfct=mcp( Bin ="Tukey"))
cld(DR5_Anova_TUK_bigModel) 

summary(anova_DR5_bigModel)

plot(TukeyHSD(aov(AvgGrayValue ~ Bin *Variables, data =image_bin_averages_dr5)))

