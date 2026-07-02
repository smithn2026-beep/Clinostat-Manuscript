
library(dplyr)
library(car)


# DR5 ground control analysis (control vs bacteria) from still clinostat
df_DR5_still<-read.csv("C:/Users/smith/OneDrive - Ohio University/literally everything/Thesis/GFP experiment/DR5_still control/r analysis/df_dr5_ground_controls_entire.csv")
df_DR5_still$ID <- as.factor(df_DR5_still$ID)
df_DR5_still$Group <- as.factor(df_DR5_still$Group)

df_norm_dr5still <- df_DR5_still %>%
  group_by(ID) %>%
  mutate(RelativePosition = (Distance - min(Distance)) / (max(Distance) - min(Distance))) %>%
  ungroup()

# df_binned_PINClino <- df_norm_PINClino %>%
#   mutate(Bin = case_when(
#     RelativePosition <= 0.125 ~ "1",
#     RelativePosition <= 0.25 ~ "2",
#     RelativePosition <= 0.375 ~ "3",
#     RelativePosition <= 0.50 ~ "4",
#     RelativePosition <= 0.625 ~ "5",
#     RelativePosition <= 0.75 ~ "6",
#     RelativePosition <= 0.875 ~ "7",
#      TRUE ~ "8"))



#1 is the left most part of the root, 4 is left center, 5 is right center, 8 is the far right
df_binned_dr5still <- df_norm_dr5still %>%
  mutate(Bin = case_when(
    RelativePosition <= 0.25 ~ "1",
    RelativePosition <= 0.5 ~ "2",
    RelativePosition <= 0.75 ~ "3",
    TRUE ~ "4"))

#df_binned$Bin <- as.factor(df_binned$Bin)
#df_binned[2,2] <- "control"

image_bin_averages_dr5still <- df_binned_dr5still %>%
  group_by(ID, Group, Bin) %>%
  summarise(AvgGrayValue = mean(Intensity, na.rm = TRUE), .groups = "drop")





image_bin_averages_dr5still$Group <- factor(image_bin_averages_dr5still$Group)
#image_bin_averages_PINCLino$Bin <- factor(image_bin_averages_PINCLino$Bin, levels = c("Left 4", "Left 3", "Left 2", "Left 1","Right 1", "Right 2", "Right 3", "Right 4"))

#graph
library(ggplot2)
ggplot(data=image_bin_averages_dr5still,aes(x=Bin,y=AvgGrayValue, color=Group))+
  geom_point(show.legend=TRUE)+
  xlab("Bin")+
  ylab("Ratio of GFP")
image_bin_averages_dr5still$factorBin  <- as.factor(image_bin_averages_dr5still$Bin)
ggplot(data=image_bin_averages_dr5still,aes(x=factorBin,y=AvgGrayValue, color=Group))+
  geom_boxplot(show.legend=TRUE)+
  xlab("Bin")+
  ylab("Average Gray Value")

#assumptions
#making the model 
dr5stilllm<-lm(AvgGrayValue ~ Group + Bin, data =image_bin_averages_dr5still)

# normality
shapiro.test(residuals(dr5stilllm))

qqnorm(residuals(dr5stilllm))
qqline(residuals(dr5stilllm))

#var
boxplot(residuals(dr5stilllm) ~ image_bin_averages_dr5still$Group)
boxplot(residuals(dr5stilllm) ~ image_bin_averages_dr5still$Bin)
boxplot(residuals(dr5stilllm) ~ image_bin_averages_dr5still$Group+image_bin_averages_dr5still$Bin)



anova_result_dr5still <- aov(AvgGrayValue ~ Group+Bin, data = image_bin_averages_dr5still)

summary(anova_result_dr5still)

Anova_TUK_PINstill <- TukeyHSD(aov(AvgGrayValue ~ Group + Bin, data = image_bin_averages_dr5still))


library(multcomp)
anova_result_dr5still <- aov(AvgGrayValue ~ Group + Bin, data = image_bin_averages_dr5still)
Anova_TUK_dr5still <- glht(anova_result_dr5still, linfct=mcp( Group ="Tukey"))
cld(Anova_TUK_dr5still) 

summary(Anova_TUK_PINClino)


TukeyHSD(aov(AvgGrayValue ~ Group + Bin, data = image_bin_averages_dr5still))
plot(TukeyHSD(aov(AvgGrayValue ~ Group * Bin, data = image_bin_averages_dr5still)))
