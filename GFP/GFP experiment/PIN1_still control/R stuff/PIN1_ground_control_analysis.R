
library(dplyr)
library(car)


# PIN1 ground control analysis (control vs bacteria) from still clinostat
df_PIN1_still<-read.csv("C:/Users/smith/OneDrive - Ohio University/literally everything/Thesis/GFP experiment/PIN1_still control/raw_data_PIN_still.csv")
df_PIN1_still$ID <- as.factor(df_PIN1_still$ID)
df_PIN1_still$group <- as.factor(df_PIN1_still$group)

df_norm_PINstill <- df_PIN1_still %>%
  group_by(ID) %>%
  mutate(RelativePosition = (distance - min(distance)) / (max(distance) - min(distance))) %>%
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
df_binned_still <- df_norm_PINstill %>%
  mutate(Bin = case_when(
    RelativePosition <= 0.25 ~ "1",
    RelativePosition <= 0.5 ~ "2",
    RelativePosition <= 0.75 ~ "3",
    TRUE ~ "4"))

#df_binned$Bin <- as.factor(df_binned$Bin)
#df_binned[2,2] <- "control"

image_bin_averages_PINstill <- df_binned_still %>%
  group_by(ID, group, Bin) %>%
  summarise(AvgGrayValue = mean(intensity, na.rm = TRUE), .groups = "drop")





image_bin_averages_PINstill$group <- factor(image_bin_averages_PINstill$group)
#image_bin_averages_PINCLino$Bin <- factor(image_bin_averages_PINCLino$Bin, levels = c("Left 4", "Left 3", "Left 2", "Left 1","Right 1", "Right 2", "Right 3", "Right 4"))

#graph
library(ggplot2)
ggplot(data=image_bin_averages_PINstill,aes(x=Bin,y=AvgGrayValue, color=group))+
  geom_point(show.legend=TRUE)+
  xlab("Treatment")+
  ylab("Ratio of GFP")
image_bin_averages_PINstill$factorBin  <- as.factor(image_bin_averages_PINstill$Bin)
ggplot(data=image_bin_averages_PINstill,aes(x=factorBin,y=AvgGrayValue, color=group))+
  geom_boxplot(show.legend=TRUE)+
  xlab("Bin")+
  ylab("Average Gray Value")

#assumptions
#making the model 
PINstilllm<-lm(AvgGrayValue ~ group + Bin, data =image_bin_averages_PINstill)

# normality
shapiro.test(residuals(PINstilllm))

qqnorm(residuals(PINstilllm))
qqline(residuals(PINstilllm))

#var
boxplot(residuals(PINstilllm) ~ image_bin_averages_PINstill$group)
boxplot(residuals(PINstilllm) ~ image_bin_averages_PINstill$Bin)
boxplot(residuals(PINstilllm) ~ image_bin_averages_PINstill$group+image_bin_averages_PINstill$Bin)



anova_result_PINstill <- aov(AvgGrayValue ~ group+Bin, data = image_bin_averages_PINstill)

summary(anova_result_PINstill)

Anova_TUK_PINstill <- TukeyHSD(aov(AvgGrayValue ~ group + Bin, data = image_bin_averages_PINstill))


library(multcomp)
anova_result_PINstill <- aov(AvgGrayValue ~ group + Bin, data = image_bin_averages_PINstill)
Anova_TUK_PINstill <- glht(anova_result_PINstill, linfct=mcp( Bin ="Tukey"))
cld(Anova_TUK_PIN) 

summary(Anova_TUK_PINClino)


TukeyHSD(aov(AvgGrayValue ~ group + Bin, data = image_bin_averages_PINCLino))
plot(TukeyHSD(aov(AvgGrayValue ~ group * Bin, data = image_bin_averages_PINCLino)))
