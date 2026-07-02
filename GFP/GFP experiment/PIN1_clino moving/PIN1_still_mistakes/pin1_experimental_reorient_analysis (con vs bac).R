
library(dplyr)
library(car)


# PIN1 reoriented analysis (control vs bacteria) from moving clinostat
df_PIN1_clino<-read.csv("C:/Users/smith/OneDrive - Ohio University/literally everything/Thesis/GFP experiment/PIN1_clino moving/pin1_data_R.csv")
df_PIN1_clino$ID <- as.factor(df_PIN1_clino$ID)
df_PIN1_clino$group <- as.factor(df_PIN1_clino$group)

df_norm_PINClino <- df_PIN1_clino %>%
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
df_binned <- df_norm_PINClino %>%
mutate(Bin = case_when(
  RelativePosition <= 0.25 ~ "1",
 RelativePosition <= 0.5 ~ "2",
RelativePosition <= 0.75 ~ "3",
TRUE ~ "4"
))

df_binned$Bin <- as.factor(df_binned$Bin)
df_binned[2,2] <- "control"

image_bin_averages_PINClino <- df_binned %>%
  group_by(ID, group, Bin) %>%
  summarise(AvgGrayValue = mean(intensity, na.rm = TRUE), .groups = "drop")

image_bin_averages_PINCLino <- df_binned %>%
  group_by(ID, group, Bin) %>%
  summarise(AvgGrayValue = mean(intensity, na.rm = TRUE), .groups = "drop")



image_bin_averages_PINCLino$group <- factor(image_bin_averages_PINCLino$group)
#image_bin_averages_PINCLino$Bin <- factor(image_bin_averages_PINCLino$Bin, levels = c("Left 4", "Left 3", "Left 2", "Left 1","Right 1", "Right 2", "Right 3", "Right 4"))

#graph
library(ggplot2)
ggplot(data=image_bin_averages_PINCLino,aes(x=Bin,y=AvgGrayValue, color=group))+
  geom_point(show.legend=TRUE)+
  xlab("Treatment")+
  ylab("Ratio of GFP")
image_bin_averages_PINCLino$factorBin  <- as.factor(image_bin_averages_PINCLino$Bin)
ggplot(data=image_bin_averages_PINCLino,aes(x=factorBin,y=AvgGrayValue, color=group))+
  geom_boxplot(show.legend=TRUE)+
  xlab("Bin")+
  ylab("Average Gray Value")

#assumptions
#making the model 
PINClinolm<-lm(AvgGrayValue ~ group + Bin, data =image_bin_averages_PINCLino)

# normality
shapiro.test(residuals(PINClinolm))

qqnorm(residuals(PINClinolm))
qqline(residuals(PINClinolm))

#var
boxplot(residuals(PINClinolm) ~ image_bin_averages_PINCLino$group)
boxplot(residuals(PINClinolm) ~ image_bin_averages_PINCLino$Bin)
boxplot(residuals(PINClinolm) ~ image_bin_averages_PINCLino$group+image_bin_averages_PINCLino$Bin)



anova_result_PINClino <- aov(AvgGrayValue ~ group+Bin, data = image_bin_averages_PINCLino)

summary(anova_result_PINClino)

Anova_TUK_PINClino <- TukeyHSD(aov(AvgGrayValue ~ group + Bin, data = image_bin_averages_PINCLino))


library(multcomp)
anova_result_PINClino <- aov(AvgGrayValue ~ group + Bin, data = image_bin_averages_PINCLino)
Anova_TUK_PINClino <- glht(anova_result_PINClino, linfct=mcp( Bin ="Tukey"))
cld(Anova_TUK_PINClino) 

summary(Anova_TUK_PINClino)


TukeyHSD(aov(AvgGrayValue ~ group + Bin, data = image_bin_averages_PINCLino))
plot(TukeyHSD(aov(AvgGrayValue ~ group * Bin, data = image_bin_averages_PINCLino)))
