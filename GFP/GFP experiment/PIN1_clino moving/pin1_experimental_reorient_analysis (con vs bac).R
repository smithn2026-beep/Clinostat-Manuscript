
library(dplyr)
library(car)


# PIN1 reoriented analysis (control vs bacteria) from moving clinostat
df_PIN1_clino<-read.csv("C:/Users/smith/OneDrive - Ohio University/literally everything/Thesis/GFP experiment/PIN1_clino moving/PIN1_experimental_R_data.csv")
df_PIN1_clino$ID <- as.factor(df_PIN1_clino$ID)
df_PIN1_clino$group <- as.factor(df_PIN1_clino$group)

df_norm_PINclino <- df_PIN1_clino %>%
  group_by(ID) %>%
  mutate(RelativePosition = (distance - min(distance)) / (max(distance) - min(distance))) %>%
  ungroup()

# df_binned_PINclino <- df_norm_PINclino %>%
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
df_binned <- df_norm_PINclino %>%
mutate(Bin = case_when(
  RelativePosition <= 0.25 ~ "1",
 RelativePosition <= 0.5 ~ "2",
RelativePosition <= 0.75 ~ "3",
TRUE ~ "4"
))

#df_binned$Bin <- as.factor(df_binned$Bin)

#df_binned[2,2] <- "control"

image_bin_averages_PINclino <- df_binned %>%
  group_by(ID, group, Bin) %>%
  summarise(AvgGrayValue = mean(intensity, na.rm = TRUE), .groups = "drop")

image_bin_averages_PINclino <- df_binned %>%
  group_by(ID, group, Bin) %>%
  summarise(AvgGrayValue = mean(intensity, na.rm = TRUE), .groups = "drop")



image_bin_averages_PINclino$group <- factor(image_bin_averages_PINclino$group)
#image_bin_averages_PINclino$Bin <- factor(image_bin_averages_PINclino$Bin, levels = c("Left 4", "Left 3", "Left 2", "Left 1","Right 1", "Right 2", "Right 3", "Right 4"))

#graph
library(ggplot2)
ggplot(data=image_bin_averages_PINclino,aes(x=Bin,y=AvgGrayValue, color=group))+
  geom_point(show.legend=TRUE)+
  xlab("Treatment")+
  ylab("Ratio of GFP")
image_bin_averages_PINclino$factorBin  <- as.factor(image_bin_averages_PINclino$Bin)
ggplot(data=image_bin_averages_PINclino,aes(x=factorBin,y=AvgGrayValue, color=group))+
  geom_boxplot(show.legend=TRUE)+
  xlab("Bin")+
  ylab("Average Gray Value")

#assumptions
#making the model 
PINclinolm<-lm(AvgGrayValue ~ group + Bin, data =image_bin_averages_PINclino)

# normality
shapiro.test(residuals(PINclinolm))

qqnorm(residuals(PINclinolm))
qqline(residuals(PINclinolm))

#var
boxplot(residuals(PINclinolm) ~ image_bin_averages_PINclino$group)
boxplot(residuals(PINclinolm) ~ image_bin_averages_PINclino$Bin)
boxplot(residuals(PINclinolm) ~ image_bin_averages_PINclino$group+image_bin_averages_PINclino$Bin)



anova_result_PINclino <- aov(AvgGrayValue ~ group+Bin, data = image_bin_averages_PINclino)

summary(anova_result_PINclino)

Anova_TUK_PINclino <- TukeyHSD(aov(AvgGrayValue ~ group + Bin, data = image_bin_averages_PINclino))


library(multcomp)
anova_result_PINclino <- aov(AvgGrayValue ~ group + Bin, data = image_bin_averages_PINclino)
Anova_TUK_PINclino <- glht(anova_result_PINclino, linfct=mcp( Bin ="Tukey"))
cld(Anova_TUK_PINclino) 

summary(Anova_TUK_PINclino)


TukeyHSD(aov(AvgGrayValue ~ group + Bin, data = image_bin_averages_PINclino))
plot(TukeyHSD(aov(AvgGrayValue ~ group * Bin, data = image_bin_averages_PINclino)))

