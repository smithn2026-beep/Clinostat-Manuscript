
library(dplyr)
library(car)


# PIN1 reoriented analysis (bacteria vs bacteria on ground and in clino) 
df_PINbac_clino_ground<-read.csv("C:/Users/smith/OneDrive - Ohio University/literally everything/Thesis/GFP experiment/R analysis/pin bac vs bac.csv")
df_PINbac_clino_ground$ID <- as.factor(df_PINbac_clino_ground$ID)
df_PINbac_clino_ground$group <- as.factor(df_PINbac_clino_ground$Treatment)

df_norm_PINbac_clino_ground <- df_PINbac_clino_ground %>%
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
df_binned_gbac_v_clibac <- df_norm_PINbac_clino_ground %>%
  mutate(Bin = case_when(
    RelativePosition <= 0.25 ~ "1",
    RelativePosition <= 0.5 ~ "2",
    RelativePosition <= 0.75 ~ "3",
    TRUE ~ "4"))

#df_binned$Bin <- as.factor(df_binned$Bin)

#df_binned[2,2] <- "control"

image_bin_averages_PINbac <- df_binned_gbac_v_clibac %>%
  group_by(ID, Treatment, Bin) %>%
  summarise(AvgGrayValue = mean(intensity, na.rm = TRUE), .groups = "drop")

#image_bin_averages_PINclino <- df_binned %>%
# group_by(ID, group, Bin) %>%
# summarise(AvgGrayValue = mean(intensity, na.rm = TRUE), .groups = "drop")



##*****image_bin_averages_PINcontrols$Treatment <- factor(image_bin_averages_PINclino$Treatment)
#image_bin_averages_PINclino$Bin <- factor(image_bin_averages_PINclino$Bin, levels = c("Left 4", "Left 3", "Left 2", "Left 1","Right 1", "Right 2", "Right 3", "Right 4"))

#graph
library(ggplot2)
ggplot(data=image_bin_averages_PINbac,aes(x=Bin,y=AvgGrayValue, color=Treatment))+
  geom_point(show.legend=TRUE)+
  xlab("Treatment")+
  ylab("Ratio of GFP")
image_bin_averages_PINbac$factorBin  <- as.factor(image_bin_averages_PINbac$Bin)
ggplot(data=image_bin_averages_PINbac,aes(x=factorBin,y=AvgGrayValue, color=Treatment))+
  geom_boxplot(show.legend=TRUE)+
  xlab("Bin")+
  ylab("Average Gray Value")

#assumptions
#making the model 
PINbaclm<-lm(AvgGrayValue ~ Treatment + Bin, data =image_bin_averages_PINbac)

# normality
shapiro.test(residuals(PINbaclm))

qqnorm(residuals(PINbaclm))
qqline(residuals(PINbaclm))

#var
boxplot(residuals(PINbaclm) ~ image_bin_averages_PINbac$Treatment)
boxplot(residuals(PINbaclm) ~ image_bin_averages_PINbac$Bin)
boxplot(residuals(PINbaclm) ~ image_bin_averages_PINbac$Treatment+image_bin_averages_PINbac$Bin)


anova_result_PINbac <- aov(AvgGrayValue ~ Treatment+Bin, data = image_bin_averages_PINbac)

summary(anova_result_PINbac)






