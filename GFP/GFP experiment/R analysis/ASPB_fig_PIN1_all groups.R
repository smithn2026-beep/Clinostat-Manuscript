pin1_aspb <- read.csv("C:/Users/smith/OneDrive - Ohio University/literally everything/Thesis/GFP experiment/PIN1_still control/R stuff/All_PIN1_Emma_updated.csv")

library(ggplot2)
library(dplyr)

pin1_aspb$Bin <- as.factor(pin1_aspb$Bin)


pin1_aspb$Variables <- paste(pin1_aspb$Treatment, pin1_aspb$group, sep="_")



ggplot(pin1_aspb, aes(x = as.factor(Bin), y = AvgGrayValue,color = Variables,))+
  geom_boxplot()+
  xlab("Bin")
 


ggplot(pin1_aspb, aes(x = as.factor(Bin), y = AvgGrayValue,color = Variables))+
  geom_boxplot()+
  xlab("Bin")

bigmodel <- lm(AvgGrayValue ~ Bin *Variables, data =pin1_aspb)
summary(bigmodel)    

summary(bigmodel)
# normality----
shapiro.test(residuals(bigmodel))

qqnorm(residuals(bigmodel))
qqline(residuals(bigmodel))

#var
boxplot(residuals(bigmodel) ~ pin1_aspb$group)
boxplot(residuals(bigmodel) ~ pin1_aspb$Bin)
boxplot(residuals(bigmodel) ~ pin1_aspb$group+pin1_aspb$Bin)



anova_bigModel <- aov(AvgGrayValue ~ Bin *Variables, data =pin1_aspb)

summary(anova_bigModel)

TukeyHSD(aov(AvgGrayValue ~ Bin *Variables, data =pin1_aspb))


library(multcompView)
anova_bigModel <- aov(AvgGrayValue ~ Bin *Variables, data =pin1_aspb)
Anova_TUK_bigModel <- glht(anova_bigModel, linfct=mcp( Bin ="Tukey"))
cld(Anova_TUK_bigModel) 

summary(anova_bigModel)

TukeyHSD(aov(AvgGrayValue ~ Bin *Variables, data =pin1_aspb))


library(writexl)

# Run the ANOVA model and save it to a variable
aov_model_pin <- aov(AvgGrayValue ~ Bin * Variables, data = pin1_aspb)

# Run Tukey HSD post-hoc test on the ANOVA model
tukey_results_pin <- TukeyHSD(aov_model_pin)

# Convert each comparison table to a dataframe and combine them into one
# TukeyHSD returns a separate table for each term (Bin, Variables, Bin:Variables)
# so we loop through each one and stack them into a single dataframe
tukey_df_pin <- do.call(rbind, lapply(names(tukey_results_pin), function(term) {
  
  # Convert the current term's results to a dataframe
  df_pin <- as.data.frame(tukey_results_pin[[term]])
  
  # Add a column to label which term (Bin, Variables, etc.) each row belongs to
  df_pin$Term <- term
  
  # Move the row names (the group comparisons) into their own column
  df_pin$Comparison <- rownames(df_pin)
  
  # Remove row names since they are now a proper column
  rownames(df_pin) <- NULL
  
  # Reorder columns so Term and Comparison appear first
  df_pin[, c("Term", "Comparison", "diff", "lwr", "upr", "p adj")]
}))

# Rename columns to be more descriptive and readable
colnames(tukey_df_pin) <- c("Term", "Comparison", "Difference", "Lower_CI", "Upper_CI", "P_Adjusted")

# Export the dataframe to an Excel file in your working directory
write_xlsx(tukey_df_pin, "tukey_pin_results.xlsx")

