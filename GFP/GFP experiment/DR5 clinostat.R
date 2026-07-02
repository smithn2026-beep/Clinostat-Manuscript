#Old DR5 clino moving bac and con----
library(ggplot2)
library(dplyr)
library(ggsignif)
DR5_group <- DR5_clino%>%
  group_by(ID)

# graphing -- boxplot
ggplot(data=DR5_clino,aes(x=ID,y=GFP.ratio, fill=ID))+
  geom_boxplot(show.legend=TRUE)+theme_grey()+
  xlab("Treatment")+
  ylab("Ratio of GFP")+
  scale_fill_manual(values=c("#ff483f","#6F7378"))
  

# subsetting groups to test for normality and var 
Con_DR5 <- subset(DR5_clino, ID == "con")
Bac_DR5 <- subset(DR5_clino, ID == "bac")


# norm
qqnorm(Con_DR5$GFP.ratio)
qqline(Con_DR5$GFP.ratio)

shapiro.test(Con_DR5$GFP.ratio)


qqnorm(Bac_DR5$GFP.ratio)
qqline(Bac_DR5$GFP.ratio)

shapiro.test(Bac_DR5$GFP.ratio)


var.test(Con_DR5$GFP.ratio, Bac_DR5$GFP.ratio)

t.test(Con_DR5$GFP.ratio ~ Bac_DR5$GFP.ratio, var=TRUE)
t.test(GFP.ratio~ID, data=DR5_clino, var.equal= TRUE)

.
.
.
# PIN1 still bac and con treatments----



