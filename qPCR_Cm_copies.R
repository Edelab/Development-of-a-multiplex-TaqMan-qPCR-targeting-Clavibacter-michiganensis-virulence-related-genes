# Cleaning of the R environment
rm(list=ls())

# Set your working directory to where your excel file is saved:
setwd("C:/Users/Edelab/Documents/R_Statistique/Statistique_qPCR/Quantification/Code_Data_final")

# Import the data
library(openxlsx)

df <- within(read.xlsx(xlsxFile="qPCR_Cm_copies_tissues.xlsx"),{
  Plant = as.factor(Plant)
  Std = as.factor(Std)
  Fluo = as.factor(Fluo)
  N_echan = as.factor(N_echan)
  Tissu = as.factor(Tissu)
  Nb_copies= as.numeric(Nb_copies)
  Log_Nb_copies = as.numeric(Log_Nb_copies)
  })
str(df)

length(levels(df$Plant))

##############################################################
#Comparison plasmid-based and bacterial-based standard curves#
##############################################################
# Investigation of a transformation
library(MASS)

boxcox = boxcox(Nb_copies ~ Fluo*Std*Tissu , data=df,
                lambda = seq(-2, 2, length = 1000))
# Transforming the copies number by log, in order to satisfy the normality assumption

# Analyse
library(nlme); library(emmeans); library(multcomp); library(multcompView)
fit = lme(log(Nb_copies) ~ Fluo*Std*Tissu, random = ~ 1 | Plant, data=df)

joint_tests(fit)

# Tukey's multiple comparisons
a1 = emmeans(fit, ~  Fluo*Std, type="response"); a1
cld(a1, Letters=letters)

a2 = emmeans(fit, ~ Tissu | Std , type="response"); a2
cld(a2, Letters=letters)

# Normality
r = residuals(fit,type="normalized", level=0)
hist(r,freq=F)
xfit<-seq(min(r),max(r),length=40)
yfit<-dnorm(xfit, mean=mean(r), sd=sd(r))
lines(xfit, yfit,col="red",lwd=2) 
shapiro.test(r) 
e1071::kurtosis(r)  


# Equality of variances
plot(fitted(fit, level=0), r, pch=16, ylab="Normalized residuals", xlab="Predicted values")
abline(h=0, lty=2)

#########################################################
#Graph plasmid-based and bacterial-based standard curves#
#########################################################
#Graph plasmid-based standard curve (Black/Gray)
library(ggplot2); library(ggpubr)
df$Tissu <- with(df, factor(Tissu, levels= c("Roots","Inoculation_site","3_cm","6_cm","Leaf_Rachis"),
                                labels =c("Roots","Inoculation site","Stem - 3 cm", "Stem - 6 cm", "Leaf")))
mydata1 = df[df$Std %in% "plasmid",]

p <- ggboxplot(mydata1, x = "Tissu", y = "Log_Nb_copies",fill = "Fluo", palette = c("gray80","gray32"))

bxp1 <- ggpar(p,
              legend = "none",
              ylab = "Log Cm copies/g of tissues",
              xlab = FALSE,
              font.xtickslab = 12,
              font.ytickslab = 10,
              font.y= 12,
              ylim = c(9,13))
bxp1

#Graph DNA-based standard curve (Black/Gray)

mydata2 = df[df$Std %in% "DNA",]

p1 <- ggboxplot(mydata2, x = "Tissu", y = "Log_Nb_copies",fill = "Fluo", palette = c("gray80","gray32"))

bxp2 <- ggpar(p1,
              ylab = "Log Cm copies/g of tissues",
              xlab = FALSE,
              font.xtickslab = 12,
              font.ytickslab = 10,
              font.y= 12,
              ylim = c(9,13),
              legend = "none")
bxp2

#Graph DNA-based standard curve (Red/blue)
p2 <- ggboxplot(mydata1, x = "Tissu", y = "Log_Nb_copies",color = "Fluo", palette =c("orangered3", "steelblue"),
                add = "jitter", shape = "Fluo")

bxp3 <- ggpar(p2,
              xlab = FALSE,
              legend = "right",
              legend.title = "Gene",
              x.text.angle = 45,
              font.legend = 11,
              ylab = "Log Cm copies/g of tissues",
              font.xtickslab = 12,
              font.ytickslab = 10,
              font.y= c(12,"bold"),
              ylim = c(9,13))
bxp3


p3 <- ggboxplot(mydata2, x = "Tissu", y = "Log_Nb_copies",color = "Fluo", palette =c("orangered3", "steelblue"),
               add = "jitter", shape = "Fluo")

bxp4 <- ggpar(p3,
              legend = "right",
              legend.title = "Gene",
              font.legend = 11,
              ylab = "Log Cm copies/g of tissues",
              xlab = FALSE,
              x.text.angle = 45,
              font.xtickslab = 12,
              font.ytickslab = 10,
              font.y= c(12,"bold"),
              ylim = c(10,13))
bxp4

