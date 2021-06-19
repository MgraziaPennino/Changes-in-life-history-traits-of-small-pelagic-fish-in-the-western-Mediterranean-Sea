# Limpieza
rm(list=ls()) 
ls()
######################
# Charge libraries
######################
library(car)
library(nlme)
library(dismo)
library(anytime)
library(ggplot2)
library(corrplot)
library(Hmisc)
library(mgcViz)
library(mgcv)
library(visreg)
library(ggplot2)
library(cairoDevice)
library(extrafont)
library(gratia)
library(FitAR)


############################################
# Read data from csv file or Rdata        #
############################################

setwd("C:/Users/Marta/Dropbox/Albo et al. spf/Age_Length_GAM")
#data<- read.csv("sardina_Age_pruebaMonth.csv")

data<- read.csv("sardina_Age.csv")
str(data)
data$date<-as.Date(with(data, paste(Year, Month, Day, sep="-")), format= "%Y-%m-%d")

data_GSA1=subset(data,GSA=="GSA01")
#data_GSA6S=subset(data,GSA=="GSA06" & Area=="2")
#data_GSA6N=subset(data,GSA=="GSA06" & Area=="1")

##### NORTH #####
data_GSA1<-data_GSA1[with(data_GSA1, order(date)),]
data_GSA1<-transform(data_GSA1, Time= as.numeric(date)/10)

plot(density(data_GSA1$TALLA_CM))

data_GSA1$Age=as.factor(data_GSA1$Age)
#data_GSA6N$Mes=as.factor(data_GSA6N$Mes)


#########################################################
############    GAM          ############################
#########################################################
############
#### M1: Only Age ####
###########

m1<-gam(TALLA_CM~Age, data= data_GSA1,method="REML") 

summary(m1)
m1$aic###### 
res <- resid(m1)
plot(res)
appraise(m1)

# check autocorrelation #
layout(matrix(1:2, ncol = 2))
acf(resid(m1), lag.max = 36, main = "ACF - AR(m) errors")
pacf(resid(m1), lag.max = 36, main = "pACF- AR(m) errors")
layout(1) 

####################
##  M2: Include Year ###
###################
m2<-gam(TALLA_CM~Age+s(Year, k=4), data= data_GSA1,method="REML") 


summary(m2)
m2$aic###### 
res <- resid(m2)
plot(res)
appraise(m2)


# check autocorrelation #
layout(matrix(1:2, ncol = 2))
acf(resid(m2), lag.max = 36, main = "ACF - AR(m) errors")
pacf(resid(m2), lag.max = 36, main = "pACF- AR(m) errors")
layout(1) 


############################
## M3:  Include Year: by Age ###
############################

m3<-gam(TALLA_CM~Age+s(Year, by= Age, bs = 'fs', k=4), data= data_GSA1,method="REML") 

summary(m3)
m3$aic######
res <- resid(m3)
plot(res)
appraise(m3)

# check autocorrelation #
layout(matrix(1:2, ncol = 2))
acf(resid(m3), lag.max = 36, main = "ACF - AR(m) errors")
pacf(resid(m3), lag.max = 36, main = "pACF- AR(m) errors")
layout(1) 


visreg_data<-visreg(m3,"Year","Age", plot=FALSE)


plot(visreg_data, overlay= TRUE, xlab="Year", ylab="Total length (cm)",
     rug=FALSE, partial=FALSE, ylim=c(11,22),
     fill=list(alpha=0.5))+
  theme_classic()

########################################
##  M4: Try removing Age as single factor ###
########################################

m4<-gam(TALLA_CM~s(Year, by= Age, bs = 'fs', k=4), data= data_GSA1,method="REML") 


summary(m4)
m4$aic###### 
res <- resid(m4)
plot(res)
appraise(m4)


# check autocorrelation #
layout(matrix(1:2, ncol = 2))
acf(resid(m4), lag.max = 36, main = "ACF - AR(m) errors")
pacf(resid(m4), lag.max = 36, main = "pACF- AR(m) errors")
layout(1) 



visreg_data<-visreg(m4,"Year","Age", plot=FALSE)

plot(visreg_data, overlay= TRUE, xlab="Year", ylab="Total length (cm)",
     rug=FALSE, partial=FALSE, ylim=c(11,30),
     fill=list(alpha=0.5))+
  theme_classic()



