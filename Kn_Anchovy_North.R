# Limpieza
rm(list=ls()) 
ls()

######################
# Charge libraries
######################
library(car)
library(nlme)
library(dismo)
library(Hmisc)
library(anytime)
library(ggplot2)
library(corrplot)
library(visreg)
library(mgcv)
library(gratia)
require(scales)
library(mgcViz)
library(visreg)
library(cairoDevice)
library(extrafont)
library(FitAR)


############################################
# Read data from csv file or Rdata        #
############################################
setwd("C:/Users/Marta/Dropbox/TimeSeries_Deseasonalyze/Decomposed_datasets")

data<- read.csv("Eenc_GSA6N_Kn_decomposed_all.csv")


str(data)
date=anydate(data$Date);date
data$Year <- as.numeric(format(date,'%Y'))
data$Month<- as.numeric(format(date,'%d'))
data$day<- as.numeric(format(date,'%m'))


###########################################################
#           Check correlation Kn and  Collinearity        #
###########################################################
# Correlations #
matrix<-rcorr(as.matrix(data[,c(7,10,13,16)]), type = "pearson")

cor.mtest <- function(mat, ...) {
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n - 1)) {
    for (j in (i + 1):n) {
      tmp <- cor.test(mat[, i], mat[, j], ...)
      p.mat[i, j] <- p.mat[j, i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  p.mat
}

corrplot(matrix$r, type="lower", tl.col = "black",method="number",
         p.mat = matrix$P, sig.level = 0.05)

# Collinearity was tested by computing the generalized variance-inflation factors (GVIF)
source("HighstatLib.r") # this file should be in the same folder as the working space 
corvif(data[,c(7,10,13,16)]) #only environmental variables



###################
#### Plot data ####
###################
# check data
data_Kn<-data[!is.na(data$Kn_trend),]

plot(density(data_Kn$Kn_trend))
data_Kn<-transform(data_Kn, Kn_log= log(Kn_trend))

plot(density(data_Kn$Kn_log))

###########################################
###########  GAM  #########################
###########################################  

#### Model with all environmental variables
m<- gam((Kn_trend)~s(CRTY_trend,k=4)+s(CRTX_trend,k=4)+s(SST_trend,k=4)+s(NPP_trend,k=4),data=data_Kn,method="REML")#
summary(m)
m$aic###### AIC=--840.0928; R-sq (adj) = 0.556;  Deviance explained= 58.7 %
res <- resid(m)
plot(res)
appraise(m)


##### From CRTY and CRTX that are correlated we remove CRTY that explained less
m0<- gam((Kn_trend)~s(CRTX_trend,k=4)+s(SST_trend,k=4)+s(NPP_trend,k=4),data=data_Kn,method="REML")#
summary(m0)

m0$aic######AIC=-831.5554; R-sq (adj) = 0.523;  Deviance explained= 55.3 %
res <- resid(m0)
plot(res)
appraise(m0)

layout(matrix(1:2, ncol = 2))
acf(res, lag.max = 36, main = "ACF - AR(0) errors")
pacf(res, lag.max = 36, main = "pACF- AR(0) errors")
layout(1) 


##### From CRTY and CRTX that are correlated we remove CRTY that explained less
m0b<- gam((Kn_trend)~s(CRTX_trend,k=4)+s(NPP_trend,k=4),data=data_Kn,method="REML")#
summary(m0b)

m0b$aic######
res <- resid(m0b)
plot(res)
appraise(m0b)

layout(matrix(1:2, ncol = 2))
acf(res, lag.max = 36, main = "ACF - AR(0) errors")
pacf(res, lag.max = 36, main = "pACF- AR(0) errors")
layout(1) 


# High autocorrelation--> include a temporal variable to account for time


data_Kn$M=c(1:135) ### create a continuous variable.

m1<- gam((Kn_trend)~s(CRTX_trend,k=4)+s(SST_trend,k=4)+s(NPP_trend,k=4)+s(M,k=4),data=data_Kn,method="REML")#
summary(m1)
m1$aic###### AIC=-968.1185; R-sq (adj) = 0.83;  Deviance explained= 84.4%
appraise(m1)
res <- resid(m1)
plot(res)

# check autocorrelation #
layout(matrix(1:2, ncol = 2))
acf(res, lag.max = 36, main = "ACF - AR(m) errors")
pacf(res, lag.max = 36, main = "pACF- AR(m) errors")

layout(1)



# Model autocorrelation persist
# Include correlation structure for the residuals in the GAM--> GAMM 
# First determine the values of the corARMA 


m1a<- gamm((Kn_trend)~s(CRTX_trend,k=4)+s(SST_trend,k=4)+s(NPP_trend,k=4),data=data_Kn,method="REML")#
summary(m1a$gam)

m1b<-gamm((Kn_trend)~s(CRTX_trend,k=4)+s(SST_trend,k=4)+s(NPP_trend,k=4)+s(M,k=4),data=data_Kn,method="REML")
summary(m1b$gam)
summary(m1b$lme)
appraise(m1b$gam)

res <- resid(m1b$gam)
plot(res)

library(forecast)
arma_res <- auto.arima(resid(m1b$gam),
                       stationary = FALSE,trace=TRUE, 
                       seasonal = FALSE,approximation = FALSE,stepwise = TRUE)  # Best model: ARIMA(2,0,3) 


Arima(res,order=c(2,0,1)) ### calculate the coefficients 

##################  GAMM ################

m2<- gamm((Kn_trend)~s(CRTX_trend,k=4)+s(SST_trend,k=4)+s(NPP_trend,k=4)+s(M,k=4),
          data=data_Kn, correlation = corARMA(value=c(0.9011,-0.3397,0.4122), p=2, q=1, fixed=F),method="REML")  #

summary(m2$gam)
summary(m2$lme)

par(mfrow=c(2,2))
acf(resid(m2$lme,type = "normalized"), main="ACF")
pacf(resid(m2$lme,type = "normalized"), main="PACF")
plot(LjungBoxTest(resid(m2$gam))[,3],ylim=c(0,1),col="black",ylab="P-valor",xlab="Retardo")
abline(h=0.05,lty=2,col="blue")
res<-resid(m2$lme,type = "normalized")
resp<-predict(m2$gam)
plot(res~resp,type="p",xlab="Predicted",ylab="Residuals")


#Compare models
anova(m1a$lme, m1b$lme, m2$lme)  ###

#################################
###           PLOTS           ###
#################################
#### We keep the Model m1 #

# plot partial effects
layout(matrix(1:2, ncol = 2))
b<- getViz(m1)
print(plot(b, allTerms =T), pages = 1) # Including also parametric effect

# to represent the partial effects rescaled we using the library "visreg"
#  m2$gam$data<-data_Kn ###we need to specify the data when plotting partial effects of GAMM

layout(matrix(1:4, ncol = 4))
par(mfrow=c(1,4))


visreg (m1,"CRTX_trend",
        xlab="ZC trend", ylab="Kn trend",legend=F,
        rug=FALSE, partial=TRUE,
        fill=list(col="#481B6DFF"),
        points = list(pch=20), line = list( col="black"))+
  theme_classic()

visreg (m1,"SST_trend",
        xlab="ST150 trend", ylab="Kn trend",legend=F,
        rug=FALSE, partial=TRUE,
        fill=list(col="#481B6DFF"),
        points = list(pch=20), line = list( col="black"))+
  theme_classic()

visreg (m1,"NPP_trend",
        xlab="NPP trend", ylab="Kn trend",legend=F,
        rug=FALSE, partial=TRUE,
        fill=list(col="#481B6DFF"),
        points = list(pch=20), line = list( col="black"))+
  theme_classic()

visreg (m1,"M",
        xlab="Time", ylab="Kn trend",legend=F,
        rug=FALSE, partial=TRUE,
        fill=list(col="#481B6DFF"),
        points = list(pch=20), line = list( col="black"))+
  theme_classic()






