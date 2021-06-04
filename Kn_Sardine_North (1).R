######################
# Charge libraries
######################
library(car)
library(nlme)
library(dismo)
library(Hmisc)
library(mgcv)
library(visreg)
library(gratia)
require(scales)
library(mgcViz)
library(visreg)
library(anytime)
library(ggplot2)
library(corrplot)
library(forecast)
library(extrafont)
library(cairoDevice)

############################################
# Read data from csv file or Rdata        #
############################################
data<- read.csv("Spil_GSA6N_Kn_decomposed_all.csv")

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
#                 GAM                     #
########################################### 

# Model with all environmental variables
m<- gam((Kn_trend)~s(CRTY_trend,k=4)+s(CRTX_trend,k=4)+s(SST_trend,k=4)+s(NPP_trend,k=4),data=data_Kn,method="REML")#
summary(m)
m$aic###### AIC=-696.8502; R-sq (adj) = 0.429;  Deviance explained= 46.6 %
appraise(m)
res <- resid(m)
plot(res)


# Variable CRTX was not significant and was removed 
m0<- gam((Kn_trend)~s(CRTY_trend,k=4)+s(SST_trend,k=4)+s(NPP_trend,k=4),data=data_Kn,method="REML")#
summary(m0)
m0$aic###### AIC=-691.2924; R-sq (adj) = 0.397;  Deviance explained= 42.6%
appraise(m0)
res <- resid(m0)
plot(res)

# check autocorrelation #
layout(matrix(1:2, ncol = 2))
acf(res, lag.max = 36, main = "ACF - AR(m) errors")
pacf(res, lag.max = 36, main = "pACF- AR(m) errors")
layout(1)

# High autocorrelation--> include a temporal variable to account for time

data_Kn$M=c(1:139) ### create a continuous variable.

m1<- gam((Kn_trend)~s(CRTY_trend,k=4)+s(SST_trend,k=4)+s(NPP_trend,k=4)+s(M,k=4),data=data_Kn,method="REML")#
summary(m1)
m1$aic###### AIC=-765.8567; R-sq (adj) = 0.653;  Deviance explained= 67.6%
appraise(m1)
res <- resid(m1)
plot(res)

# check autocorrelation #
layout(matrix(1:2, ncol = 2))
acf(res, lag.max = 36, main = "ACF - AR(m) errors")
pacf(res, lag.max = 36, main = "pACF- AR(m) errors")
layout(1)

# Autocorrelation was improved, but not eliminated. The model diagnosis plots improved as well. 
# Include correlation structure for the residuals in the GAM--> GAMM 
# First determine the values of the corARMA  
m1b<- gamm((Kn_trend)~s(CRTY_trend,k=4)+s(SST_trend,k=4)+s(NPP_trend,k=4)+s(M,k=4),data=data_Kn,method="REML")#
summary(m1b$gam)
summary(m1b$lme)
appraise(m1b$gam)

res <- resid(m1b$gam)
plot(res)

#check best ARIMA model for our data
arma_res <- auto.arima(resid(m1b$gam),
                       stationary = FALSE,trace=TRUE, 
                       seasonal = FALSE,approximation = FALSE,stepwise = TRUE)  # Best model: ARIMA(2,0,3) 


Arima(res,order=c(3,0,0)) ### calculate the coefficients 

##################  GAMM ################

m2<- gamm((Kn_trend)~s(CRTY_trend,k=4)+s(SST_trend,k=4)+s(NPP_trend,k=4)+s(M,k=4),
          data=data_Kn, correlation = corARMA(value=c(0.99,-0.8562,0.1830), p=3, q=0, fixed=F),method="REML")  #

#Plot diagnostic
par(mfrow=c(2,2))
acf(resid(m2$lme,type = "normalized"), main="ACF")
pacf(resid(m2$lme,type = "normalized"), main="PACF")
plot(LjungBoxTest(resid(m2$gam))[,3],ylim=c(0,1),col="black",ylab="P-valor",xlab="Retardo")
abline(h=0.05,lty=2,col="blue")
res<-resid(m2$lme,type = "normalized")
resp<-predict(m2$gam)
plot(res~resp,type="p",xlab="Predichos",ylab="Residuos")

#Compare models
anova(m1b$lme, m2$lme)  

# the autocorrelation in the residuals improved but the model diagnosis was not good and the R-sq (adj.) = -3.79

#################################
###           PLOTS           ###
#################################
#### We keep the Model m1 #

# plot partial effects


layout(matrix(1:2, ncol = 2))
b<- getViz(m1)
print(plot(b, allTerms =T), pages = 1) # Including also parametric effect

# to represent the partial effects rescaled we using the library "visreg"
#m2$gam$data<-data_Kn ###we need to specify the data when plotting partial effects of GAMM

layout(matrix(1:4, ncol = 4))
par(mfrow=c(1,4))

visreg (m1,"CRTY_trend",
        xlab="MC trend", ylab="Kn trend",legend=F,
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

