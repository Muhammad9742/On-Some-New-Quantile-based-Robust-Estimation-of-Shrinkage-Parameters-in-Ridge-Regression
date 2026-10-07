rm(list=ls())
set.seed(1990)
library(MASS)
p=4                                      #No of explanatory variables
I=diag(p)
setwd("E:/PHD/Real data")      #To change working directory
data=read.csv("live stock pop and production1.csv",header=TRUE)
x1=data$Milk   # use $ to add variables
x2=data$goat
x3=data$cattle
x4=data$Bones
n=length(x1)
x=cbind(x1,x2,x3,x4)
x=scale(x,center = TRUE,scale = TRUE)   #Scaling x,,, standardizing x and y
y=data$Blood
y=(y-mean(y))/sd(y)
y[3]=y[3]+15
y[8]=y[8]-15
#y[8]=y[8]-50

#y=scale(x,center = TRUE,scale = TRUE)
eigen(cor(x))
ev=eigen(cor(x))$values
vec=eigen(cor(x))$vectors
which.max(ev)
CN=max(ev)/min(ev)        #Computing condition number
CI=sqrt(max(ev)/min(ev))
beta=vec[,which.max(ev)]
C=t(x)%*%(x)
D=eigen(C)$vectors
Z=x%*%D
lam=t(Z)%*%Z
lam=round(lam,4)
lamda=diag(lam)
alpha=t(D)%*%beta

betahat=solve(t(x)%*%x)%*%t(x)%*%y
yhat=x%*%betahat
sigmahat=(sum((y-yhat)^2))/(n-p)
alphahat=solve(lam)%*%t(Z)%*%y
#y[5]=y[5]+10*sigmahat
#y[15]=y[15]-5*sigmahat
#y[20]=y[20]+20*sigmahat
#y[35]=y[35]-750*sigmahat
#y[40]=y[40]+1000*sigmahat

#Estimation of ridge parameter
#OLS estimator
alphahat.ols=alphahat
alphahat.ols=c(alphahat.ols)

#Ridge Regression estimator of HK method (1970)

HK=sigmahat/(max(alphahat^2))
K1=HK
T1=solve(lam+I*K1)%*%lam
alphahat.RR=T1%*%alphahat.ols
alphahat.RR = c(alphahat.RR)

##################Modified KL Estimator (2023)


IKL=((3*ev*alphahat^2+sigmahat)^2+(4*alphahat^2*sigmahat*ev)-(3*ev*alphahat^2)+sigmahat)
#AIKL=abs(IKL)
A=sqrt(abs(IKL))/(2*alphahat^2)
B=((3*ev*alphahat^2)+sigmahat)
IKLE=(A-B)

MKL1=min(IKLE)
MKL2=max(IKLE)

#K2=MKL1
#alphahat.MKL1=solve(lam+I*K2)%*%(lam-I*K2)%*%solve(lam+I*K2)%*%lam%*%alphahat
#alphahat.MKL1=c(alphahat.MKL1)

#K3=MKL2
#alphahat.MKL2=solve(lam+I*K3)%*%(lam-I*K3)%*%solve(lam+I*K3)%*%lam%*%alphahat
#alphahat.MKL2=c(alphahat.MKL2)

# M-estimator denoted by alphahat.m of Huber 1981
r.fit=rlm(y~x-1,psi=psi.huber,k=1.345,maxit = 1000)
alphahat.m=as.vector(r.fit$coefficients)
alphahat.m=c(alphahat.m)
a.huber=r.fit$k

# To obtain A2
e=r.fit$residuals
s=r.fit$s
u=e/s
w.psi=psi.huber(u)
u.psi=u*w.psi
deriv.psi=psi.huber(u,deriv = 1)
num.a2=((s^2)*sum(u.psi^2))/(n-p)
denum.a2=(sum(deriv.psi)/n)^2
A2=num.a2/denum.a2

# Silvapulle (1991) Proposed M-estimator of HKB(1975) & LW(1976)
 
  K2=p/sum(alphahat.m^2/(sqrt((lamda*(alphahat.m^2)*A2)+(A2^2))+A2))

  alphahat.MKL1=solve(lam+I*K2)%*%(lam-I*K2)%*%alphahat.m
  alphahat.MKL1=c(alphahat.MKL1)
  K3=quantile(A2/(alphahat.m^2),probs = 0.95)
  alphahat.MKL2=solve(lam+I*K3)%*%(lam-I*K3)%*%alphahat.m
  alphahat.MKL2=c(alphahat.MKL2)
  

  MHKB=(p*A2)/(sum(alphahat.m^2))
  K4=MHKB
  T2=solve(lam+I*K4)%*%lam
  alphahat.MHKB=T2%*%alphahat.m
  alphahat.MHKB=c(alphahat.MHKB)
  
  MLW=(p*A2)/(sum(ev*(alphahat.m^2)))
  K5=MLW
  T3=solve(lam+I*K5)%*%lam
  alphahat.MLW=T3%*%alphahat.m
  alphahat.MLW=c(alphahat.MLW)
  
  # Majid Kibria-Lukman M-estimator (2022) M-estimator  
  
  KLME=A2/((2*alphahat.m^2)+(A2/ev))
  KLME1=sum(KLME)/p
  KLME2=(prod(KLME))^(1/p)
  KLME3=p/(sum(1/KLME))
  
  K6=KLME1
  alphahat.KLME1=solve(lam+I*K6)%*%(lam-I*K6)%*%alphahat.m
  alphahat.KLME1=c(alphahat.KLME1)
  
  K7=KLME2
  alphahat.KLME2=solve(lam+I*K7)%*%(lam-I*K7)%*%alphahat.m
  alphahat.KLME2=c(alphahat.KLME2)
  
  K8=KLME3
  alphahat.KLME3=solve(lam+I*K8)%*%(lam-I*K8)%*%alphahat.m
  alphahat.KLME3=c(alphahat.KLME3)
  
  # Danish Proposed Robust improve Kibria-Lukman Estimator(2023)
  
  MKL=((3*ev*alphahat.m^2)^2+4*alphahat.m^2*A2*ev-3*ev*alphahat.m^2+A2)
  AMKL=abs(MKL)                         # only absolute is taken 
  A1=sqrt(AMKL)/(2*alphahat.m^2)
  B1=((3*ev*alphahat.m^2)+A2)
  IKLME=abs(A1-B1)                     # only absolute is taken
  
  DSM1=median(IKLME)
  DSM2=(prod(IKLME))^(1/p)
  DSM3=p/(sum(1/IKLME))
       W=max(ev)*sd(ev)/sum(abs(alphahat.m))  
  MUS=W*B1 
  DSM4=quantile(MUS,prob=0.10)
  DSM5=quantile(MUS,prob=0.25)
  DSM6=quantile(MUS,prob=0.50)
  DSM7=quantile(MUS,prob=0.75)
  DSM8=quantile(MUS,prob=0.95)
  DSM9=quantile(MUS,prob=0.99)
  
  K9=DSM1
  alphahat.DSM1=solve(lam+I*K9)%*%(lam-I*K9)%*%solve(lam+I*K9)%*%lam%*%alphahat.m
  alphahat.DSM1=c(alphahat.DSM1)
  
  K10=DSM2
  alphahat.DSM2=solve(lam+I*K10)%*%(lam-I*K10)%*%solve(lam+I*K10)%*%lam%*%alphahat.m
  alphahat.DSM2=c(alphahat.DSM2)
  
  K11=DSM3
  alphahat.DSM3=solve(lam+I*K11)%*%(lam-I*K11)%*%solve(lam+I*K11)%*%lam%*%alphahat.m
  alphahat.DSM3=c(alphahat.DSM3)
  
  K12=DSM4
  alphahat.DSM4=solve(lam+I*K12)%*%(lam-I*K12)%*%solve(lam+I*K12)%*%lam%*%alphahat.m
  alphahat.DSM4=c(alphahat.DSM4)
  
  K13=DSM5
  alphahat.DSM5=solve(lam+I*K13)%*%(lam-I*K13)%*%solve(lam+I*K13)%*%lam%*%alphahat.m
  alphahat.DSM5=c(alphahat.DSM5)
  
  K14=DSM6
  alphahat.DSM6=solve(lam+I*K14)%*%(lam-I*K14)%*%solve(lam+I*K14)%*%lam%*%alphahat.m
  alphahat.DSM6=c(alphahat.DSM6)
  
  K15=DSM7
  alphahat.DSM7=solve(lam+I*K15)%*%(lam-I*K15)%*%solve(lam+I*K15)%*%lam%*%alphahat.m
  alphahat.DSM7=c(alphahat.DSM7)
  K16=DSM8
  alphahat.DSM8=solve(lam+I*K16)%*%(lam-I*K16)%*%solve(lam+I*K16)%*%lam%*%alphahat.m
  alphahat.DSM8=c(alphahat.DSM8)
  
  K17=DSM9
  alphahat.DSM9=solve(lam+I*K17)%*%(lam-I*K17)%*%solve(lam+I*K17)%*%lam%*%alphahat.m
  alphahat.DSM9=c(alphahat.DSM9)
  
  
comb.coeff=rbind(alphahat.ols,alphahat.RR,alphahat.m,alphahat.MKL1,alphahat.MKL2,alphahat.MHKB,
                 alphahat.MLW,alphahat.KLME1,alphahat.KLME2,alphahat.KLME3,alphahat.DSM1,alphahat.DSM2,
                 alphahat.DSM3,alphahat.DSM4,alphahat.DSM5,alphahat.DSM6,alphahat.DSM7,alphahat.DSM8,alphahat.DSM9)
comb.coeff=round(comb.coeff,6)

###########################################################################################################
#MSE of all estimators
#MSE of OLS and RR

K<-c(0,K1) 
a=rep(0,length(K)); b=rep(0,length(K))                            #to obtain MSE of ridge estimators
for(i in 1:length(K)){
  
  a[i]=sum(ev/(ev+K[i])^2)
  b[i]=sum(((K[i]^2)*(alphahat.ols^2))/(ev+K[i])^2)
  
}
MSE1=sigmahat*a+b
#MSE of M and Ridge M-estimators
Km=c(0)
a.m=rep(0,length(Km)); b.m=rep(0,length(Km))
for (i in 1:length(Km)) {
  a.m[i]=sum(((ev^2)*A2)/(ev*(ev+Km[i])^2))
  b.m[i]=sum(((Km[i]^2)*(alphahat.m^2))/(ev+Km[i])^2)
}
MSE2=a.m+b.m

#MSE of Kibria-Lukman Ridge M-estimators

Km1=c(K4,K5,K6)
a.m1=rep(1,length(Km1)); b.m1=rep(1,length(Km1))
for (i in 1:length(Km1)) {
  a.m1[i]=sum(((ev-Km1[i])^2*A2)/(ev*(ev+Km1[i])^2))
  b.m1[i]=sum(((Km1[i]^2)*(alphahat.m^2))/(ev+Km1[i])^2)
}
MSE3=a.m1+4*b.m1

#MSE of Robust improve Kibria-Lukman estimators (2023)

Km2=c(K2,K3,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,K17)
a.m2=rep(1,length(Km2)); b.m2=rep(1,length(Km2))
for (i in 1:length(Km2)) {
  a.m2[i]=sum((ev*(ev-Km2[i])^2*A2)/(ev+Km2[i])^4)
  b.m2[i]=sum((Km2[i]^2)*(3*ev+Km2[i])^2*alphahat.m^2)/((ev+Km2[i])^4)
}
MSE4=a.m2+b.m2

  #a.m1=sum(((ev^2)*A2)/(ev*(ev+KLME)^2))
  #b.m1=sum(((KLME^2)*(alphahat.m^2))/(ev+KLME)^2)
  #MSE3=a.m1+b.m1

  #a.m2=sum(((ev^2)*A2)/(ev*(ev+DSS11)^2))
  #b.m2=sum(((DSS11^2)*(alphahat.m^2))/(ev+DSS11)^2)
  
  # MSE4=a.m2+b.m2

MSE=c(MSE1,MSE2,MSE3,MSE4)
MSE=as.matrix(MSE)
MSE=round(MSE,8)
K.hat=c(K1,K2,K3,K4,K5,K6,K7,K8,K9,K10,K11,K12,K13,K14,K15,K16,K17)
write.csv(K.hat,file ="D:/PHD/Kvalues.csv")
row.names(MSE)=c("OLS","RR","HUB.M","S1","S2","MHKB","MLW","KLME1","KLME2","KLME3","DSM1","DSM2","DSM3","DSM4","DSM5","DSM6","DSM7","DSM8","DSM9")
colnames(MSE)=c("MSE")
write.csv(MSE,file ="D:/PHD/MSER.csv")
write.csv(comb.coeff,file ="D:/PHD/coefficientsR.csv")
############################################
# 95% Mean Confidence Interval at x0
# for all 19 estimators (SINGLE CODE)
############################################

# x0 (original scale)
x0=c(46700,78,50,950)

# standardize x0 (model fitted on standardized X)
x.mean <- apply(x, 2, mean)
x.sd   <- apply(x, 2, sd)
x0s <- (x0 - x.mean) / x.sd

z <- qnorm(0.975)
XtX <- t(x) %*% x

# function for mean CI
meanCI <- function(beta_hat, K, sigma2){
  M <- solve(XtX + K*I) %*% XtX %*% solve(XtX + K*I)
  yhat0 <- as.numeric(t(x0s) %*% beta_hat)
  var0  <- as.numeric(sigma2 * t(x0s) %*% M %*% x0s)
  lower <- yhat0 - z * sqrt(var0)
  upper <- yhat0 + z * sqrt(var0)
  width <- upper - lower
  c(Mean=yhat0, Lower95=lower, Upper95=upper, Width95=width)
}

# Mean CI for all 19 estimators
CI.mean <- rbind(
  OLS   = meanCI(alphahat.ols, 0,  sigmahat),
  RR    = meanCI(alphahat.RR,  K1, sigmahat),
  HUB.M = meanCI(alphahat.m,   0,  A2),

  S1    = meanCI(alphahat.MKL1, K2, A2),
  S2    = meanCI(alphahat.MKL2, K3, A2),

  MHKB  = meanCI(alphahat.MHKB, K4, A2),
  MLW   = meanCI(alphahat.MLW,  K5, A2),

  KLME1 = meanCI(alphahat.KLME1, K6, A2),
  KLME2 = meanCI(alphahat.KLME2, K7, A2),
  KLME3 = meanCI(alphahat.KLME3, K8, A2),

  DSM1  = meanCI(alphahat.DSM1, K9,  A2),
  DSM2  = meanCI(alphahat.DSM2, K10, A2),
  DSM3  = meanCI(alphahat.DSM3, K11, A2),
  DSM4  = meanCI(alphahat.DSM4, K12, A2),
  DSM5  = meanCI(alphahat.DSM5, K13, A2),
  DSM6  = meanCI(alphahat.DSM6, K14, A2),
  DSM7  = meanCI(alphahat.DSM7, K15, A2),
  DSM8  = meanCI(alphahat.DSM8, K16, A2),
  DSM9  = meanCI(alphahat.DSM9, K17, A2)
)

CI.mean <- round(CI.mean, 8)

# save result
write.csv(CI.mean, "D:/PHD/MeanCI_95_at_x0.csv")

