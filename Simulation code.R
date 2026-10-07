
# Authors: Muhammad Khan, Dr. Muhammad Suhail
###################################################################

rm(list=ls())
set.seed(7860)
library(MASS)

sigma1 = c(0.1,1,2,5)
ssize=c(25,50,100)
pred=c(4,10)
corr=c(0.85,0.95,0.99,0.999)
N=10000                                 #Simulation runs

MSE.r=matrix(0,19,length(corr))
MSE.s=array(0,dim = c(19,length(corr),length(sigma1)))
MSE.n=array(0,dim = c(19,length(corr),length(sigma1),length(ssize)))
MSE.p=array(0,dim = c(19,length(corr),length(sigma1),length(ssize),length(pred)))

#-----------------------
#Predictors loop
for(b in 1:length(pred)){
  p=pred[b]
  I=diag(p)
  
  #------------------------
  #Sample size loop
  for(a in 1:length(ssize)){
    n=ssize[a]
    
z=matrix(0,n,p)
for (i in 1:p) {
  z[,i]=rnorm(n,0,1)                     #Random no generation   
}
x=matrix(0,n,p)
#--------------------------
#Error variance loop
for(S in 1:length(sigma1)){
  sigma=sigma1[S]
for(j in 1:length(corr)){
  r=corr[j]
  
  for (i in 1:p) {
    x[,i]=sqrt(1-r^2)*z[,i]+r*z[,p]    #Generating x-matrix
  }
  
  x=scale(x,center = TRUE,scale = TRUE)    #standardizing x-matrix
  c=t(x)%*%x
  D=eigen(c)$vectors
  ev=eigen(c)$values
  beta=D[,which.max(ev)]
  Z=x%*%D
  lam=t(Z)%*%Z
  lam=round(lam,5)
  lamda=diag(lam)
  alpha=t(D)%*%beta
  
  #Defining Vectors
  K1=rep(0,N)
  K2=rep(0,N)
  K3=rep(0,N)
  K4=rep(0,N)
  K5=rep(0,N)
  K6=rep(0,N)
  K7=rep(0,N)
  K8=rep(0,N)
  K9=rep(0,N)
  K10=rep(0,N)
  K11=rep(0,N)
  K12=rep(0,N)
  K13=rep(0,N)
  K14=rep(0,N)
  K15=rep(0,N)
  K16=rep(0,N)
  K17=rep(0,N)
  #Defining matrices
  
  alphahat.ols=matrix(0,nrow = p,ncol = N)
  alphahatk.RR=matrix(0,nrow = p,ncol = N)
  alphahat.MKL1=matrix(0,nrow = p,ncol = N)
  alphahat.MKL2=matrix(0,nrow = p,ncol = N)
  alphahat.m=matrix(0,nrow = p,ncol = N)
  alphahat.MHKB=matrix(0,nrow = p,ncol = N)
  alphahat.MLW=matrix(0,nrow = p,ncol = N)
  alphahat.KLME1=matrix(0,nrow = p,ncol = N)
  alphahat.KLME2=matrix(0,nrow = p,ncol = N)
  alphahat.KLME3=matrix(0,nrow = p,ncol = N)
  alphahat.DSM1=matrix(0,nrow = p,ncol = N)
  alphahat.DSM2=matrix(0,nrow = p,ncol = N)
  alphahat.DSM3=matrix(0,nrow = p,ncol = N)
  alphahat.DSM4=matrix(0,nrow = p,ncol = N)
  alphahat.DSM5=matrix(0,nrow = p,ncol = N)
  alphahat.DSM6=matrix(0,nrow = p,ncol = N)
  alphahat.DSM7=matrix(0,nrow = p,ncol = N)
  alphahat.DSM8=matrix(0,nrow = p,ncol = N)
  alphahat.DSM9=matrix(0,nrow = p,ncol = N)
  
  #Simulation loop starts here
  
  for (i in 1:N) {
    e1=rnorm(ceiling(0.80*n),0,sigma)
    e2=rnorm(floor(0.20*n),0,sigma+10)
    e=c(e1,e2)
    y=Z%*%alpha+e
    betahat=solve(t(x)%*%x)%*%t(x)%*%y
    yhat=x%*%betahat
    sigmahat=sum((y-yhat)^2)/(n-p)
    alphahat=solve(lam)%*%t(Z)%*%y
    
    #OLS estimator
    
    alphahat.ols[,i]=alphahat
    
    #Ridge Regression estimator of HKB method (1975)
    
    HK=sigmahat/(max(alphahat^2))
    K1[i]=HK
    alphahatk=solve(t(Z)%*%(Z)+I*K1[i])%*%t(Z)%*%y
    alphahatk=c(alphahatk)
    alphahatk.RR[,i] = alphahatk
    
    
    ##################Modified KL Estimator (2023)
    
    
    IKL=((3*ev*alphahat.ols[,i]^2+sigmahat)^2+(4*alphahat.ols[,i]^2*sigmahat*ev)-(3*ev*alphahat.ols[,i]^2)+sigmahat)
    #AIKL=abs(IKL)
    A=sqrt(abs(IKL))/(2*alphahat.ols[,i]^2)
    B=((3*ev*alphahat.ols[,i]^2)+sigmahat)
    IKLE=(A-B)
    
    #MKL1=min(IKLE)
    #MKL2=max(IKLE)
    
    #K2[i]=MKL1
    #alphahat.MKL1[,i]=solve(lam+I*K2[i])%*%(lam-I*K2[i])%*%solve(lam+I*K2[i])%*%lam%*%alphahat.ols[,i]
    
    #K3[i]=MKL2
    #alphahat.MKL2[,i]=solve(lam+I*K3[i])%*%(lam-I*K3[i])%*%solve(lam+I*K3[i])%*%lam%*%alphahat.ols[,i]
    
    # M-estimator denoted by alphahat.m of Huber 1981
    
    r.fit=rlm(y~x-1,psi=psi.huber,k=1.345,maxit = 1000)
    alphahat.m[,i]=as.vector(r.fit$coefficients)
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
       
    QS=p/sum(alphahat.m[,i]^2/(sqrt((lamda*(alphahat.m[,i]^2)*A2)+(A2^2))+A2))
    K2[i]=QS
    alphahat.MKL1[,i]=solve(lam+I*K2[i])%*%lam%*%alphahat.m[,i]
    ######
    QS1=quantile(A2/(alphahat.m[,i]^2),probs = 0.95)

    K3[i]=QS1
    alphahat.MKL2[,i]=solve(lam+I*K3[i])%*%lam%*%alphahat.m[,i]
    
    # Silvapulle (1991) Proposed M-estimator of HKB(1975) & LW(1976) 
    
    MHKB=(p*A2)/(sum(alphahat.m[,i]^2))
    K5[i]=MHKB
    alphahat.MHKB[,i]=solve(lam+I*K5[i])%*%lam%*%alphahat.m[,i]
    
    MLW=(p*A2)/(sum(ev*(alphahat.m[,i]^2)))
    K6[i]=MLW
    alphahat.MLW[,i]=solve(lam+I*K6[i])%*%lam%*%alphahat.m[,i]      
    
    # Majid Kibria-Lukman M-estimator (2022) M-estimator A.M, GM & HM  
    
    KLME=A2/((2*alphahat.m[,i]^2)+(A2/ev))
    
    KLME1=sum(KLME)/p
    KLME2=(prod(KLME))^(1/p)
    KLME3=p/(sum(1/KLME))
    
    K6[i]=KLME1
    alphahat.KLME1[,i]=solve(lam+I*K6[i])%*%(lam-I*K4[i])%*%alphahat.m[,i]
    
    K7[i]=KLME2
    alphahat.KLME2[,i]=solve(lam+I*K7[i])%*%(lam-I*K5[i])%*%alphahat.m[,i]
    
    K8[i]=KLME3
    alphahat.KLME3[,i]=solve(lam+I*K8[i])%*%(lam-I*K8[i])%*%alphahat.m[,i]
    
    
    # Danish Proposed Robust Improve Kibria-Lukman Estimator(2023)
    
    # MKL=(9*ev^2*alphahat.m[,i]^4+10*alphahat.m[,i]^2*A2*ev+A2^2)
    
    MKL=(((3*ev*alphahat.m[,i]^2)+A2)^2+(4*alphahat.m[,i]^2*A2*ev)-(3*ev*alphahat.m[,i]^2+A2))
                      
    A1=sqrt(abs(MKL))/(2*alphahat.m[,i]^2) # only absolute is taken 
    B1=((3*ev*alphahat.m[,i]^2)+A2)
    IKLME=abs(A1-B1) 
     W=max(ev)*sd(ev)/sum(abs(alphahat.m[,i]))  

    MUS=W*B1                   # only absolute is taken
    
    DSM1=median(IKLME)
    DSM2=(prod(IKLME))^(1/p)
    DSM3=p/(sum(1/IKLME))
 #Proposed estimator by M khan
    DSM4=quantile(MUS,prob=0.10)
    DSM5=quantile(MUS,prob=0.25)
    DSM6=quantile(MUS,prob=0.50)
    DSM7=quantile(MUS,prob=0.75)
    DSM8=quantile(MUS,prob=0.95)
    DSM9=quantile(MUS,prob=0.99)
    
    
    K9[i]=DSM1
    alphahat.DSM1[,i]=solve(lam+I*K9[i])%*%(lam-I*K9[i])%*%solve(lam+I*K9[i])%*%lam%*%alphahat.m[,i]
    
    K10[i]=DSM2
    alphahat.DSM2[,i]=solve(lam+I*K10[i])%*%(lam-I*K10[i])%*%solve(lam+I*K10[i])%*%lam%*%alphahat.m[,i]
    
    K11[i]=DSM3
    alphahat.DSM3[,i]=solve(lam+I*K11[i])%*%(lam-I*K11[i])%*%solve(lam+I*K11[i])%*%lam%*%alphahat.m[,i]
    
    K12[i]=DSM4
    alphahat.DSM4[,i]=solve(lam+I*K12[i])%*%(lam-I*K12[i])%*%solve(lam+I*K12[i])%*%lam%*%alphahat.m[,i]
    
    K13[i]=DSM5
    alphahat.DSM5[,i]=solve(lam+I*K13[i])%*%(lam-I*K13[i])%*%solve(lam+I*K13[i])%*%lam%*%alphahat.m[,i]
    
    K14[i]=DSM6
    alphahat.DSM6[,i]=solve(lam+I*K14[i])%*%(lam-I*K14[i])%*%solve(lam+I*K14[i])%*%lam%*%alphahat.m[,i]
    
    K15[i]=DSM7
    alphahat.DSM7[,i]=solve(lam+I*K15[i])%*%(lam-I*K15[i])%*%solve(lam+I*K15[i])%*%lam%*%alphahat.m[,i]
    K16[i]=DSM8
    alphahat.DSM8[,i]=solve(lam+I*K16[i])%*%(lam-I*K16[i])%*%solve(lam+I*K16[i])%*%lam%*%alphahat.m[,i]
    
    K17[i]=DSM9 
    alphahat.DSM9[,i]=solve(lam+I*K17[i])%*%(lam-I*K17[i])%*%solve(lam+I*K17[i])%*%lam%*%alphahat.m[,i]
  } # End of simulation loop for N=5000 runs.
  
  MSE.ols=sum((alphahat.ols-c(alpha))^2)/N
  MSE.RR=sum((alphahatk.RR-c(alpha))^2)/N
  MSE.MKL1=sum((alphahat.MKL1-c(alpha))^2)/N
  MSE.MKL2=sum((alphahat.MKL2-c(alpha))^2)/N
  MSE.m=sum((alphahat.m-c(alpha))^2)/N
  MSE.MHKB=sum((alphahat.MHKB-c(alpha))^2)/N
  MSE.MLW=sum((alphahat.MLW-c(alpha))^2)/N
  MSE.KLME1=sum((alphahat.KLME1-c(alpha))^2)/N
  MSE.KLME2=sum((alphahat.KLME2-c(alpha))^2)/N
  MSE.KLME3=sum((alphahat.KLME3-c(alpha))^2)/N
  MSE.DSM1=sum((alphahat.DSM1-c(alpha))^2)/N
  MSE.DSM2=sum((alphahat.DSM2-c(alpha))^2)/N
  MSE.DSM3=sum((alphahat.DSM3-c(alpha))^2)/N
  MSE.DSM4=sum((alphahat.DSM4-c(alpha))^2)/N
  MSE.DSM5=sum((alphahat.DSM5-c(alpha))^2)/N
  MSE.DSM6=sum((alphahat.DSM6-c(alpha))^2)/N
  MSE.DSM7=sum((alphahat.DSM7-c(alpha))^2)/N
  MSE.DSM8=sum((alphahat.DSM8-c(alpha))^2)/N
  MSE.DSM9=sum((alphahat.DSM9-c(alpha))^2)/N
  
  
  MSE=c(MSE.ols,MSE.RR,MSE.m,MSE.MHKB,MSE.MLW,MSE.MKL1,MSE.MKL2,MSE.KLME1,MSE.KLME2,MSE.KLME3,MSE.DSM1,MSE.DSM2,
        MSE.DSM3,MSE.DSM4,MSE.DSM5,MSE.DSM6,MSE.DSM7,MSE.DSM8,MSE.DSM9)
  MSE=round(MSE,7)
  as.matrix(MSE)
  MSE.r[,j]=MSE
  
} #End of correlation loop
  MSE.s[,,S]=MSE.r
} #End of error variance loop
MSE.n[,,,a]=MSE.s
  }  #End of sample size loop
  MSE.p[,,,,b]=MSE.n
}   #End of predictors loop


col.names=c("0.85","0.95","0.99","0.999")
row.names=c("OLS","RR","HUB.M","MHKB","MLW","MKL1","MKL2","KLME1","KLME2","KLME3","DSM1","DSM2","DSM3","DSM4",
            "DSM5","DSM6","DSM7","DSM8","DSM9")
matrix.names=c("0.1","1","2","5")
matrix.names2=c("25","50","100")
matrix.names3=c("4","10")
dimnames(MSE.p)=list(row.names,col.names,matrix.names,matrix.names2,matrix.names3)
write.csv(MSE.p,file = "D:/PHD/MSEQNORMAL.csv")


