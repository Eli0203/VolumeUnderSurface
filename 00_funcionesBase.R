###############################################
####ALGORITMO-VUS 3 CLASES
##ELIANA VALLEJO
###LOAD LIBRARIES
library(pracma)
###END LIBRARIES
#####Necessary functions
##functions
###box-cox transformation
transformBoxCox<-function(ldata){
  x<-ldata[[1]];y<-ldata[[2]];z<-ldata[[3]]
  minimum<-min(min(x),min(y),min(z))
 
  ###re-escaling dataset
  data<-lapply(ldata,function(i){
    if(minimum<0){
      i<-(i-minimum)+1
    }else{i=i}
  })
  
  ###################Likelihood function for lambda estimation###################
  likbox=function(h,data){
    x<-data[[1]];y<-data[[2]];z<-data[[3]]
     n<-length(x);m<-length(y);l<-length(z)
    if (h==0){
      xh=log(x)
      yh=log(y)
      zh=log(z)
    } else {
      xh=((x^h)-1)/h
      yh=((y^h)-1)/h
      zh=((z^h)-1)/h
    }
    oout=-n/2*log(sum((xh-sum(xh)/n)^2)/n)-m/2*log(sum((yh-sum(yh)/m)^2)/m)-
      l/2*log(sum((zh-sum(zh)/l)^2)/l)+(h-1)*(sum(log(x))+sum(log(y))+sum(log(z)))
    return(-oout)
  }
  h_ini<-0
  #Optimizing the Log-Likelihood
  hhat=optim(h_ini,likbox,data=data,method="BFGS")$par
  result<-lapply(data,function(i){
    if (hhat==0){
     i<-ifelse(!is.na(i),log(i),i)}
    else{
     i<-ifelse(!is.na(i),(i^hhat-1)/hhat,i)
    }

  })
  return(result)
}

#############################################################################
############EMPIRICAL METHOD################################################
####estimate function for empirical method####
F_emp<-function(fd,x){ecdf(fd)(x)}
q_emp<-function(fd,p){as.numeric(quantile(fd,p,na.rm=T))}
####INTEGRATION LIMITS####
lim_int<-function(p1){1-F_emp(f3x,F_emp(f1x,q_emp(f1x,p1)))}
####Empirical ROC Surface
roc_emp<-function(p3,p1){ F_emp(f2x,q_emp(f3x,1-p3))-F_emp(f2x,q_emp(f1x,p1))}
#############################################################################
############NORMAL METHOD################################################
####APROX VUS FOR NORMAL DISTRIBUTION####
lim_int_n<-function(p1){1-pnorm(qnorm(p1,mean(f1x),sd(f1x)),mean(f3x,na.rm=T),sd(f3x,na.rm=T))}
####Normal ROC Surface
roc_normal_est<-function(p3,p1){pnorm(qnorm(1-p3,mean(f3x,na.rm=T),sd(f3x,na.rm=T)),mean(f2x,na.rm=T),sd(f2x,na.rm=T))-
    pnorm(qnorm(p1,mean(f1x,na.rm=T),sd(f1x,na.rm=T)),mean(f2x,na.rm=T),sd(f2x,na.rm=T))}


##Article option
f0 <- function(s,a,b,c,d,p,q)
{
  ###integrate for s over (-Inf, Inf) given a,b,c,d to obtain VUS
  (pnorm(a*s-b)*pnorm(-c*s+d)-p*pnorm(-c*s+d)-q*pnorm(a*s-b)+p*q)*dnorm(s)     
}

#################################################################################
#####VUS ESTIMATION####
vusVar <-function(f1x,f2x,f3x,alpha=0.05,NBOOT=50,type="emp")
{
  ###puntual estimation
  #f1x=x;f2x=y;f3x=z
  if(type=="normal"){
    a <- sd(f2x)/sd(f1x)
    b <- (mean(f1x)-mean(f2x))/sd(f1x)
    c <- sd(f2x)/sd(f3x)
    d <- (mean(f3x)-mean(f2x))/sd(f3x)
    puntual<- integrate(f0,a=a,b=b,c=c,d=d,p=0,q=0,lower=-Inf,upper=Inf,subdivisions=10000000)$value
    }
  if(type=="kernel") {
      bwf<-function(x){((4/(3*length(x)))^(1/5))*min(sd(x),IQR(x)/1.39)}
      f1k<-density(f1x,bw=bwf(f1x))
      f1k<-sample(f1k$x,length(f1x),prob=f1k$y,replace=F)
      f2k<-density(f2x,bw=bwf(f2x))
      #f2k<-sample(f2k$x,length(f2x),prob=f2k$y,replace=F)
      f3k<-density(f3x,bw=bwf(f3x))
      f3k<-sample(f3k$x,length(f3x),prob=f3k$y,replace=F)
      f2k<-approxfun(f2k,yleft=0,yright = 0)
      kernel_surface<-function(y){F_emp(f1k,y)*(1-F_emp(f3k,y))*f2k(y)}
      puntual<-integrate(kernel_surface,lower=-Inf,upper=Inf,subdivisions = 1000,rel.tol = 0.0005)$value
  }
  if(type=="emp"){
    lim_int<-function(p1){1-F_emp(f3x,F_emp(f1x,q_emp(f1x,p1)))}
    puntual <- pracma::integral2(roc_emp,0,1,0,lim_int)$Q}
  
  #boot.VUS <- sapply(1:NBOOT,function(jj){
  res0<-numeric(NBOOT)
  x<-matrix(NA,nrow=length(f1x)*2,ncol=NBOOT)
  y<-matrix(NA,nrow=length(f2x)*2,ncol=NBOOT)
  z<-matrix(NA,nrow=length(f3x)*2,ncol=NBOOT)
  for(i in 1:NBOOT){
  x[,i]<-sample(f1x,length(f1x)*2,replace=T);
  y[,i]<-sample(f2x,length(f2x)*2,replace=T);
  z[,i]<-sample(f3x,length(f3x)*2,replace=T);
  if(type=="normal"){
     # lim_int_n<-function(p1){1-pnorm(pnorm(qnorm(p1,mean(x[,i],na.rm=T),
      #                                            sd(x[,i],na.rm=T)),mean(x[,i],na.rm=T),sd(x[,i],na.rm=T)),
          #                           mean(z[,i],na.rm=T),sd(z[,i],na.rm=T))}
      #res0<-pracma::integral2(roc_normal_est,0,1,0,lim_int_n)$Q
      ###params
      a <- sd(y[,i])/sd(x[,i])
      b <- (mean(x[,i])-mean(y[,i]))/sd(x[,i])
      c <- sd(y[,i])/sd(z[,i])
      d <- (mean(z[,i])-mean(y[,i]))/sd(z[,i])
      res0[i] <- integrate(f0,a=a,b=b,c=c,d=d,p=0,q=0,lower=-Inf,upper=Inf,subdivisions=10000000)$value
      
    }
    if(type=="kernel"){
      bwf<-function(x){((4/(3*length(x)))^(1/5))*min(sd(x),IQR(x)/1.39)}
      f1k<-density(x[,i],bw=bwf(x[,i]))
      f1k<-sample(f1k$x,length(x[,i]),prob=f1k$y,replace=F)
      f2k<-density(y[,i],bw=bwf(y[,i]))
      f3k<-density(z[,i],bw=bwf(z[,i]))
      f3k<-sample(f3k$x,length(z[,i]),prob=f3k$y,replace=F)
      f2k<-approxfun(f2k,yleft=0,yright = 0)
      kernel_surface<-function(y){F_emp(f1k,y)*(1-F_emp(f3k,y))*f2k(y)}
      res0[i]<-integrate(kernel_surface,lower=-Inf,upper=Inf,subdivisions = 1000,rel.tol = 0.0005)$value
      
    }
    if(type=="emp"){
    
      lim_int<-function(p1){1-F_emp(z[,i],F_emp(x[,i],q_emp(x[,i],p1)))}
      roc_emp<-function(p3,p1){ F_emp(y[,i],q_emp(z[,i],1-p3))-F_emp(y[,i],q_emp(x[,i],p1))}
      res0[i] <- pracma::integral2(roc_emp,0,1,0,lim_int)$Q
    }
  }
  ###caculate variance
  var0 <- var(res0,na.rm=T)
  prob0 <- alpha/2
  ###CI
  lower <- quantile(boot.VUS,prob=prob0)
  upper <- quantile(boot.VUS,prob=1-prob0)
  CI <- c(lower,upper);
  names(CI) <- c(paste(prob0*100,"%",sep=""),paste(100-prob0*100,"%",sep=""))
  return(list(variance=var0,CI=CI,estimate=puntual))
}

