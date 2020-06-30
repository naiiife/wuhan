

#### Real data

data0 <- read.table('forwardtime.csv',sep=',',head=T)
m = nrow(data0)
x = data0$forward

##########################################gamma incubation
f <- function(y,a,b){
 return(dgamma(y,shape=a,rate=b))
}
h <- function(y,a,b){
 RE = b/a*(1-pgamma(y,shape=a,rate=b))
 return(RE)
}
FI <- function(y,a,b) return(pgamma(y,shape=a,rate=b))
F <- function(k,a,b,p){
  RE = rep(0,length(k))
  posit = which(k>0)
    REf = pgamma(k[posit],shape=a,rate=b)
    REh = pgamma(k[posit],shape=a+1,rate=b)+k[posit]*b/a*(1-pgamma(k[posit],shape=a,rate=b))
    RE[posit] = REf*(1-p)+REh*p
    return(RE)
}
##########################################weibull incubation
f <- function(y,a,b){
  return(dweibull(y,shape=a,scale=b))
}
h <- function(y,a,b){
  RE = (1-pweibull(y,shape=a,scale=b))/gamma(1+1/a)/b
  return(RE)
}
FI <- function(y,a,b) return(pweibull(y,shape=a,scale=b))
F <- function(k,a,b,p){
  RE = rep(0,length(k))
  posit = which(k>0)
  REh = pgamma((k[posit]/b)^a,shape=1/a,rate=1)
  if (p==1) {
    RE[posit] = REh
    return(RE)
  }
  else{
    REf = pweibull(k[posit],shape=a,scale=b)
    RE[posit] = REf*(1-p)+REh*p
    return(RE)
  }
}
##########################################lognormal incubation
f <- function(y,u,s){
  RE = rep(0,length(y))
  posit = which(y>0)
  RE[posit]=dnorm(log(y[posit]),u,s)/y[posit]
  return(RE)
}
h <- function(y,u,s){
  RE = exp(-u-s^2/2)*(1-pnorm(log(y),u,s))
  return(RE)
}
FI <- function(y,u,s){
  RE = rep(0,length(y))
  posit = which(y>0)
  RE[posit]=pnorm(log(y[posit]),u,s)
  return(RE)
}
F <- function(k,u,s,p){
  RE = rep(0,length(k))
  posit = which(k>0)
    REf = pnorm(log(k[posit]),u,s)
    REh = pnorm(log(k[posit]),u+s^2,s)+exp(-u-s^2/2)*k[posit]*(1-pnorm(log(k[posit]),u,s))
    RE[posit] = REf*(1-p)+REh*p
    return(RE)
}

    loglik <- function(pa){
      a = pa[1]
      b = pa[2]
      p = pa[3]
      if (p>0.999) P=h(x,a,b)
      else P = p*h(x,a,b)+(1-p)*f(x,a,b)
      P[P<0.00001]=0.00001
      RE = - sum( log(P) )
      if (RE>10000) RE=10000
      return(RE)
    }
    loglik0 <- function(pa){
      a = pa[1]
      b = pa[2]
      P = h(x,a,b)
      P[P<0.00001]=0.00001
      RE = - sum( log(P) )
      if (RE>10000) RE=10000
      return(RE)
    }
    
    Loglik <- function(pa){
      a = pa[1]
      b = pa[2]
      p = pa[3]
      P = F(x+0.5,a,b,p)-F(x-0.5,a,b,p)
      P[P<0.00001]=0.00001
      RE = - sum( log(P) )
      if (RE>10000) RE=10000
      return(RE)
    }
    Loglik0 <- function(pa){
      a = pa[1]
      b = pa[2]
      P = F(x+0.5,a,b,1)-F(x-0.5,a,b,1)
      P[P<0.00001]=0.00001
      RE = - sum( log(P) )
      if (RE>10000) RE=10000
      return(RE)
    }


#########################gamma
par1 <- optim(par=c(4,0.5,0.8),Loglik,method='L-BFGS-B',
                  lower=c(1.1,0.1,0),upper=c(10,2,1))$par
-round(Loglik(par1),2)
a = par1[1]
b = par1[2]
p = par1[3]
round(a/b,2)
round(qgamma(c(0.25,0.5,0.75,0.9,0.99),a,b),2)
########################weibull
par1 <- optim(par=c(2,10,0.8),Loglik,method='L-BFGS-B',
                  lower=c(1.1,2,0),upper=c(5,15,1))$par
-round(Loglik(par1),2)
a = par1[1]
b = par1[2]
p = par1[3]
round(b*gamma(1+1/a),2)
round(qweibull(c(0.25,0.5,0.75,0.9,0.99),shape=a,scale=b),2)
########################lognormal
par1 <- optim(par=c(2,0.4,0.8),Loglik,method='L-BFGS-B',
                lower=c(1,0.1,0),upper=c(5,1,1))$par
-round(Loglik(par1),2)
a = par1[1]
b = par1[2]
p = par1[3]
round(exp(a+b^2/2),2)
round(exp(qnorm(c(0.25,0.5,0.75,0.9,0.99),a,b)),2)

###goodness-of-fit
Ei = rep(NA,17)
Oi = rep(NA,17)
for (k in 1:16){
 Ei[k] = (F(k-0.5,a,b,p)-F(k-1.5,a,b,p))*m
 Oi[k] = (sum(x<k-0.5)-sum(x<k-1.5))
}
Ei[17] = (1-F(16-0.5,a,b,p))*m
Oi[17] = sum(x>=16-0.5)
sum((Oi-Ei)^2/Ei)
1-pchisq(sum((Oi-Ei)^2/Ei),13)
#qchisq(0.95,17-3-1)

seq0=seq(0,22,0.1)
hist(x+0.1,breaks=20,prob=T,ylim=c(0,0.16),
     xlab = "Time from departure to symptoms onset",
     main='log-normal incubation')
#lines(seq0,F(seq0+0.5,a,b,p)-F(seq0-0.5,a,b,p),lwd=1.5)
lines(seq0,p*h(seq0,a,b)+(1-p)*f(seq0,a,b),lwd=1.5)
lines(seq0,f(seq0,a,b),col='blue',lwd=2)
lines(seq0,h(seq0,a,b),col='red',lwd=1.5)
savePlot('Lognorm',type='png')
savePlot('Lognorm',type='eps')


