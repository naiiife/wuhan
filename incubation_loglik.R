#### Real data

data0 <- read.table('forwardtime.csv',sep=',',head=T)
m = nrow(data0)
x = data0$forward

###############gamma incubation
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
    RE[posit] = REf*p+REh*(1-p)
    return(RE)
}
###############weibull incubation
f <- function(y,a,b){
  return(dweibull(y,shape=a,scale=b))
}	##rate=1/b
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
    RE[posit] = REf*p+REh*(1-p)
    return(RE)
  }
}
###############lognormal incubation
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
    RE[posit] = REf*p+REh*(1-p)
    return(RE)
}



####minus log-likelihood
loglik <- function(pa,p){
      a = pa[1]
      b = pa[2]
      if (p<0.001) P=h(x,a,b)
      else P=(1-p)*h(x,a,b)+p*f(x,a,b)
      P[P<0.000001]=0.000001
      RE = - sum( log(P) )
      if (RE>100000) RE=100000
      return(RE)
}
Loglik <- function(pa,p){
      a = pa[1]
      b = pa[2]
      P = F(x+0.5,a,b,p)-F(x-0.5,a,b,p)
      P[P<0.000001]=0.000001
      RE = - sum( log(P) )
      if (RE>100000) RE=100000
      return(RE)
}


##The difference of loglik is qchisq(0.9,1)/2=1.353

#########################gamma
#par1 <- optim(par=c(4,0.5),loglik,method='L-BFGS-B',p=0,
#                  lower=c(1.1,0.1),upper=c(10,3))$par
#-loglik(par1,0)		#-3203.45
#par0 <- optim(par=c(4,0.5),loglik,method='L-BFGS-B',p=0.188,
#                  lower=c(1.1,0.1),upper=c(10,3))$par
#-loglik(par0,0.188)	#-3204.80
#pi=0.188
########################weibull
#par1 <- optim(par=c(2,10),loglik,method='L-BFGS-B',p=0,
#                  lower=c(1.1,2),upper=c(5,15))$par
#-loglik(par1,0)		#-3202.96
#par0 <- optim(par=c(2,10),loglik,method='L-BFGS-B',p=0.174,
#                  lower=c(1.1,2),upper=c(5,15))$par
#-loglik(par0,0.174)	#-3204.3
#pi=0.174
########################lognormal
#par1 <- optim(par=c(2,0.4),loglik,method='L-BFGS-B',p=0,
#                lower=c(1,0.1),upper=c(5,1))$par
#-loglik(par1,0)		#-3205.92
#par0 <- optim(par=c(2,0.4),loglik,method='L-BFGS-B',p=0.104,
#                lower=c(1,0.1),upper=c(5,1))$par
#-loglik(par0,0.104)	#-3207.27
#pi=0.104



L =rep(NA,20)
for (pp in 1:20){
p=(pp-1)/40
par1 <- optim(par=c(2,0.4),Loglik,method='L-BFGS-B',p=p,
                  lower=c(1.1,0.1),upper=c(10,2))$par
L[pp]=-Loglik(par1,p)
}
plot((0:19)/40,L,type='l',ylab='log-likelihood',
    main='log-normal incubation',xlab=expression(pi))
abline(h=max(L)-qchisq(0.90,1)/2,col='blue',lty=2)
savePlot('Lognorm_l',type='png')
savePlot('Lognorm_l',type='eps')
