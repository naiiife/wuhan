
#### Real data
data0 <- read.table('forwardtime.csv',sep=',',head=T)
m = nrow(data0)
x0 = data0$forward

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



RUN = 1000
a.vec = rep(NA,RUN)
b.vec = rep(NA,RUN)
q.vec = rep(NA,RUN)
ave.vec = rep(NA,RUN)
Q1.vec = rep(NA,RUN)
Q2.vec = rep(NA,RUN)
Q3.vec = rep(NA,RUN)
Q4.vec = rep(NA,RUN)
Q5.vec = rep(NA,RUN)


for (simu in 1:RUN){
set.seed(simu)
samp = sample(length(x0),replace=T)
x = x0[samp]

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

#par1 <- optim(par=c(4,0.5,0.8),Loglik,method='L-BFGS-B',
#              lower=c(1.1,0.1,0),upper=c(10,2,1))$par
#par1 <- optim(par=c(2,10,0.8),Loglik,method='L-BFGS-B',
#              lower=c(1.1,2,0),upper=c(5,15,1))$par
par1 <- optim(par=c(2,0.4,0.8),Loglik,method='L-BFGS-B',
              lower=c(1,0.1,0),upper=c(5,1,1))$par
a = par1[1]
b = par1[2]
p = par1[3]

a.vec[simu] = a
b.vec[simu] = b
q.vec[simu] = 1-p
#ave.vec[simu] = b*gamma(1+1/a)
ave.vec[simu] = exp(a+b^2/2)
#Q1.vec[simu] = qweibull(0.25,a,b)
#Q2.vec[simu] = qweibull(0.5,a,b)
#Q3.vec[simu] = qweibull(0.75,a,b)
#Q4.vec[simu] = qweibull(0.9,a,b)
#Q5.vec[simu] = qweibull(0.99,a,b)
Q1.vec[simu] = exp(qnorm(0.25,a,b))
Q2.vec[simu] = exp(qnorm(0.5,a,b))
Q3.vec[simu] = exp(qnorm(0.75,a,b))
Q4.vec[simu] = exp(qnorm(0.9,a,b))
Q5.vec[simu] = exp(qnorm(0.99,a,b))
}
#########################gamma
#round(c(a,b,1-p),2)
#round(a/b,2)
#round(qgamma(c(0.25,0.5,0.75,0.9),a,b),2)
#-round(loglik(c(a,b)),2)
#########################weibull
#round(c(a,b,1-p),2)
#round(b*gamma(1+1/a),2)
#round(qweibull(c(0.25,0.5,0.75,0.9),shape=a,scale=b),2)
#-round(loglik(c(a,b)),2)
########################lognormal
#round(c(a,b,1-p),2)
#round(exp(a+b^2/2),2)
#round(exp(qnorm(c(0.25,0.5,0.75,0.9),a,b)),2)
#-round(loglik(c(a,b)),2)

a.vec = sort(a.vec)
b.vec = sort(b.vec)
q.vec = sort(q.vec)
ave.vec = sort(ave.vec)
Q1.vec = sort(Q1.vec)
Q2.vec = sort(Q2.vec)
Q3.vec = sort(Q3.vec)
Q4.vec = sort(Q4.vec)
Q5.vec = sort(Q5.vec)
round(c(a.vec[26],a.vec[975]), 2)
round(c(b.vec[26],b.vec[975]), 2)
round(c(q.vec[1],q.vec[950]), 2)
round(c(ave.vec[26],ave.vec[975]), 2)
round(c(Q1.vec[26],Q1.vec[975]), 2)
round(c(Q2.vec[26],Q2.vec[975]), 2)
round(c(Q3.vec[26],Q3.vec[975]), 2)
round(c(Q4.vec[26],Q4.vec[975]), 2)
round(c(Q5.vec[26],Q5.vec[975]), 2)

