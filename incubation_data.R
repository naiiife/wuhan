

#### Real data
data0 <- read.table('forward.csv',sep=',',head=T)
attach(data0)
n = nrow(data0)
m = n-sum(label)
x = forward[label==0]
x[(m+1):n] = forward[label==1]
x = x+0.1


##########################################gamma incubation
f <- function(y,a,b){
 return(dgamma(y,shape=a,rate=b))
}
h <- function(y,a,b){
 RE = b/a*(1-pgamma(y,shape=a,rate=b))
 return(RE)
}
##########################################weibull incubatiom
#f <- function(y,a,b){
#  return(dweibull(y,shape=a,scale=b))
#}
#h <- function(y,a,b){
#  RE = (1-pweibull(y,shape=a,scale=b))/gamma(1+1/a)/b
#  return(RE)
#}
##########################################lognormal incubation
f <- function(y,u,s){
  return(dnorm(log(y),u,s)/y)
}
h <- function(y,u,s){
  RE = exp(-u-s^2/2)*(1-pnorm(log(y),u,s))
  RE =return(RE)
}

## Inital for gamma
#a = 6
#b = 0.8
#p = 0.8
## Init for weibull
#a = 2
#b = 7
#p = 0.8
## Inital for lognormal
a = 1.8
b = 0.4
p = 0.8


## Iter
g <- rep(0,n)
g[1:m] = p
for (k in 1:1000){
 #g[1:m] = 1  ##use this line if no mixture  
 g[1:m] = p*h(x[1:m],a,b)/ (p*h(x[1:m],a,b)+(1-p)*f(x[1:m],a,b))
 loglik <- function(pa){
   a = pa[1]
   b = pa[2]
   RE = - sum( (1-g[1:m])*log(f(x[1:m],a,b)) + g[1:m]*log(h(x[1:m],a,b)) )
   #RE = - sum( (1-g)*log(f(x,a,b)) + g*log(h(x,a,b)) )
   return(RE)
 }
 a0 = a
 b0 = b
 par0 <- optim(par=c(a0,b0),loglik)$par
 a = par0[1]
 b = par0[2]
 p = mean(g[1:m])
 if (sum((a-a0)^2+(b-b0)^2)<1e-6) break
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


#seq0=seq(0.1,22,0.1)
#hist(x[1:m],breaks=20,prob=T,ylim=c(0,0.2),main='log-normal incubation',
#     xlab = "Time from departure to symptom onset")
#lines(seq0,p*h(seq0,a,b)+(1-p)*f(seq0,a,b),lwd=1.5)
#lines(seq0,f(seq0,a,b),col='blue',lwd=2)
#lines(seq0,h(seq0,a,b),col='red',lwd=1.5)


