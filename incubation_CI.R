
#### Real data
data0 <- read.table('forward.csv',sep=',',head=T)
attach(data0)
n = nrow(data0)
m = n-sum(label)
x = forward[label==0]
x[(m+1):n] = forward[label==1]
#n = m
x = x+0.1

##########################################gamma incubation

#f <- function(y,a,b){
# return(dgamma(y,shape=a,rate=b))
#}
#h <- function(y,a,b){
# RE = b/a*(1-pgamma(y,shape=a,rate=b))
# return(RE)
#}
##########################################weibull incubation
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
#a = 1.5
#b = 0.5
#p = 0.8

x0 = x
label0 = c(rep(0,m),rep(1,n-m))
RUN = 1000
a.vec = rep(NA,RUN)
b.vec = rep(NA,RUN)
q.vec = rep(NA,RUN)
ave.vec =  rep(NA,RUN)
Q1.vec = rep(NA,RUN)
Q2.vec = rep(NA,RUN)
Q3.vec = rep(NA,RUN)
Q4.vec = rep(NA,RUN)
logL = rep(NA,RUN)


for (simu in 1:RUN){
samp = sample(length(x0),replace=T)
x = x0[samp]
lab = label0[samp]
## Iter
a = 1.5
b = 0.5
p = 0.8
g <- rep(0,n)
g[lab==0] = p
for (k in 1:1000){  
 g[lab==0] = p*h(x[lab==0],a,b)/ (p*h(x[lab==0],a,b)+(1-p)*f(x[lab==0],a,b))
 loglik <- function(pa){
   a = pa[1]
   b = pa[2]
   #RE = - sum( (1-g[1:m])*log(f(x[1:m],a,b)) + g[1:m]*log(h(x[1:m],a,b)) )
   RE = - sum( (1-g)*log(f(x,a,b)) + g*log(h(x,a,b)) )
   return(RE)
 }
 a0 = a
 b0 = b
 par0 <- optim(par=c(a0,b0),loglik)$par
 a = par0[1]
 b = par0[2]
 p = mean(g[lab==0])
 if (sum((a-a0)^2+(b-b0)^2)<1e-6) break
}
a.vec[simu] = a
b.vec[simu] = b
q.vec[simu] = 1-p
ave.vec[simu] = exp(a+b^2/2)
Q1.vec[simu] = exp(qnorm(0.25,a,b))
Q2.vec[simu] = exp(qnorm(0.5,a,b))
Q3.vec[simu] = exp(qnorm(0.75,a,b))
Q4.vec[simu] = exp(qnorm(0.9,a,b))
logL[simu] = -loglik(c(a,b))
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
logL = sort(logL)
round(c(a.vec[51],a.vec[950]), 2)
round(c(b.vec[51],b.vec[950]), 2)
round(c(q.vec[51],q.vec[950]), 2)
round(c(ave.vec[51],ave.vec[950]), 2)
round(c(Q1.vec[51],Q1.vec[950]), 2)
round(c(Q2.vec[51],Q2.vec[950]), 2)
round(c(Q3.vec[51],Q3.vec[950]), 2)
round(c(Q4.vec[51],Q4.vec[950]), 2)
round(c(logL[51],logL[950]), 2)


