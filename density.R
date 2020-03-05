
#### Y = X + Z
set.seed(1)
#n = 100
alpha = 4.97
beta = 0.55
#alpha = 3
#beta = 0.5
#X = rgamma(n,6.4,0.8)
#Z1 = rgamma(n,alpha,beta)
#Z2 = rgamma(n,alpha,beta)

#Y = X + Z1 - Z2

interval <- c(4,  6,  2,  5, -3,  4,  2,  8,  7,  8,  4,  1, 13,  2,  1,  
              5,  6,  4,  9,  4,  3,  7,  1,  5,  3,  6, -5, 11,
              4,  6, 11,  8,  4,  6,  5,  0,  6,  7,  9, 10, -2,  5,  2,  
              7,  2,  4,  4,  1,  0,  7,  4,  2,  0,  6,  2,  2,
              6, 13,  1,  7,  6,  2,  3,  0,  0,  2,  2,  6,  3,  4,  3)
Y <- interval
n <- length(Y)


hh1 = 1.6
hh5 = 1.9
## Sample sample size: hh should be chosen be hand!
## Larger hh, flatter.

phi.Z <- function(t,a=alpha,b=beta){
num = b^2
den = b^2+t^2
return( (num/den)^a )
}

phi.K1 <- function(x){	# Vallee Poussin
return( max(0, 1-abs(x)) )
}

phi.K5 <- function(x){ 	# 
return( max(0, 1-x^4) )
}

phi.K2 <- function(x){ 	# 
return( phi.Z(x)^2 )
}


M = 2
dt = 0.0001
t0 = seq(0,M,dt)
phiKt1 = rep(0,length(t0))
phiKt5 = rep(0,length(t0))
for (i in 1:length(t0)){
phiKt1[i] = phi.K1(t0[i]*hh1)
phiKt5[i] = phi.K5(t0[i]*hh5)
}
phiZt = phi.Z(t0)

f.est1 <- function(x){
termj = rep(0,n)
for (j in 1:n){
termj[j] = dt * sum(cos(t0*(Y[j]-x)) * phiKt1 / phiZt)
}
return(mean(termj)/pi)
}
f.est5 <- function(x){
termj = rep(0,n)
for (j in 1:n){
termj[j] = dt * sum(cos(t0*(Y[j]-x)) * phiKt5 / phiZt)
}
return(mean(termj)/pi)
}

xrange = rep(NA,60)
dens1 = rep(NA,60)
dens5 = rep(NA,60)
for (k in 1:60){
dens1[k] = f.est1(k/2-5)
dens5[k] = f.est5(k/2-5)
xrange[k] = k/2-5
}
loglik <- function(pa){
a = pa[1]
b = pa[2]
RE = - sum( log(dgamma(Y[Y>0],shape=a,rate=b)) )
return(RE)
}
par0 <- optim(par=c(4,2),loglik)$par
hist(Y,prob=T,ylim=c(0,0.15),xlab='Day',main='Histogram of S')
#points(xrange,dgamma(xrange,shape=par0[1],rate=par0[2]),type='l',col='brown')
#points(xrange,dgamma(xrange,6.4,0.8),type='l',lwd=2)
#points(xrange,dens5,type='l',col='red',lwd=2)
#points(xrange,dens1,type='l',col='blue',lwd=2)
#mtext('E(G)=8, V(G)=10, E(I)=6, V(I)=12')
#savePlot('2b',type='png')
#savePlot('2b',type='eps')

xrange = rep(NA,55)
dens5 = rep(NA,55)
for (k in 1:55){
dens5[k] = f.est5(k/5-1.55)
xrange[k] = k/5-1.55
}
dens5 = dens5/sum(dens5*0.2)
hist(Y,prob=T,ylim=c(0,0.21),main='Histogram of serial intervals',xlab='Day')
points(density(Y),type='l',lwd=2)
points(xrange,dens5,type='l',col='red',lwd=2)
mtext('Red line: Estimated density of generation time')
#savePlot('GenerationTime',type='png')
#savePlot('GenerationTime',type='eps')

# From Jan 20 to Jan 27, Data source: China CDC
Hubei <- c(270, 375, 444, 549, 729, 1052, 1423, 2714)
China <- c(291, 440, 571, 830, 1287, 1975, 2744, 4535)
cases <- China-Hubei
newcase <- cases[2:8]-cases[1:7]
tvec = 1:7
m0 = lm(log(newcase)~1+tvec)
r = m0$coef[2]
se = summary(m0)$coef[2,2]
dof = 5
R0 = 1 / sum(exp(-r*xrange)*dens5*0.2)
R0l = 1 / sum(exp(-(r-qt(0.975,dof)*se)*xrange)*dens5*0.2)
R0u = 1 / sum(exp(-(r+qt(0.975,dof)*se)*xrange)*dens5*0.2)
R0b = 1 / sum(exp(-r*Y)/length(Y))

