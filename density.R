
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

serial <- read.table('serial.csv',sep=',',head=T)
Y <- serial$serial
n <- length(Y)


phi.Z <- function(t,a=alpha,b=beta){
num = b^2
den = b^2+t^2
return( (num/den)^a )
}

#### Parametric ####

#Y.id <- as.numeric(names(table(Y)))
#Y.num <- as.numeric(table(Y))
#dt = 0.001
#t <- seq(-5,5,dt)
#loglik <- function(pa){
#a = pa[1]
#b = pa[2]
#fS = rep(NA,length(Y.id))
#for (j in 1:length(Y.id)){
#int.cal = phi.Z(t)*complex(1,1,-t/b)^(-a)*exp(complex(1,0,-t*Y.id[j]))
#fS[j] = log(max(sum(Re(int.cal))*dt/2/pi,0.0001))
#}
#return(-sum(fS*Y.num))
#}
#par0 <- optim(par=c(2,1),loglik)$par
#a = par0[1]
#b = par0[2]
#par0


#### Nonparametric ####

hh1 = 1.6
hh4 = 6
hh5 = 3
## hh should be chosen be hand!
## Larger hh, flatter.

phi.K1 <- function(x){	# Fejer / Vallee Poussin
return( max(0, 1-abs(x)) )
}
phi.K4 <- function(x){	# trapezoid
if (abs(x)<1) return(1)
else return( max(0, 2-abs(x)) )
}
phi.K5 <- function(x){ 	# Cesaro
return( max(0, 1-x^4) )
}

M = 2
dt = 0.0001
t0 = seq(0,M,dt)
#phiKt1 = rep(0,length(t0))
phiKt4 = rep(0,length(t0))
phiKt5 = rep(0,length(t0))
for (i in 1:length(t0)){
#phiKt1[i] = phi.K1(t0[i]*hh1)
phiKt4[i] = phi.K4(t0[i]*hh4)
phiKt5[i] = phi.K5(t0[i]*hh5)
}
phiZt = phi.Z(t0)

f.est4 <- function(x){
termj = rep(0,n)
for (j in 1:n){
termj[j] = dt * sum(cos(t0*(Y[j]-x)) * phiKt4 / phiZt)
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

#xrange = rep(NA,60)
#dens1 = rep(NA,60)
#dens5 = rep(NA,60)
#for (k in 1:60){
#dens1[k] = f.est1(k/2-5)
#dens5[k] = f.est5(k/2-5)
#xrange[k] = k/2-5
#}
#loglik <- function(pa){
#a = pa[1]
#b = pa[2]
#RE = - sum( log(dgamma(Y[Y>0],shape=a,rate=b)) )
#return(RE)
#}
#par0 <- optim(par=c(4,2),loglik)$par
#hist(Y,prob=T,ylim=c(0,0.15),xlab='Day',main='Histogram of S')
#points(xrange,dgamma(xrange,shape=par0[1],rate=par0[2]),type='l',col='brown')
#points(xrange,dgamma(xrange,6.4,0.8),type='l',lwd=2)
#points(xrange,dens5,type='l',col='red',lwd=2)
#points(xrange,dens1,type='l',col='blue',lwd=2)
#mtext('E(G)=8, V(G)=10, E(I)=6, V(I)=12')
#savePlot('2b',type='png')
#savePlot('2b',type='eps')



