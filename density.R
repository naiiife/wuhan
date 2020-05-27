
#### Y = X + Z
set.seed(1)
n = 200
#(a) EG=5,VG=8.33,EI=7.5,VI=12.5
alphaG = 3
betaG = 0.6
alphaI = 4.5
betaI = 0.6
h = 1.51
#(b) EG=5,VG=8.33,EI=7.5,VI=6.25
alphaG = 3
betaG = 0.6
alphaI = 9
betaI = 1.2
h = 1.20
#(c) EG=5,VG=4.17,EI=7.5,VI=12.5
alphaG = 6
betaG = 1.2
alphaI = 4.5
betaI = 0.6
h = 1.38
#(d) EG=5,VG=4.17,EI=7.5,VI=6.25
alphaG = 6
betaG = 1.2
alphaI = 9
betaI = 1.2
h = 1.25

X = rgamma(n,alphaG,betaG)
Z1 = rgamma(n,alphaI,betaI)
Z2 = rgamma(n,alphaI,betaI)

Y = X + Z1 - Z2



phi.Z <- function(t,a=alphaI,b=betaI){
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

## hh should be chosen be hand!
## Larger hh, flatter.

sinx <- function(x)  sin(x)/x
sini <- function(c)  integrate(sinx,0,c)
phi.K <- function(x,c){
if (c==0){
sinic = 0
S11 = 0.5
S21 = -16/5/pi
S12 = 16/35/pi
S22 = -0.5
}
else {
sinic = sini(c)$value
S11 = (2*c*(120-2*c^2+c^4)*cos(c)+2*(-120+42*c^2+c^4)*sin(c)+c^6*(pi+2*sinic))/(2*c^6*pi)
S21 = 48*(3*c*cos(c)+(-3+c^2)*sin(c))/(c^5*pi)
S12 = -48*(c*(-15+c^2)*cos(c)+3*(5-2*c^2)*sin(c))/(c^7*pi)
S22 = -(2*c*(840-50*c^2+c^4)*cos(c)+2*(-840+330*c^2+c^4)*sin(c)+c^6*(pi+2*sinic))/(2*c^6*pi)
}
a = solve(matrix(c(S11,S21,S12,S22),2,2))[,1]
return( complex(1,a[1],-a[2]*x)*(1-x^2)^3*I(abs(x)<=1) )
}

M = 1
dt = 0.01
t0 = seq(-M,M,dt)
phiZt = phi.Z(t0)

f.est <- function(x){
termj = rep(0,n)
for (j in 1:n){
yj = Y[j]
phiKtc = rep(0,length(t0))
for (i in 1:length(t0)){
phiKtc[i] = phi.K(t0[i]*hh,yj/hh)
}
termj[j] = dt * sum(Re(exp(complex(1,0,t0*(yj-x)))*phiKtc) / phiZt)
}
return(mean(termj)/pi/2)
}


xrange = rep(NA,70)
dens = rep(NA,70)
for (k in 1:70){
print(k)
dens[k] = f.est(k/3-1/3)
xrange[k] = k/3-1/3
}

dens = dens/sum(dens/3)

loglik <- function(pa){
a = pa[1]
b = pa[2]
RE = - sum( log(dgamma(Y[Y>0],shape=a,rate=b)) )
return(RE)
}
par0 <- optim(par=c(4,2),loglik)$par

hist(Y,prob=T,ylim=c(0,0.22),xlab='Day',main='Histogram of S')
points(xrange,dgamma(xrange,shape=par0[1],rate=par0[2]),type='l',lwd=2,col='darkcyan')
points(xrange,dgamma(xrange,alphaG,betaG),type='l',lwd=2)
points(xrange,dens,type='l',col='red',lwd=2)
#mtext('E(G)=5, V(G)=4.17, E(I)=7.5, V(I)=12.5')
savePlot('simuG4',type='png')
savePlot('simuG4',type='eps')


