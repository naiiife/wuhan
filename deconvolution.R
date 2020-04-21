
#### Y = X + Z = X + Z1 - Z2
alpha = 4.97
beta = 0.55
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

hh = 2.25
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
dt = 0.001
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


xrange = rep(NA,75)
dens = rep(NA,75)
for (k in 1:75){
print(k)
dens[k] = f.est(k/3-5)
xrange[k] = k/3-5
}

pdfG = dens[15:57]/sum(dens[15:57]/3)
Grange = xrange[15:57]
hist(Y,prob=T,ylim=c(0,0.13),breaks=15,main='Histogram of serial intervals',xlab='Day')
abline(v=0,lty=4,col='yellow',lwd=1.5)
points(density(Y),type='l',lwd=2)
points(Grange,pdfG,type='l',col='red',lwd=2)
mtext('Red line: Estimated density of generation time')
#savePlot('GenerationTime',type='png')
#savePlot('GenerationTime',type='eps')


# From Jan 21 to Jan 30, Data source: China CDC
cases <- c(65,127,281,558,923,1321,1801,2420,3125,3886)
newcase <- cases[-1]-cases[-length(cases)]
tvec = 1:length(newcase)
m0 = lm(log(newcase)~1+tvec)
r = m0$coef[2]
se = summary(m0)$coef[2,2]
dof = length(newcase)-2
R0 = 1 / sum(exp(-r*Grange)*pdfG/3)
R0l = 1 / sum(exp(-(r-qt(0.975,dof)*se)*Grange)*pdfG/3)
R0u = 1 / sum(exp(-(r+qt(0.975,dof)*se)*Grange)*pdfG/3)
R0b = 1 / sum(exp(-r*Y[Y>=0])/length(Y[Y>=0]))


##hh=2  2.9626 (2.15, 3.86)  R0b=2.1804
