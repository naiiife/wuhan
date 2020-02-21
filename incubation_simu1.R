
result = NULL

aa = 6
bb = 0.8

##########################################gamma incubation
f <- function(y,a,b){
  return(dgamma(y,shape=a,rate=b))
}
h <- function(y,a,b){
  RE = b/a*(1-pgamma(y,shape=a,rate=b))
  return(RE)
}
##########################################weibull incubation
#f <- function(y,a,b){
#  return(dweibull(y,shape=a,scale=b))
#}
#h <- function(y,a,b){
#  RE = (1-pweibull(y,shape=a,scale=b))/gamma(1+1/a)/b
#  return(RE)
#}

##########################################lognormal incubation
#f <- function(y,u,s){
#  return(dnorm(log(y),u,s)/y)
#}
#h <- function(y,u,s){
#  RE = exp(-u-s^2/2)*(1-pnorm(log(y),u,s))
#  return(RE)
#}


for (w in 1:12){
if (w%%6==1){
m = 400
r = 0
}
if (w%%6==2){
m = 400
r = 20
}
if (w%%6==3){
m = 800
r = 0
}
if (w%%6==4){
m = 800
r = 40
}
if (w%%6==5){
r = 1200
n = 0
}
if (w%%6==0){
m = 1200
r = 60
}
if (w<=6) pp=0.9
if (w>=7) pp=0.7

a.vec = NULL
b.vec = NULL
q.vec = NULL

for (run in 1:1000){

#### Genarating data
set.seed(run)
n = m + r
x = rep(NA,n)
for (i in 1:m){
  select = rbinom(1,1,pp)
  if (select==0) x[i]=round(rgamma(1,shape=aa,rate=bb),4)
  else{
    while (TRUE){
      C = runif(1,0,30)
      Y = rgamma(1,shape=aa,rate=bb)
      if (Y>C) break
    }
    x[i] = round(Y-C,4)
  }
}
if (r!=0){
x[(m+1):n] = rgamma(r,shape=aa,rate=bb)
}
#### p: censor probobility
#### g=1: censored  g=0: exact

x[x<0.01] = 0.01
x[x>20] = 20


## Init
a = aa
b = bb
p = 0.8
g <- rep(0,n)
g[1:m] = p
## Iter
for (k in 1:1000){
g[1:m] = p*h(x[1:m],a,b)/ (p*h(x[1:m],a,b)+(1-p)*f(x[1:m],a,b))
loglik <- function(pa){
a = pa[1]
b = pa[2]
RE = - sum( (1-g)*log(f(x,a,b)) + g*log(h(x,a,b)) )
return(RE)
}
a0 = a
b0 = b
par0 <- optim(par=c(a0,b0),loglik)$par
#par0 <- optim(par=c(a0,b0),loglik,method='L-BFGS-B',lower=c(0.001,0.001))$par
a = par0[1]
b = par0[2]
p = mean(g[1:m])
if (sum((a-a0)^2+(b-b0)^2)<1e-6) break
}

a.vec=append(a.vec,a)
b.vec=append(b.vec,b)
q.vec=append(q.vec,1-p)
}

res=c(paste(round(mean(a.vec),3), ' (', round(sd(a.vec),3), ')', sep=''),
paste(round(mean(b.vec),3), ' (', round(sd(b.vec),3), ')', sep=''),
paste(round(mean(q.vec),3), ' (', round(sd(q.vec),3), ')', sep=''))
result = rbind(result, res)
colnames(result) = c('alpha','beta','pi')
print(result)
}

write.table(result, 'setting_gamma.csv')
