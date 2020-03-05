
result = NULL

aa = 2
bb = 9

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
  if (p>0.999) RE[posit] = REh
  else{
    REf = pweibull(k[posit],shape=a,scale=b)
    RE[posit] = REf*(1-p)+REh*p
  }
  return(RE)
}

for (w in 1:6){
  if (w%%3==1) m=600
  if (w%%3==2) m=1200
  if (w%%3==0) m=1800
  if (w<=3) pp=1
  if (w>=4) pp=0.8
  
  a.vec = NULL
  b.vec = NULL
  q.vec = NULL
  a.vec.em = NULL
  b.vec.em = NULL
  q.vec.em = NULL
  
  for (run in 1:1000){
    
    #### Genarating data
    set.seed(run)
    x = rep(NA,m)
    for (i in 1:m){
      select = rbinom(1,1,pp)
      if (select==0) x[i]=rweibull(1,aa,bb)
      else{
        while (TRUE){
          C = runif(1,0,30)
          Y = rweibull(1,aa,bb)
          if (Y>C) break
        }
        x[i] = round(Y-C,4)
      }
    }
    #### p: censor probobility
    #### g=1: censored  g=0: exact
    
    x = round(x,0)
    x[x>25] = 25
    
    
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
    
    ## Init
    a = 1.5
    b = 10
    p = 0.9
    g <- rep(p,m)
    ## Iter
    for (k in 1:1000){
      g = p*(F(x+0.5,a,b,1)-F(x-0.5,a,b,1)) / (F(x+0.5,a,b,p)-F(x-0.5,a,b,p))
      logl <- function(pa){
        a = pa[1]
        b = pa[2]
        RE = - sum( (1-g)*log(FI(x+0.5,a,b)-FI(x-0.5,a,b)) + g*log(F(x+0.5,a,b,1)-F(x-0.5,a,b,1)) )
        return(RE)
      }
      a0 = a
      b0 = b
      par1 <- optim(par=c(a0,b0),logl)$par
      a = par1[1]
      b = par1[2]
      p = mean(g)
      if (sum((a-a0)^2+(b-b0)^2)<1e-6) break
    }
    par0 <- optim(par=c(2,8),Loglik0)$par
    if (Loglik(c(par1,p))>Loglik0(par0)) {
      a = par0[1]
      b = par0[2]
      p = 1
    }
    a.vec.em=append(a.vec.em,a)
    b.vec.em=append(b.vec.em,b)
    q.vec.em=append(q.vec.em,1-p)
    
    #par0 <- optim(par=c(a,b,p),Loglik)$par
    par1 <- optim(par=c(2,10,0.8),Loglik,method='L-BFGS-B',
                  lower=c(1.1,2,0),upper=c(5,15,1))$par
    #par0 <- optim(par=c(2,10),loglik0,method='L-BFGS-B',
    #               lower=c(1.1,2),upper=c(5,15))$par
    #loglik(par1)
    #loglik0(par0)
    a = par1[1]
    b = par1[2]
    p = par1[3]
    a.vec=append(a.vec,a)
    b.vec=append(b.vec,b)
    q.vec=append(q.vec,1-p)
  }
  
  res=c(paste(round(mean(a.vec),3), ' (', round(sd(a.vec),3), ')', sep=''),
        paste(round(mean(b.vec),3), ' (', round(sd(b.vec),3), ')', sep=''),
        paste(round(mean(q.vec),3), ' (', round(sd(q.vec),3), ')', sep=''),
        paste(round(mean(a.vec.em),3), ' (', round(sd(a.vec.em),3), ')', sep=''),
        paste(round(mean(b.vec.em),3), ' (', round(sd(b.vec.em),3), ')', sep=''),
        paste(round(mean(q.vec.em),3), ' (', round(sd(q.vec.em),3), ')', sep=''))
  result = rbind(result, res)
  colnames(result) = c('mu','sigma','pi','mu','sigma','pi')
  print(res)
}
write.table(result, 'setting_weibull.csv')
