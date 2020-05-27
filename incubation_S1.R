
result = NULL

aa = 5
bb = 0.8

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
  REh = pgamma(k[posit],shape=a+1,rate=b)+k[posit]*b/a*(1-pgamma(k[posit],shape=a,rate=b))
  if (p>0.999) RE[posit] = REh
  else{
    REf = pgamma(k[posit],shape=a,rate=b)
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
  a.vec.qin = NULL
  b.vec.qin = NULL
  q.vec.qin = NULL
  a.vec.ic = NULL
  b.vec.ic = NULL
  q.vec.ic = NULL

  
  for (run in 1:1000){
    
    #### Genarating data
    set.seed(run)
    x = rep(NA,m)
    c = rep(NA,m)
    for (i in 1:m){
      select = rbinom(1,1,pp)
      if (select==0){
	x[i] = rgamma(1,aa,bb)
	c[i] = round(runif(1,0,30),0)
	}
      else{
        while (TRUE){
          C = runif(1,0,30)		#departure time
          Y = rgamma(1,aa,bb)		#symptom time
          if (Y>C) break
        }
        x[i] = round(Y-C,4)		#observed forward time
	  c[i] = round(runif(1,0,30),0)	#interval length
      }
    }
    #### p: censor probobility
    
    x=round(x,0)
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

    loglikc <- function(pa){
	a = pa[1]
	b = pa[2]
	p = pa[3]
	#P = p*(FI(x+c,a,b)-FI(x,a,b))+(1-p)*f(x,a,b)
	P1 = FI(x+c,a,b)-FI(x,a,b)
	P2 = f(x,a,b)
	P1[P1<0.00001]=0.00001
	P2[P2<0.00001]=0.00001
      RE = - sum( p*log(P1)+(1-p)*log(P2) )
      if (RE>10000) RE=10000
      return(RE)
    }

    
    ## Init
#    a = 4
#    b = 0.5
#    p = 0.9
#    g <- rep(p,m)
#    ## Iter
#    for (k in 1:1000){
#      g = p*(F(x+0.5,a,b,1)-F(x-0.5,a,b,1)) / (F(x+0.5,a,b,p)-F(x-0.5,a,b,p))
#      logl <- function(pa){
#        a = pa[1]
#        b = pa[2]
#        RE = - sum( (1-g)*log(FI(x+0.5,a,b)-FI(x-0.5,a,b)) + g*log(F(x+0.5,a,b,1)-F(x-0.5,a,b,1)) )
#        return(RE)
#      }
#      a0 = a
#      b0 = b
#      par1 <- optim(par=c(a0,b0),logl)$par
#      a = par1[1]
#      b = par1[2]
#      p = mean(g)
#      if (sum((a-a0)^2+(b-b0)^2)<1e-6) break
#    }
#    par0 <- optim(par=c(4,0.5),Loglik0)$par
#    if (Loglik(c(par1,p))>Loglik0(par0)) {
#      a = par0[1]
#      b = par0[2]
#      p = 1
#    }
#    a.vec.em=append(a.vec.em,a)
#    b.vec.em=append(b.vec.em,b)
#    q.vec.em=append(q.vec.em,1-p)
    
    #par0 <- optim(par=c(5,0.5,0.8),Loglik)$par
    par1 <- optim(par=c(4,0.5,0.8),Loglik,method='L-BFGS-B',
                  lower=c(1.1,0.1,0),upper=c(10,2,1))$par
    #par0 <- optim(par=c(5,0.5),Loglik0,method='L-BFGS-B',
    #              lower=c(1.1,0.1),upper=c(10,2))$par
    #loglik(par1)
    #loglik0(par0)
    a = par1[1]
    b = par1[2]
    p = par1[3]
    a.vec=append(a.vec,a)
    b.vec=append(b.vec,b)
    q.vec=append(q.vec,1-p)
    par2 <- optim(par=c(4,0.5,0.8),loglik,method='L-BFGS-B',
                  lower=c(1.1,0.1,0),upper=c(10,2,1))$par
    a = par2[1]
    b = par2[2]
    p = par2[3]
    a.vec.qin=append(a.vec.qin,a)
    b.vec.qin=append(b.vec.qin,b)
    q.vec.qin=append(q.vec.qin,1-p)
    par3 <- optim(par=c(4,0.5,0.8),loglikc,method='L-BFGS-B',
                  lower=c(1.1,0.1,0),upper=c(20,5,1))$par
    a = par3[1]
    b = par3[2]
    p = par3[3]
    a.vec.ic=append(a.vec.ic,a)
    b.vec.ic=append(b.vec.ic,b)
    q.vec.ic=append(q.vec.ic,1-p)

  }
  
  res=c(paste(round(mean(a.vec),2), ' (', round(sd(a.vec),2), ')', sep=''),
        paste(round(mean(b.vec),2), ' (', round(sd(b.vec),2), ')', sep=''),
        paste(round(mean(q.vec),2), ' (', round(sd(q.vec),2), ')', sep=''),
        paste(round(mean(a.vec.qin),2), ' (', round(sd(a.vec.qin),2), ')', sep=''),
        paste(round(mean(b.vec.qin),2), ' (', round(sd(b.vec.qin),2), ')', sep=''),
        paste(round(mean(q.vec.qin),2), ' (', round(sd(q.vec.qin),2), ')', sep=''),
        paste(round(mean(a.vec.ic),2), ' (', round(sd(a.vec.ic),2), ')', sep=''),
        paste(round(mean(b.vec.ic),2), ' (', round(sd(b.vec.ic),2), ')', sep=''),
        paste(round(mean(q.vec.ic),2), ' (', round(sd(q.vec.ic),2), ')', sep=''))

  result = rbind(result, res)
  colnames(result) = c('alpha','beta','pi','alpha_q','beta_q','pi_q',
	'alpha_ic','beta_ic','pi_ic')
  print(res)
}
write.table(result, 'setting_gamma.csv')
