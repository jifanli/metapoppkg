
library(metapoppkg)
i <- 1 # set to 2 for larger test, 3 for big test

DEBUG=FALSE
CHECK_LIK=TRUE

U=switch(i,5,20,373)

set.seed(99)
p <- li23v4(U=U)
s <- simulate(p)
obs(s)[1:3,1:5]

theta <- coef(p)

theta_trans <- partrans(p,theta,dir="toEst")
theta_trans

if(CHECK_LIK){
  b <- bpfilter(p,block_size=1,Np=switch(i,10,100,1000))
  cat("bpfilter loglik:",round(logLik(b),5),"\n")
}

f <- flow(p,rinit(p))


round(f[c("C1","C2"),1,1:5],4)


if(DEBUG){
  # plot to check that skeleton is close to a simulation
  par(mfcol=c(2,1))
  matplot(y=t(s@states[paste0("C",1:U),])+1,x=1:30,log="y",ty="l")
  matplot(y=t(f[paste0("C",1:U),1,])+1,x=1:30,log="y",ty="l")
}


  

