
i <- 1
U <- 5
N <- 5
Np <- switch(i,20,500)

library(metapoppkg)
m1 <- li23(U=U, version="li20period3", for_ibpf = F, days=N)

set.seed(526)
e1 <- enkf(m1,Np=Np)

round(logLik(e1),4)

set.seed(527)
s1 <- simulate(m1)
v1 <- vunit_measure(s1,states(s1),1:U,1:N,coef(s1))
v1[,,1]

set.seed(528)
i1 <- ienkf(s1,Np=Np,cooling.fraction.50=0.5,Nenkf=2,
  rw.sd = rw_sd(alpha_be1=0.02,Beta_af1=0.02))

round(logLik(i1),4)

