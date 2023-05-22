
library(metapoppkg)
i <- 1 # set to 2 for larger test, 3 for big test

set.seed(99)
p <- li20v3(U=switch(i,5,20,373))
s <- simulate(p)
obs(s)[1:3,1:5]

theta <- coef(p)

theta_trans <- partrans(p,theta,dir="toEst")
theta_trans

b <- bpfilter(p,block_size=1,Np=switch(i,10,100,1000))
logLik(b)

# to investigate NA results, look at
# b@block.cond.loglik