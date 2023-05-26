
library(metapoppkg)
m1 <- li23(U=5, for_ibpf = F, days=5, version="li20period3")

# test bootstrap girf
set.seed(530)
bg1 <- girf(m1,kind="bootstrap",Np=20,Ninter=2,Nguide=10,lookahead=1,tol=1e-5)
round(logLik(bg1),4)

if(0){ # NOT YET DEBUGGED
library(metapoppkg)
m1 <- li23(U=5, for_ibpf = F, days=5, version="li20period3")

# test moment girf
mg1 <- girf(m1,kind="moment",Np=20,Ninter=2,Nguide=10,lookahead=1,tol=1e-5,verbose=TRUE)
round(logLik(mg1),4)
}



