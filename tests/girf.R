
library(metapoppkg)
m1 <- li23(U=5, for_ibpf = F, days=5)
bg1 <- girf(m1,kind="bootstrap",Np=20,Ninter=2,Nguide=10,lookahead=1,tol=1e-5,verbose=TRUE)
round(logLik(bg1),4)




