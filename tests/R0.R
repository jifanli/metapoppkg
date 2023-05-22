library(metapoppkg)
## R0 for mle
po_mle <- li20v3(U=5,version="MLEperiod3")
p_mle <- coef(po_mle)
R0(p_mle,be=T)
R0(p_mle,be=F)

## R0 for li20period2
po_li <- li20v3(U=5,version="li20period2")
p_li <- coef(po_li)
R0(p_li,be=T)
R0(p_li,be=F)

if(any(abs(c(R0(p_li,be=T),R0(p_li,be=F)) - c(2.38,1.34))>0.1)) stop("R0 not matching Li et al")





