library(metapoppkg)
## R0 for mle
po_mle <- li23v3(U=5,version="MLEperiod3")
p_mle <- coef(po_mle)
R0v3(p_mle,be=T)
R0v3(p_mle,be=F)

## R0 for li20period2
po_li <- li23v3(U=5,version="li20period2")
p_li <- coef(po_li)
R0v3(p_li,be=T)
R0v3(p_li,be=F)
