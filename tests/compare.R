
# plotting the metapoppkg version of the original li20 model and some variations
# to assess the level of agreement.

library(metapoppkg)

i <- 1

set.seed(99)
p1 <- li23(U=switch(i,5,100))
s1 <- simulate(p1)
#pdf(file="s1%03d.pdf",onefile=F)
pdf(file="compare-s1.pdf",onefile=T)
plot(pomp(s1))
dev.off()

set.seed(99)
p2 <- li20(U=switch(i,5,100))
s2 <- simulate(p2)
#pdf(file="s2%03d.pdf",onefile=F)
pdf(file="compare-s2.pdf",onefile=T)
plot(pomp(s2))
dev.off()

set.seed(99)
p3 <- li23(U=5,version="li20period2")
s3 <- simulate(p3)
#pdf(file="s3%03d.pdf",onefile=F)
pdf(file="compare-s3.pdf",onefile=T)
plot(pomp(s3))
dev.off()

set.seed(99)
p4 <- li23(U=5,version="li20period1")
s4 <- simulate(p4)
#pdf(file="s4%03d.pdf",onefile=F)
pdf(file="compare-s4.pdf",onefile=T)
plot(pomp(s4))
dev.off()

