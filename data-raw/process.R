library(spatPomp)
coordinates=read.csv("data-raw/city_coordinates.csv")
wujiaqushi_index=which(coordinates$Name=="Wujiaqushi (Xinjiang)")
hainan_index=which(coordinates$Name=="Hainan")
remove_index=c(wujiaqushi_index,hainan_index)

mobility=read.csv("data-raw/M.csv",header = F)[-remove_index,-c(wujiaqushi_index+(0:13)*375,hainan_index+(0:13)*375)]

population0=read.csv("data-raw/pop.csv")[-remove_index,]
population0[population0$City=="Ezhou",2] <- 1079353
population0[population0$City=="Panjin",2] <- 1389691
population0[population0$City=="Nanchang",2] <- 6255007
population0[population0$City=="Jindezhen",2] <- 1669057
population0[population0$City=="Pingxiang",2] <- 1804805
population0[population0$City=="Jiujiang",2] <- 4600276
population0[population0$City=="Xinyu",2] <- 1138873
population0[population0$City=="Yingtan",2] <- 1124906
population0[population0$City=="Ganzhou",2] <- 8970014
population0[population0$City=="Ji.an",2] <- 4469176
population0[population0$City=="Yichun.Jiangxi",2] <- 5419575
population0[population0$City=="Fuzhou.Jiangxi",2] <- 3912312
population0[population0$City=="Shangrao",2] <- 6491088
population0[population0$City=="Yinchuan",2] <- 2859074
population0[population0$City=="Shizuishan",2] <- 730400
population0[population0$City=="Urumuqi",2] <- 4054000
population0[population0$City=="Kelamayi",2] <- 462347

# Maybe later I can write a function to make it shorter
population0$City[population0$City=="Xintai"] <- "Xingtai"
population0$City[population0$City=="Yinkou"] <- "Yingkou"
population0$City[population0$City=="Zhaoyang"] <- "Chaoyang"
population0$City[population0$City=="Liushui"] <- "Lishui"
population0$City[population0$City=="Wuhui"] <- "Wuhu"
population0$City[population0$City=="Saming"] <- "Sanming"
population0$City[population0$City=="Nannig"] <- "Nanning"
population0$City[population0$City=="Zhanzhou"] <- "Danzhou"
population0$City[population0$City=="Dengmai"] <- "Chengmai"
population0$City[population0$City=="Yuedong"] <- "Ledong"
population0$City[population0$City=="Changdu"] <- "Chamdo"
population0$City[population0$City=="Hanzhou"] <- "Hanzhong"
population0$City[population0$City=="Dongxi"] <- "Dingxi"
population0$City[population0$City=="Wuzhou.Ningxia"] <- "Wuzhong.Ningxia"

incidence=read.csv("data-raw/Incidence.csv")[,-(remove_index+1)]
colnames(incidence)[colnames(incidence)=="Xintai"] <- "Xingtai"
colnames(incidence)[colnames(incidence)=="Yinkou"] <- "Yingkou"
colnames(incidence)[colnames(incidence)=="Zhaoyang"] <- "Chaoyang"
colnames(incidence)[colnames(incidence)=="Liushui"] <- "Lishui"
colnames(incidence)[colnames(incidence)=="Wuhui"] <- "Wuhu"
colnames(incidence)[colnames(incidence)=="Saming"] <- "Sanming"
colnames(incidence)[colnames(incidence)=="Nannig"] <- "Nanning"
colnames(incidence)[colnames(incidence)=="Zhanzhou"] <- "Danzhou"
colnames(incidence)[colnames(incidence)=="Dengmai"] <- "Chengmai"
colnames(incidence)[colnames(incidence)=="Yuedong"] <- "Ledong"
colnames(incidence)[colnames(incidence)=="Changdu"] <- "Chamdo"
colnames(incidence)[colnames(incidence)=="Hanzhou"] <- "Hanzhong"
colnames(incidence)[colnames(incidence)=="Dongxi"] <- "Dingxi"
colnames(incidence)[colnames(incidence)=="Wuzhou.Ningxia"] <- "Wuzhong.Ningxia"

coordinates=coordinates[-remove_index,]
coordinates$Name[coordinates$Name=="Xintai"] <- "Xingtai"
coordinates$Name[coordinates$Name=="Yinkou"] <- "Yingkou"
coordinates$Name[coordinates$Name=="Zhaoyang"] <- "Chaoyang"
coordinates$Name[coordinates$Name=="Liushui"] <- "Lishui"
coordinates$Name[coordinates$Name=="Wuhui"] <- "Wuhu"
coordinates$Name[coordinates$Name=="Saming"] <- "Sanming"
coordinates$Name[coordinates$Name=="Nannig"] <- "Nanning"
coordinates$Name[coordinates$Name=="Zhanzhou"] <- "Danzhou"
coordinates$Name[coordinates$Name=="Dengmai"] <- "Chengmai"
coordinates$Name[coordinates$Name=="Yuedong"] <- "Ledong"
coordinates$Name[coordinates$Name=="Changdu"] <- "Chamdo"
coordinates$Name[coordinates$Name=="Hanzhou"] <- "Hanzhong"
coordinates$Name[coordinates$Name=="Dongxi"] <- "Dingxi"
coordinates$Name[coordinates$Name=="Wuzhou.Ningxia"] <- "Wuzhong.Ningxia"

coordinates[coordinates$Name=="Wujiaqu (Xinjiang)",2:3]=c(87.543, 44.167)
coordinates[coordinates$Name=="Beitun",2:3]=c(87.818056, 47.358333)

U=nrow(coordinates)

mobility=t(mobility)
mobility=as.data.frame(mobility)
mobility=melt(mobility)$value
dim(mobility)=c(U,U*14)
mobility=t(mobility)

colnames(mobility)=population0$City
rownames(mobility)=paste0(rep(population0$City,each=14),1:14)

date <- seq(from=1,to=30,by=1) # populations are changing due to human movement
incidence$Date = date
daterep = rep(date,U)
mobin = rowSums(mobility)
dim(mobin) = c(14,U)
mobin = t(mobin)
mobout=read.csv("data-raw/M.csv",header = F)[-remove_index,-c(wujiaqushi_index*14-14+1:14,hainan_index*14-14+1:14)]
mobout = colSums(mobout)
dim(mobout)=c(U,14)

population_value=population0$Population
population=rep(0,15*U)
for (l in 1:U) {
  population[l*15-14]=population_value[l]
}
for (q in 2:15) {
  for (l in 1:U) {
    # 1.36 is "theta", the multiplicative factor to adjust mobility data, and we use Li et al.'s result 1.36 here
    population[l*15-15+q]=population[l*15-16+q]+1.36*mobin[l,q-1]-1.36*mobout[l,q-1]
    if(population[l*15-15+q]<0.6*population_value[l]){ # same as Li et al. did
      population[l*15-15+q]=0.6*population_value[l]
    }
  }
}

dim(population)=c(15,U)
popuaf = rep(population[15,],each=15) # we assume after the lockdown, there is no human movement
dim(popuaf)=c(15,U)
population=rbind(population,popuaf)
dim(population)=c(30*U,1)
population=cbind(population0[rep(1:nrow(population0),each=30),1],population)

coordinates=coordinates[,-1]

# Distance between two points on a sphere radius R
# Adapted from geosphere package, which has been cited in the package
distHaversine <- function (p1, p2, r = 6378137)
{
  toRad <- pi/180
  p1 <- p1 * toRad
  p2 <- p2 * toRad
  p <- cbind(p1[, 1], p1[, 2], p2[, 1], p2[, 2], as.vector(r))
  dLat <- p[, 4] - p[, 2]
  dLon <- p[, 3] - p[, 1]
  a <- sin(dLat/2) * sin(dLat/2) + cos(p[, 2]) * cos(p[, 4]) *
    sin(dLon/2) * sin(dLon/2)
  a <- pmin(a, 1)
  dist <- 2 * atan2(sqrt(a), sqrt(1 - a)) * p[, 5]
  return(as.vector(dist))
}

long_lat <- coordinates[1:U,c("Long","Lat")]
dmat <- matrix(0,U,U)
for(u1 in 1:U) {
  for(u2 in 1:U) {
    dmat[u1,u2] <- round(distHaversine(long_lat[u1,],long_lat[u2,]) / 1609.344,1)
  }
}

p <- population0$Population[1:U]
v_by_g <- matrix(0,U,U)
dist_mean <- sum(dmat)/(U*(U-1))
p_mean <- mean(p)
for(u1 in 2:U){
  for(u2 in 1:(u1-1)){
    v_by_g[u1,u2] <- (dist_mean*p[u1]*p[u2]) / (dmat[u1,u2] * p_mean^2)
    v_by_g[u2,u1] <- v_by_g[u1,u2]
  }
}

v_by_g_day <- matrix(nrow = U*14, ncol = U)
for (i in 1:U) {
  for (j in 1:14) {
    v_by_g_day[14*i-14+j,] = v_by_g[i,]
  }
}

usethis::use_data(
  mobility,
  population,
  incidence,
  coordinates,
  population0,
  v_by_g_day,
  overwrite = TRUE
)
