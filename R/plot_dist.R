#' Mobility data plot
#'
#' Generate a plot of the mobility data.
#'
#' @param day indicate which day will be presented on the plot
#' @param U number of the cities
#' @param city_index indicate which cities will be presented, if U is not given
#' @param linewidth width of the lines
#' @param arrowsize size of the arrows
#' @param pointsize size of the points
#' @param mob_modify_factor A multiplicative factor, which is used to control the magnitude of adjustment for the mobility data
#' @importFrom ggplot2 ggplot geom_point geom_curve aes_string
#' @importFrom pomp melt
#' @return a plot of the mobility data
#'
#' @export


plot_dist <- function(day = 1, city_index=1:373, U=NULL, linewidth=1/50000, arrowsize=0.03, pointsize=1, mob_modify_factor=0) {
  mobi=metapoppkg::mobility+metapoppkg::v_by_g_day*mob_modify_factor
  coordinates=metapoppkg::coordinates
  day_index=1:373*14-14+day
  mobi=mobi[day_index,]

  if(!is.null(U)){
    # choose the cities with largest number of people who traveled from Wuhan
    index=order(mobi[,170], decreasing = T)[1:U-1]
    index=c(170,index) # 170 is Wuhan
  } else if(is.null(U)){
    index=city_index
    U=length(city_index)
  }
  
  mobi=mobi[index,index]
  rownames(mobi)=colnames(mobi)
  coordinates=coordinates[index,]
  mobi=pomp::melt(mobi)
  colnames(mobi)[1:2]=c("to","from")
  mobi$Long_from=rep(coordinates$Long,each=U)
  mobi$Lat_from=rep(coordinates$Lat,each=U)
  mobi$Long_to=rep(coordinates$Long,U)
  mobi$Lat_to=rep(coordinates$Lat,U)
  identical_index=which(mobi$to==mobi$from)
  mobi=mobi[-identical_index,]
  in_index=which(mobi$to=="Wuhan")
  out_index=which(mobi$from=="Wuhan")
  mobi$dir=rep("excluding Wuhan",nrow(mobi))
  mobi$dir[in_index]="toward Wuhan"
  mobi$dir[out_index]="away from Wuhan"

  ggplot2::ggplot(coordinates,aes_string("Long","Lat"))+
    ggplot2::geom_point(size=pointsize) + ggplot2::geom_curve(
      curvature = 0.1, 
      aes_string(x = "Long_from", y = "Lat_from", xend = "Long_to", yend = "Lat_to", col = "dir"),
      linewidth = mobi[-c(out_index,in_index),]$value*linewidth,
      data = mobi[-c(out_index,in_index),]
    ) + ggplot2::geom_curve(
      curvature = 0.1, 
      aes_string(x = "Long_from", y = "Lat_from", xend = "Long_to", yend = "Lat_to", col = "dir"),
      linewidth = mobi[out_index,]$value*linewidth,
      data = mobi[in_index,]
    ) + ggplot2::geom_curve(
      curvature = 0.1, 
      aes_string(x = "Long_from", y = "Lat_from", xend = "Long_to", yend = "Lat_to", col = "dir"),
      linewidth = mobi[out_index,]$value*linewidth,
      data = mobi[out_index,]
    )
    
}
