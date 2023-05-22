#' Plotting li20 data
#'
#' Visualize li20 data, set the higher case numbers to black and the lower ones to white, 
#' sort the units by population size. And set the color of Wuhan to blue, Shenyang to red.
#'
#' @param U number of the cities
#' @param x a \code{spatPomp} object
#' @param limit the upper limit of case numbers in the plot
#' @param order the order of the units in the plot
#' @import ggplot2
#' @import ggnewscale
#' @importFrom rlang sym
#' @return a \code{ggplot} plot of class \sQuote{gg} and \sQuote{ggplot} visualizing
#' the time series data over multiple spatial units via a tile-plot.
#'
#' @export


plot_li20 <- function(x, U=373, limit=4, order="population") {
  df <- as.data.frame(x)
  df[x@unit_obsnames] <- log10(df[x@unit_obsnames]+1)
  # change the order, sort by population
  if (order=="population"){
    popu <- df$pop[1:U]
    index <- order(order(popu, decreasing = T))
  } else if (order=="cases"){
    index <- 1:U
  }
  wuhan_index <- which(df$city=="Wuhan")
  df$Wuhan <- NA
  df$Wuhan[wuhan_index] <- df$cases[wuhan_index]
  shenyang_index <- which(df$city=="Shenyang")
  df$Shenyang <- NA
  df$Shenyang[shenyang_index] <- df$cases[shenyang_index]
  df$y_ticks <- rep(index,30)

  ggplot2::ggplot(data = df,
                  mapping = ggplot2::aes(x = !!rlang::sym(x@timename),
                                         y = .data$y_ticks)) +
    ggplot2::scale_x_continuous(expand=c(0,0)) +
    ggplot2::scale_y_continuous(expand=c(0,0), sec.axis = sec_axis(~ .)) +
    ggplot2::geom_tile(mapping = ggplot2::aes(fill = !!rlang::sym(x@unit_obsnames))) +
    ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#000000",limits=c(0,limit)) +
    ggplot2::labs(x = "time",
                  y = "unit",
                  fill = paste("log10\n(",x@unit_obsnames,"+1)",sep="")) +
    ggnewscale::new_scale("fill") +
    ggplot2::geom_tile(mapping = ggplot2::aes(fill = .data$Wuhan), data = ~subset(., !is.na(Wuhan))) +
    ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#0000FF",limits=c(0,limit)) + 
    ggplot2::labs(fill = "log10\n(Wuhan+1)") +
    ggnewscale::new_scale("fill") +
    ggplot2::geom_tile(mapping = ggplot2::aes(fill = .data$Shenyang), data = ~subset(., !is.na(Shenyang))) +
    ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#FF0000",limits=c(0,limit)) + 
    ggplot2::labs(fill = "log10\n(Shenyang+1)") +
    ggplot2::theme(
      axis.text.y=ggplot2::element_blank(),
      axis.ticks.y.left = element_blank(),
      panel.border=ggplot2::element_blank())
  
}
