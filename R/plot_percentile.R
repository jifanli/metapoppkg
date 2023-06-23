#' Make 10, 50, 90 percentile plots for li20 or li23
#'
#' This function makes spatiotemporal plots, with a 3-plot grid giving the 10, 50, 90 percentiles.
#'
#' @param U number of the cities
#' @param Nsim number of simulations
#' @param order the order of the units in the plot
#' @param model choose which model will be used
#' @param seed the seed used for simulations
#' @param registerDoRNG the seed used to register the doRNG foreach backend
#' @param params a named numeric vector or a matrix with rownames
#' @importFrom ggplot2 ggplot aes sec_axis element_blank
#' @importFrom ggnewscale new_scale
#' @importFrom doParallel registerDoParallel
#' @importFrom doRNG registerDoRNG
#' @importFrom dplyr .data
#' @import patchwork
#' @import spatPomp
#' @importFrom rlang sym
#' @importFrom stats quantile
#' @importFrom foreach %dopar%
#' @return a \code{ggplot} plot of class \sQuote{gg} and \sQuote{ggplot} visualizing
#' the time series data over multiple spatial units via a tile-plot.
#'
#' @export

plot_percentile <- function(U = 373, Nsim=100, order="population", model = "li23", seed = 12315, registerDoRNG = 3123465, params) {
  cores <-  as.numeric(Sys.getenv('SLURM_NTASKS_PER_NODE', unset = NA))
  if(is.na(cores)) cores <- parallel::detectCores()
  doParallel::registerDoParallel(cores)
  
  if (model == "li23"){
    x <- li23(U=U)
  } else if (model == "li20"){
    x <- li20(U=U)
  }
  
  set.seed(seed)
  doRNG::registerDoRNG(registerDoRNG)
  sim_out <- foreach::foreach(i = 1:Nsim, .packages = 'spatPomp') %dopar% obs(simulate(x,params=params))
  
  sim_array <- array(unlist(sim_out), c(U, 30, length(sim_out)))
  quantiles1 <- log10(apply(sim_array, c(1, 2), function(x) quantile(x, probs = 0.1))+1)
  quantiles2 <- log10(apply(sim_array, c(1, 2), function(x) quantile(x, probs = 0.5))+1)
  quantiles3 <- log10(apply(sim_array, c(1, 2), function(x) quantile(x, probs = 0.9))+1)
  
  quantiles1 <- data.frame(value = matrix(as.vector(quantiles1), nrow = U*30, ncol = 1))
  quantiles2 <- data.frame(value = matrix(as.vector(quantiles2), nrow = U*30, ncol = 1))
  quantiles3 <- data.frame(value = matrix(as.vector(quantiles3), nrow = U*30, ncol = 1))
  limit <- max(quantiles1,quantiles2,quantiles3)
  
  df <- as.data.frame(x)
  # change the order, sort by population
  if (order=="population"){
    popu <- df$pop[1:U]
    index <- order(order(popu, decreasing = T))
  } else if (order=="cases"){
    index <- 1:U
  }
  df$y_ticks <- rep(index,30)
  
  df$cases <- quantiles1$value
  wuhan_index <- which(df$city=="Wuhan")
  df$Wuhan <- NA
  df$Wuhan[wuhan_index] <- df$cases[wuhan_index]
  shenyang_index <- which(df$city=="Shenyang")
  df$Shenyang <- NA
  df$Shenyang[shenyang_index] <- df$cases[shenyang_index]
  
  plot1 <- ggplot2::ggplot(data = df,
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
      panel.border=ggplot2::element_blank(),
      legend.position = "none",
      axis.title.x=ggplot2::element_blank(), 
      axis.ticks.x=ggplot2::element_blank(),
      axis.ticks.y=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_blank(),
      axis.line.y = ggplot2::element_line(linewidth = 0.5, linetype = "solid")
    ) + 
    ggplot2::ggtitle("10th percentile")
  
  df$cases <- quantiles2$value
  df$Wuhan[wuhan_index] <- df$cases[wuhan_index]
  df$Shenyang[shenyang_index] <- df$cases[shenyang_index]
  
  plot2 <- ggplot2::ggplot(data = df,
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
       panel.border=ggplot2::element_blank(),
      legend.position = "none",
      axis.title.x=ggplot2::element_blank(), 
      axis.ticks.x=ggplot2::element_blank(), 
      axis.ticks.y=ggplot2::element_blank(),
      axis.text.x=ggplot2::element_blank(),
      axis.line.y = ggplot2::element_line(linewidth = 0.5, linetype = "solid")
   ) + 
    ggplot2::ggtitle("50th percentile")
  
  df$cases <- quantiles3$value
  df$Wuhan[wuhan_index] <- df$cases[wuhan_index]
  df$Shenyang[shenyang_index] <- df$cases[shenyang_index]
  
  plot3 <- ggplot2::ggplot(data = df,
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
      axis.ticks.y=ggplot2::element_blank(),
      panel.border=ggplot2::element_blank(),
      legend.position = "bottom",
      axis.line.y = ggplot2::element_line(linewidth = 0.5, linetype = "solid")
    ) + 
    ggplot2::ggtitle("90th percentile")

  
  panel_plot=plot1/plot2/plot3
  panel_plot
  
}
