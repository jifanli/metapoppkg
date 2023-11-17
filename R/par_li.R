#' Parameter transform for li20 and li23 model
#'
#' Transform li20 or li23 parameters to and from ranges used by Li et al (2020)
#'
#' @param pars names parameter vector
#' @param unitOneInterval subset of "theta","sigma_SE","tau","E_0","A_0" which are unit-specific
#' @param unitTwoInterval subset of "alpha","Beta","mu","Z","D","Td" which are unit-specific
#' @param U number of units
#' @param dir direction, "from_li" or "to_li"
#' @return transformed parameter vector
#'
#' @examples
#' par_li1 <- c(
#' alpha_be1 = 0.09389817,
#' Beta_be1  = 0.77943975,
#' mu_be1    = 0.96682356,
#' Z_be1     = 0.11650832,
#' D_be1     = 0.93467362,
#' Td_be1    = 0.80000000,
#' alpha_af1 = 0.46182906,
#' Beta_af1  = 0.13459503,
#' mu_af1    = 0.68640173,
#' Z_af1     = 0.04444878,
#' D_af1     = 0.12646338,
#' Td_af1    = 0.20000000,
#' theta1    = 0.89626861,
#' tau1      = 0.91718058,
#' sigma_SE1 = 0.13483710,
#' E_01      = 0.44362243,
#' A_01     = 0.98409838
#' )
#' par_nat <- par_li(par_li1,U=5,dir="from_li")
#' par_li2 <- par_li(par_nat,U=5,dir="to_li")
#' @export

par_li <- function(pars,unitOneInterval=NULL,unitTwoInterval=NULL,U,
  dir=c("from_li","to_li")){

  range <- data.frame(
    theta = c(1, 1.75),
    Beta = c(0.2, 1.5), # the lower bound of Beta_af is 0.2
    alpha = c(0.02, 0.98),
    Z = c(2, 5),
    D = c(2, 5),
    Td = c(5, 10),
    mu = c(0.2, 1),
    E_0 = c(0,2000),
    A_0 = c(0,2000),
    sigma_SE = c(0,5),
    tau = c(0,3),
    row.names=c("lo","hi")
  )
  liTransformed <- colnames(range)
  allOneInterval <- intersect(c("theta","sigma_SE","tau","E_0","A_0"),liTransformed)
  allTwoInterval <- intersect(c("alpha","Beta","mu","Z","D","Td"),liTransformed)

  sharedOneInterval <- setdiff(allOneInterval,unitOneInterval)
  sharedTwoInterval <- setdiff(allTwoInterval,unitTwoInterval)
  lo <- hi <- pars
  for(p in unitOneInterval){
    nn <- paste0(p,1:U)
    lo[nn] <- range["lo",p]
    hi[nn] <- range["hi",p]
  }
  for(p in unitTwoInterval){
    nn <- c(paste0(p,"_be",1:U),paste0(p,"_af",1:U))
    lo[nn] <- range["lo",p]
    hi[nn] <- range["hi",p]
  }
  if("Beta" %in% unitTwoInterval){
    nn <- paste0("Beta_be",1:U)
    lo[nn] <- 0.8
    hi[nn] <- 1.5
    nn <- paste0("Beta_af",1:U)
    lo[nn] <- 0.2
    hi[nn] <- 1.2
  }
  for(p in sharedOneInterval){
    nn <- paste0(p,1)
    lo[nn] <- range["lo",p]
    hi[nn] <- range["hi",p]
  }
  for(p in sharedTwoInterval){
    nn <- c(paste0(p,"_be",1),paste0(p,"_af",1))
    lo[nn] <- range["lo",p]
    hi[nn] <- range["hi",p]
  }
  if("Beta" %in% sharedTwoInterval){
    nn <- paste0("Beta_be",1)
    lo[nn] <- 0.8
    hi[nn] <- 1.5
    nn <- paste0("Beta_af",1)
    lo[nn] <- 0.2
    hi[nn] <- 1.2
  }
  if(dir[1]=="to_li"){
    if(any(pars<lo[names(pars)] | pars> hi[names(pars)])) stop(
     "one or more parameters out of range specified by Li et al (2020)")
    tol <- 1e-10
    pars[hi-lo>tol] = ((pars-lo)/(hi-lo))[hi-lo>tol]
    return(pars)
  } else if(dir[1]=="from_li"){
    return(lo + pars*(hi-lo))
  } else stop("dir should be 'to_li' or 'from_li'")
}

