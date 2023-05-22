#' R0 for li20 model
#'
#' Calculate  R_0 = alpha beta D + (1-alpha) mu beta D
#' for the model of Li et al
#'
#' @param params parameter vector
#' @param be logical vector to calculate using parameters before Jan 23
#' @return R0 value
#'
#' @examples
#' ## R0 for mle
#' po_mle <- li23(U=5,version="MLEperiod3")
#' p_mle <- coef(po_mle)
#' R0(p_mle,be=TRUE)
#' R0(p_mle,be=FALSE)
#' ## R0 for li20
#' po_li <- li23(U=5,version="li20period3")
#' p_li <- coef(po_li)
#' R0(p_li,be=TRUE)
#' R0(p_li,be=FALSE)
#'
#' @export
R0 <- function (params, be=T){
 a = params[ifelse(be,"alpha_be1","alpha_af1")]
 b = params[ifelse(be,"Beta_be1","Beta_af1")]
 D = params[ifelse(be,"D_be1","D_af1")]
 mu = params[ifelse(be,"mu_be1","mu_af1")]
 unname((a + (1-a) * mu) * b * D)
}
