#' Covid spatPomp generator for Jan 10 - 23
#'
#' Generate a \sQuote{spatPomp} object for Covid in \code{U} cities in China.
#'
#' @param li20_version choose which version of the model will be used
#' @param U A length-one numeric signifying the number of cities to be represented in the spatPomp object
#' @param dt step size, in days, for the euler approximation
#' @param sharedOneInterval estimated parameters that are equal for each unit and constant throughout time
#' @param sharedTwoInterval estimated parameters that are equal for each unit but change on day 15
#' @param version choose which set of parameters will be used for code testing
#' @param sharedParNames for li20v1, estimated parameters that are equal for each unit
#' @param unitParNames for li20v1, estimated parameters that are different for each unit
#' @param for_ibpf whether this object is generated for ibpf
#' @param mob_modify_factor A multiplicative factor, which is used to control the magnitude of adjustment for the mobility data
#' @import pomp
#' @import spatPomp
#' @return An object of class \sQuote{spatPomp} representing a \code{U}-dimensional spatially coupled Covid POMP model.
#'
#' @export

li20period1 <- function(li20_version = "v5", U = 373, dt = 1/4,
                        sharedOneInterval = c("theta","tau","sigma_SE","E_0","Iu_0"),
                        sharedTwoInterval = c("alpha","Beta","mu","Z","D","Td"),
                        version = c("MLEperiod3","li20period1","li20period2","li20period3"), 
                        sharedParNames = c("alpha_be","Beta_be","alpha_af","Beta_af",
                                           "mu_be","Z_be","D_be","mu_af","Z_af","D_af",
                                           "theta","tau","Td_be","Td_af","size_a",
                                           "size_b","E_0","Iu_0"),
                        unitParNames = NULL, for_ibpf = T, mob_modify_factor = 20
) {
  if (li20_version == "v5"){
    li20v5(U = U, dt = dt,
           sharedOneInterval = sharedOneInterval,
           sharedTwoInterval = sharedTwoInterval,
           version = version, mob_modify_factor=mob_modify_factor, days = 14, for_ibpf = for_ibpf)
  } else if (li20_version == "v3"){
    li20v3(U = U, dt = dt,
           sharedOneInterval = sharedOneInterval,
           sharedTwoInterval = sharedTwoInterval,
           version = version, mob_modify_factor=mob_modify_factor, days = 14, for_ibpf = for_ibpf)
  } else if (li20_version == "v1"){
    li20(U = U, dt = dt, 
         sharedParNames = sharedParNames,
         unitParNames = unitParNames, for_ibpf = for_ibpf, 
         mob_modify_factor=mob_modify_factor, days = 14)
  }
}


