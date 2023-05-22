#' Round values, returning a character string keeping trailing zeros
#'
#' Fixes the formatting issue that 'round' and 'signif' drop trailing zeros.
#'  
#' @param x a numeric vector
#' @param digits integer indicating the number of decimal places
#' ('myround') or significant digits ('signif') to be used.  For 'myround',
#' negative values are allowed
#'
#' @return a character vector
#'
#' @references 
#' 'myround' is adapted from 'broman::myround'. Similar functionality
#' is provided by 'finalfit::round_tidy'.
#'
#' @family utilities
#' @rdname myround
#' @examples
#' myround(1.3001,3)
#' mysignif(1.3001,3)
#'
#' @export
myround <- function (x, digits = 1) {
  if (length(digits) > 1) {
    digits <- digits[1]
    warning("Using only digits[1]")
  }
  if (digits < 1) {
    as.character(round(x,digits))
  } else {
    tmp <- sprintf(paste("%.", digits, "f", sep = ""), x)
    zero <- paste0("0.", paste(rep("0", digits), collapse = ""))
    tmp[tmp == paste0("-", zero)] <- zero
    tmp
  }
}
#' @rdname myround
#' @export
mysignif <- function (x, digits = 1) {
  myround(x, digits - ceiling(log10(abs(x))))
}
