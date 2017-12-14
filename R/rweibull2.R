#' Generating observations from a reparameterized Weibull distribution.
#'
#' @param n       number of observations. If \code{length(n) > 1}, then the length is taken to be the number required.
#' @param alpha   a positive-valued scale parameter  (Default is 1)
#' @param lambda  a positive-valued shape parameter (Default is 1)
#' @details Given the scale parameter, \code{alpha}, and shape parameter, \code{lambda},
#' the density is given by 
#' \deqn{f(x|\alpha, \lambda) = \alpha\lambdax^(\alpha-1)exp(-\lambdax^\alpha),}
#' for \eqn{x > 0}. Note that this parameterization is different than the one provided in \code{rweibull}.
#' @examples
#' library(weibull)
#' x <- rweibull2(10000, alpha = 2, lambda = 5)
#' mean(x)
#' gamma(1 + 1 / 2) / 5^(1 / 2)
#' @export
#' @importFrom stats rweibull
#'
rweibull2 <- function(n, alpha = 1, lambda = 1) {
  if (alpha <= 0) {
    stop("Shape parameter must be positive.")
  } else if (lambda <= 0) {
    stop("Scale parameter must be positive.")
  }
  a <- alpha
  b <- lambda^(-1 / alpha)
  rweibull(n, shape = a, scale = b)
}



