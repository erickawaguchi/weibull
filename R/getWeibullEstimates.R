#' Tranforming Log-Linear Estimates to Weibull Estimates
#'
#' @param fit     object of class \code{survreg} and \code{distribution = "weibull"}. 
#' @return Returns ANOVA tables for both the Weibull representation, \code{WeibullModel}, and log-linear
#' representation, \code{LogLinearModel}. Note, the log-linear estimates are the same as
#' returned from \code{survreg}. The variance-covariance matrix for the Weibull estimates
#' can be extracted from \code{var}.
#' @examples
#' library(weibull)
#' library(survival)
#' library(KMsurv)
#' data(alloauto, packge = "KMsurv")
#' fit <- survreg(Surv(time, delta) ~ 1, data = alloauto, dist = "weibull", 
#' subset = (type == 1))
#' getWeibullEstimates(fit)
#' @export
#' @importFrom stats pnorm
#' 
getWeibullEstimates <- function(fit) {
  
  #- Make sure object is of class survreg and weibull model
  if (attr(fit, "class") != "survreg" |  fit$dist != "weibull") {
    stop("Object must be of class 'survreg' with distribution = 'weibull'.")
  }
  
  n.est <- length(fit$coefficients) + 1
  mu    <- fit$coefficients[1]
  sigma <- fit$scale
  Var   <- fit$var
  
  weibull.est <- c(exp(-mu / sigma), 1 / sigma)
  
  param <- c("lambda", "alpha")
  if(n.est > 2) { #If model has covariates
    param <- c(param , names(fit$coefficients)[-1])
    g.est <- fit$coefficients[-1]
    weibull.est <- c(weibull.est, -g.est / sigma)
  }
  
  names(weibull.est) <- param
  
  #- Get variance estimates (via delta method)
  g1 <- c(-exp(- 1 / sigma * mu - log(sigma)), rep(0, n.est - 2),
          mu * exp(- 1 / sigma * mu - log(sigma)))
  g2 <- c(rep(0, n.est - 1), - 1 / sigma)
  g  <- rbind(g1, g2)
  
  if (n.est > 2) {
    g3 <- cbind(rep(0, n.est - 2), diag(- 1 / sigma, n.est - 2), g.est / sigma)
    g  <- rbind(g, g3)
  }
  var.est <- g %*% Var %*% t(g)
  colnames(var.est) <- rownames(var.est) <- param
  
  se <- sqrt(diag(var.est))
  
  if (n.est > 2) {
    z  <- (weibull.est - c(0, 1, rep(0, n.est - 2))) / se
    logHR <- weibull.est[3:n.est]
  } else {
    z  <- (weibull.est - c(0, 1)) / se
  }
  out <- list()
  
  out$WeibullModel <- round(data.frame(Value = weibull.est, se = se, 
                                       z = z, p = 2 * (1 - abs(pnorm(z)))), 3)
  if (n.est > 2) {
    out$estimates <- round(data.frame(logHR = logHR, 
                                      HR = exp(logHR)), 3)
  }
  out$var <- var.est
  
  out$LogLinearModel <- summary(fit)$table
  return(out)
}
