#' Graphical Diagnostic Checks for Weibull AFT Model
#'
#' @param fit     object of class \code{survreg} and \code{distribution = "weibull"}. 
#' @return  
#' Returns the Cox-Snell residual plot and a plot of the deviance residuals against the study time.
#' Residuals (Cox-Snell, Martingale, and deviance) are also returned.
#' @examples
#' library(weibull)
#' library(survival)
#' library(KMsurv)
#' data(alloauto, packge = "KMsurv")
#' fit <- survreg(Surv(time, delta) ~ 1, data = alloauto, dist = "weibull", 
#' subset = (type == 1))
#' getWeibullDiagnostics(fit)
#' @export
#' @importFrom stats pnorm
#' @importFrom graphics plot lines abline
#' @importFrom survival survfit resid

getWeibullDiagnostics <- function(fit) {
  
  #- Make sure object is of class survreg and weibull model
  if (attr(fit, "class") != "survreg" |  fit$dist != "weibull") {
    stop("Object must be of class 'survreg' with distribution = 'weibull'.")
  }
  
  #lambda <- exp(-fit$coefficients[1] / fit$scale)
  alpha  <- 1 / fit$scale
  eta    <- -fit$linear.predictors * alpha # XBeta
  ti     <- exp(fit$y[, 1]) #Observed survival times
  di     <- fit$y[, 2] #Censoring Indicators
  
  #Cox-Snell Residuals
  ri <- ti^alpha * exp(eta)
  
  #Create Nelson-Aalen Estimator
  fit2 <- survfit(Surv(ri, di) ~ 1)
  H.na <- cumsum(fit2$n.event / fit2$n.risk) #Nelson Aalen Estimator
  
  plot(H.na ~ fit2$time, type = "l", main = "Cox-Snell Residual Plot",
       ylab = "Estimated Cumulative Hazard Rates", 
       xlab = "Cox-Snell Residual")
  abline(a = 0,  b = 1, col = "blue", lty = 2)
  
  #Deviance Residuals
  dev <- resid(fit, "deviance")
  plot(dev ~ ti, type = "p", main = "Deviance Residual Plot",
       ylab = "Deviance Residuals", 
       xlab = "Study Time", pch = 19, cex = 0.5)
  abline(h = 0, col = "blue", lty = 2)
  
  out <- list()
  out$coxsnell   <- ri
  out$martingale <- di - ri
  out$deviance   <- dev
  invisible(out)
}

