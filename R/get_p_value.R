

#' get_p_value
#'
#' A wrapper function used to call the correct function for calculating the p-value.
#'
#' @param p_value_method one of \code{c("sub-Gaussian", "sub-Gaussian-normalmixture", "Gaussian", "classification", "quantile")} specifying
#' which p-value construction should be used. See Definitions 1, 18, 19 and 21 and Lemma 24 by \insertCite{MRCS2023;textual}{ISS} respectively.
#' For \code{p_value_method == "quantile"}, the version with the p-value from Definition 19 is implemented.
#' @param X a numeric matrix specifying the covariates.
#' @param y a numeric vector with \code{length(y) == nrow(X)} specifying the responses.
#' @param x0 a numeric vector specifying the point of interest, such that \code{length(x0) == ncol(X)}.
#' @param tau a single numeric value specifying the threshold of interest.
#' @param sigma2 a single positive numeric value specifying the variance parameter (required only if \code{p_value_method \%in\% c("sub-Gaussian", "sub-Gaussian-normalmixture"}).
#' @param rho a single positive numeric value serving as hyperparameter (required only if \code{p_value_method == "sub-Gaussian-normalmixture"}).
#' @param theta a single numeric value in (0, 1) specifying the quantile of interest when \code{p_value_method == "quantile"}. Defaults to 1/2, i.e.~the median.
#'
#' @return A single numeric value in (0, 1].
#' @export
#'
#' @references \insertRef{MRCS2023}{ISS}
#'
#' @examples
#' set.seed(123)
#' n <- 100
#' d <- 2
#' X <- matrix(runif(d * n), ncol = d)
#' eta <- function(x) sum(x)
#' X_eta <- apply(X, MARGIN = 1, FUN = function(x) 1 / (1 + exp(-eta(x))))
#' y <- as.numeric(runif(n) < X_eta)
#' get_p_value(p_value_method = "classification", X, y, x0 = c(1, 1), tau = 0.6)
#' get_p_value(p_value_method = "classification", X, y, x0 = c(1, 1), tau = 0.9)
#'
#' X_eta <- apply(X, MARGIN = 1, FUN = eta)
#' y <- X_eta + rcauchy(n)
#' get_p_value(p_value_method = "quantile", X, y, x0 = c(1, 1), tau = 1/2)
#' get_p_value(p_value_method = "quantile", X, y, x0 = c(1, 1), tau = 3)
#' get_p_value(p_value_method = "quantile", X, y, x0 = c(1, 1), tau = 3, theta = 0.95)
#'
get_p_value <- function(p_value_method = c("sub-Gaussian-normalmixture", "sub-Gaussian", "Gaussian", "classification", "quantile"),
                        X, y, x0, tau, sigma2, rho = 1/2, theta = 1/2) {

  p_value_method <- match.arg(p_value_method)

  if (p_value_method == "sub-Gaussian") {
      p <- get_p_subGaussian(X = X, y = y, x0 = x0, sigma2 = sigma2, tau = tau)
  } else if (p_value_method == "sub-Gaussian-normalmixture") {
      p <- get_p_subGaussian_NM(X = X, y = y, x0 = x0, sigma2 = sigma2, tau = tau, rho = rho)
  } else if (p_value_method == "classification") {
      p <- get_p_classification(X = X, y = y, x0 = x0, tau = tau)
  } else if (p_value_method == "Gaussian") {
      p <- get_p_Gaussian(X = X, y = y, x0 = x0, tau = tau)
  } else if (p_value_method == "quantile") {
      y_transformed <- as.numeric(y > tau)
      p <- get_p_classification(X = X, y = y_transformed, x0 = x0, tau = 1 - theta)
  }

  p
}


