
#' ISS
#'
#' The function implements the combination of p-value calculation and familywise
#' error rate control through DAG testing procedures described in \insertCite{MRCS2023;textual}{ISS}.
#'
#'
#' @param X a numeric matrix specifying the covariates.
#' @param y a numeric vector with \code{length(y) == nrow(X)} specifying the responses.
#' @param tau a single numeric value specifying the threshold of interest.
#' @param alpha a numeric value in (0, 1] specifying the Type I error rate.
#' @param m an integer value between 1 and \code{nrow(X)} specifying the size of
#' the subsample of \code{X} at which the hypotheses should be tested.
#' @param p_value one of \code{c("sub-Gaussian", "sub-Gaussian-normalmixture", "classification")} specifying
#' which p-value construction should be used. See Definitions 1, 18 and 19 by \insertCite{MRCS2023;textual}{ISS} respectively.
#' @param DAG_testing_procedure one of \code{c("ISS", "Holm", "MG all", "MG any")}, specifying how the
#' familywise error rate is controlled. The first corresponds to Algorithm 1 by \insertCite{MRCS2023;textual}{ISS},
#' the second is Holm's procedure and the final two correspond to the procedures due to \insertCite{meijer2015multiple;textual}{ISS}
#' for one-way logical relationships.
#' @param sigma2 a single positive numeric value specifying the variance parameter (only needed if \code{p_value == "sub-Gaussian"} or \code{p_value == "sub-Gaussian-normalmixture"}).
#' @param rho a single positive numeric value serving as hyperparameter (only used if \code{p_value == "sub-Gaussian-normalmixture"}).
#' @param minimal a logical value determining whether the output should be reduced to the minimal number of points leading to the same selected set.
#'
#' @return A numeric matrix giving the points in \code{X} determined to lie in the \code{tau}-superlevel set of the regression function with probability at least 1 - \code{alpha} or, if \code{minimal == TRUE}, a subset of points thereof that have the same upper hull.
#' @export
#'
#' @references
#' \insertRef{meijer2015multiple}{ISS}
#'
#' \insertRef{MRCS2023}{ISS}
#'
#' @examples
#' d <- 2
#' n <- 1000
#' m <- 100
#' sigma2 <- (1 / 4)^2
#' tau <- 0.5
#' alpha <- 0.05
#'
#' X <- matrix(runif(n * d), nrow = n)
#' eta_X <- apply(X, MARGIN = 1, max)
#' y <- eta_X + rnorm(n, sd = sqrt(sigma2))
#' X_rej <- ISS(X = X, y = y, tau = tau, alpha = alpha, m = m, sigma2 = sigma2)
#'
#' if (d == 2) {
#'   plot(0, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = NA, ylab = NA)
#'   for (i in 1:nrow(X_rej)) {
#'     rect(
#'       xleft = X_rej[i, 1], xright = 1, ybottom = X_rej[i, 2], ytop = 1,
#'       border = NA, col = "indianred"
#'     )
#'   }
#'
#'   points(X, pch = 16, cex = 0.5, col = "gray")
#'   points(X[1:m, ], pch = 16, cex = 0.5, col = "black")
#'   lines(x = c(0, tau), y = c(tau, tau), lty = 2)
#'   lines(x = c(tau, tau), y = c(tau, 0), lty = 2)
#'
#'   legend(
#'     x = "bottomleft",
#'     legend = c(
#'       "superlevel set boundary",
#'       "untested covariate points",
#'       "tested covariate points",
#'       "selected set"
#'     ),
#'     col = c("black", "gray", "black", "indianred"),
#'     lty = c(2, NA, NA, NA),
#'     lwd = c(1, NA, NA, NA),
#'     pch = c(NA, 16, 16, NA),
#'     fill = c(NA, NA, NA, "indianred"),
#'     border = c(NA, NA, NA, "indianred")
#'   )
#' }
#'
ISS <- function(X, y, tau, alpha = 0.05, m = nrow(X),
                p_value = c("sub-Gaussian-normalmixture", "sub-Gaussian", "classification"),
                sigma2, rho = 1 / 2,
                DAG_testing_procedure = c("ISS", "Holm", "MG all", "MG any"),
                minimal = FALSE) {
  n <- nrow(X)
  d <- ncol(X)

  p_value <- match.arg(p_value)
  DAG_testing_procedure <- match.arg(DAG_testing_procedure)

  Xm <- X[1:m, , drop = FALSE]

  # calculation of the p-value according to the chosen method
  p <- rep(NA, m)

  if (p_value == "sub-Gaussian") {
    for (i in 1:m) {
      x0 <- Xm[i, ]
      p[i] <- get_p_subGaussian(X = X, y = y, x0 = x0, sigma2 = sigma2, tau = tau)
    }
  } else if (p_value == "sub-Gaussian-normalmixture") {
    for (i in 1:m) {
      x0 <- Xm[i, ]
      p[i] <- get_p_subGaussian_NM(X = X, y = y, x0 = x0, sigma2 = sigma2, tau = tau, rho = rho)
    }
  } else if (p_value == "classification") {
    for (i in 1:m) {
      x0 <- Xm[i, ]
      p[i] <- get_p_classification(X = X, y = y, x0 = x0, tau = tau)
    }
  }

  # control FWER
  if (DAG_testing_procedure == "Holm") {
    rejection_set <- dag_test_Holm(X0 = Xm, p = p, alpha = alpha)
  } else if (DAG_testing_procedure == "MG all") {
    rejection_set <- dag_test_MG(X0 = Xm, p = p, alpha = alpha, version = "all", sparse = FALSE)
  } else if (DAG_testing_procedure == "MG any") {
    rejection_set <- dag_test_MG(X0 = Xm, p = p, alpha = alpha, version = "any", sparse = FALSE)
  } else if (DAG_testing_procedure == "ISS") {
    rejection_set <- dag_test_ISS(X0 = Xm, p = p, alpha = alpha)
  }

  # return the points whose upper hull is the rejection set
  # (i.e. the points define the boundary of the selected subgroup)
  rejected_points <- Xm[rejection_set, , drop = FALSE]
  if (minimal) {
    get_boundary_points(rejected_points)
  } else {
    rejected_points
  }
}
