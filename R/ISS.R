
#' ISS
#'
#' The function implements the combination of p-value calculation and familywise
#' error rate control through DAG testing procedures described in \insertCite{MRCS2023v2;textual}{ISS}.
#'
#'
#' @param X a numeric matrix specifying the covariates.
#' @param y a numeric vector with \code{length(y) == nrow(X)} specifying the responses.
#' @param tau a single numeric value specifying the threshold of interest.
#' @param alpha a numeric value in (0, 1] specifying the Type I error rate.
#' @param m an integer value between 1 and \code{nrow(X)} specifying the size of
#' the subsample of \code{X} at which the hypotheses should be tested.
#' @param p_value one of \code{c("sub-Gaussian", "sub-Gaussian-normalmixture", "Gaussian", "classification", "quantile")} specifying
#' which p-value construction should be used. See Definitions 1, 18, 19 and 21 and Lemma 24 by \insertCite{MRCS2023v2;textual}{ISS} respectively.
#' For \code{p_value == "quantile"}, the version with the p-value from Definition 19 is implemented.
#' @param FWER_control one of \code{c("ISS", "Holm", "MG all", "MG any", "split", "split oracle")}, specifying how the
#' familywise error rate is controlled. The first corresponds to Algorithm 1 by \insertCite{MRCS2023v2;textual}{ISS},
#' the second is Holm's procedure, the two starting with "MG" correspond to the procedures by \insertCite{meijer2015multiple;textual}{ISS}
#' for one-way logical relationships, and the final two containing "split" to the sample splitting techniques in Appendix B of \insertCite{MRCS2023v2;textual}{ISS}.
#' @param sigma2 a single positive numeric value specifying the variance parameter (only needed if \code{p_value \%in\% c("sub-Gaussian", "sub-Gaussian-normalmixture")}).
#' @param rho a single positive numeric value serving as hyperparameter (only used if \code{p_value == "sub-Gaussian-normalmixture"}).
#' @param theta a single numeric value in (0, 1) specifying the quantile of interest when \code{p_value_method == "quantile"}. Defaults to 1/2, i.e.~the median.
#' @param minimal a logical value determining whether the output should be reduced to the minimal number of points leading to the same selected set.
#' @param split_proportion when \code{FWER_control \%in\% c("split", "split oracle")}, the number of data points in the first split of the data is \code{ceiling(split_proportion * nrow(X))}.
#' @param eta when \code{FWER_control == "split oracle"}, this parameter needs to be used to provide the true regression function, which should take a vector of covariates as inputs and output a single numeric value.
#'
#' @return A numeric matrix giving the points in \code{X} determined to lie in the \code{tau}-superlevel set of the regression function with probability at least 1 - \code{alpha} or, if \code{minimal == TRUE}, a subset of points thereof that have the same upper hull.
#' @export
#'
#' @importFrom Rdpack reprompt
#'
#' @references
#' \insertAllCited{}
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
                p_value = c("sub-Gaussian-normalmixture", "sub-Gaussian", "Gaussian", "classification", "quantile"),
                sigma2, rho = 1 / 2,
                FWER_control = c("ISS", "Holm", "MG all", "MG any", "split", "split oracle"),
                minimal = FALSE, split_proportion = 1/2, eta = NA, theta = 1/2) {
  n <- nrow(X)
  d <- ncol(X)

  p_value <- match.arg(p_value)
  FWER_control <- match.arg(FWER_control)

  Xm <- X[1:m, , drop = FALSE]
  is_sample_split_based <- FWER_control %in% c("split", "split oracle")

  if (!is_sample_split_based) {

    # calculation of the p-value according to the chosen method
    p <- rep(NA, m)
    for (i in 1:m) {
      x0 <- Xm[i, ]
      p[i] <- get_p_value(p_value_method = p_value,
                          X = X, y = y, x0 = x0,
                          tau = tau, sigma2 = sigma2, rho = rho, theta = theta)
    }

    # control FWER
    if (FWER_control == "Holm") {
      rejection_set <- dag_test_Holm(X0 = Xm, p = p, alpha = alpha)
    } else if (FWER_control == "MG all") {
      rejection_set <- dag_test_MG(X0 = Xm, p = p, alpha = alpha, version = "all", sparse = FALSE)
    } else if (FWER_control == "MG any") {
      rejection_set <- dag_test_MG(X0 = Xm, p = p, alpha = alpha, version = "any", sparse = FALSE)
    } else if (FWER_control == "ISS") {
      rejection_set <- dag_test_ISS(X0 = Xm, p = p, alpha = alpha)
    }

  } else {

    # split the data
    D0 <- sort(sample(1:n, size = ceiling(n * split_proportion)))
    D1 <- (1:n)[-D0]

    # determine ordering of tests (larger value indicates later position in sequence)
    if (FWER_control == "split") {
      # calculate p-values based on first fold of data to determine order of hypotheses
      # for fixed sequence testing
      p_order <- rep(NA, m)
      for (i in 1:m) {
        x0 <- Xm[i, ]
        p_order[i] <- get_p_value(p_value_method = p_value,
                            X = X[D0, , drop = FALSE], y = y[D0], x0 = x0,
                            sigma2 = sigma2, tau = tau, rho = rho, theta = theta)
      }
    } else if (FWER_control == "split oracle") {
      # calculate value of eta at the m points and transform to values in (0, 1]
      # to determine order of hypotheses for fixed sequence testing. Don't use
      # first fold of data.
      Xm_eta <- apply(Xm, MARGIN = 1, FUN = eta)
      p_order <- rank(-Xm_eta)/length(Xm_eta)
    }

    # re-calculate p-values based on second fold of the data
    p <- rep(NA, m)
    for (i in 1:m) {
      x0 <- Xm[i, ]
      p[i] <- get_p_value(p_value_method = p_value,
                          X = X[D1, , drop = FALSE], y = y[D1], x0 = x0,
                          tau = tau, sigma2 = sigma2, rho = rho, theta = theta)
    }

    # apply fixed sequence testing
    rejection_set <- dag_test_FS(p_order = p_order, p = p, alpha = alpha)
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
