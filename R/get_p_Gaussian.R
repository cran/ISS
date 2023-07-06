
#' get_p_Gaussian
#'
#' Calculate the p-value in Definition 19 of \insertCite{MRCS2023;textual}{ISS}.
#'
#' @param X a numeric matrix specifying the covariates.
#' @param y a numeric vector with \code{length(y) == nrow(X)} specifying the
#' responses.
#' @param x0 a numeric vector specifying the point of interest, such that
#' \code{length(x0) == ncol(X)}.
#' @param tau a single numeric value specifying the threshold of interest.
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
#' y <- apply(X, MARGIN = 1, FUN = eta) + rnorm(n, sd = 1)
#' get_p_Gaussian(X, y, x0 = c(1, 1), tau = 1)
#' get_p_Gaussian(X, y, x0 = c(1, 1), tau = -1)
#'
#'

get_p_Gaussian <- function(X, y, x0, tau) {
  d <- ncol(X)

  # Reduce data to those observations with covariates coordinate-wise smaller than x0
  subset_ind <- colSums(t(X) <= x0) == d
  X_subset <- X[subset_ind, , drop = FALSE]
  y_subset <- y[subset_ind]
  nx <- length(y_subset)

  # Determine the values and ordering of the sup-norm for all covariates
  X_dist <- apply(abs(t(X_subset) - x0), MARGIN = 2, FUN = max)
  X_order <- order(X_dist)
  y_seq <- y_subset[X_order]

  # calculate sequence of MLE under the null hypothesis
  sigma_hat_sq_0 <- cumsum(pmax(y_seq - tau, 0)^2)/seq_along(y_seq)

  # calculate predictable sequence of unconstrained estimators
  # first vector entry corresponds to 0-index here - need to shift future
  # subsetting calls by +1
  y_bar_1 <- c(0, cumsum(y_seq))
  sigma_hat_sq_1 <- rep(NA, nx + 1)
  sigma_hat_sq_1[c(1,2)] <- 1
  for (k in 3:(nx+1)) {
    sigma_hat_sq_1[k] <- mean((y_seq[1:(k-1)] - y_bar_1[k])^2)
  }

  # calculate p_k vector
  p_sequence <- rep(NA, nx)
  for (k in 1:nx) {
    p_sequence[k] <- 1/(sigma_hat_sq_0[k]^(k/2) * exp(k/2)) *
      prod(sqrt(sigma_hat_sq_1[1:k]) *
             exp((y_seq[1:k] - y_bar_1[1:k])^2/(2 * sigma_hat_sq_1[1:k])))
  }

  # Then calculate the overall p-value
  p <- min(c(1, p_sequence))

  # return the p-value
  p
}
