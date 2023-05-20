
#' get_p_subGaussian
#'
#' Calculate the p-value in Definition 1 of \insertCite{MRCS2023;textual}{ISS}.
#'
#' @param X a numeric matrix specifying the covariates.
#' @param y a numeric vector with \code{length(y) == nrow(X)} specifying the
#' responses.
#' @param x0 a numeric vector specifying the point of interest, such that
#' \code{length(x0) == ncol(X)}.
#' @param sigma2 a single positive numeric value specifying the variance parameter.
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
#' X <- matrix(runif(d*n), ncol = d)
#' eta <- function(x) sum(x)
#' y <- apply(X, MARGIN = 1, FUN = eta) + rnorm(n, sd = 0.5)
#' get_p_subGaussian(X, y, x0 = c(1, 1), sigma2 = 0.25, tau = 1)
#' get_p_subGaussian(X, y, x0 = c(1, 1), sigma2 = 0.25, tau = 3)


get_p_subGaussian <- function(X, y, x0, sigma2, tau) {

  d <- ncol(X)

  # Reduce data to those observations with covariates coordinate-wise at most x0
  subset_ind <- colSums(t(X) <= x0) == d
  X_subset <- X[subset_ind, , drop = FALSE]
  y_subset <- y[subset_ind]

  # Determine the values and ordering of the sup-norm for
  X_dist <- apply(abs(t(X_subset) - x0), MARGIN = 2, FUN = max)
  X_order <- order(X_dist)

  # Calculate the martingale test sequence
  S <- cumsum((y_subset[X_order] - tau) / sqrt(sigma2))

  # Convert the martingale test sequence to p-values
  S_plus <- pmax(S, 0)
  p_sequence <- 5.2 * exp(-S_plus^2 / (2.0808 * seq_along(S_plus)) +
    log(log(2 * seq_along(S_plus))) / 0.72)

  # Then calculate the overall p-value
  p <- min(c(1, p_sequence))

  # return the p-value
  p
}


