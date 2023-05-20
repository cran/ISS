#' get_p_classification
#'
#' Calculate the p-value in Definition 19 of \insertCite{MRCS2023;textual}{ISS}.
#'
#' @param X a numeric matrix specifying the covariates.
#' @param y a numeric vector with \code{length(y) == nrow(X)} and \code{all((y >= 0) & (y <= 1))} specifying the responses.
#' @param x0 a numeric vector specifying the point of interest, such that
#' \code{length(x0) == ncol(X)}.
#' @param tau a single numeric value in [0,1) specifying the threshold of interest.
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
#' get_p_classification(X, y, x0 = c(1, 1), tau = 0.6)
#' get_p_classification(X, y, x0 = c(1, 1), tau = 0.9)
get_p_classification <- function(X, y, x0, tau) {
  if (!all((y >= 0) & (y <= 1))) stop("All elements of the response vector y need to take values in [0, 1].")

  d <- ncol(X)

  # Reduce data to those observations with covariates coordinate-wise smaller than x0
  subset_ind <- colSums(t(X) <= x0) == d
  X_subset <- X[subset_ind, , drop = FALSE]
  y_subset <- y[subset_ind]

  # Determine the values and ordering of the sup-norm for all covariates
  X_dist <- apply(abs(t(X_subset) - x0), MARGIN = 2, FUN = max)
  X_order <- order(X_dist)


  # Calculate the martingale test sequence
  S <- cumsum(y_subset[X_order])
  nvec <- seq_along(S)

  # Convert the martingale test sequence to p-values
  # incomplete beta function of x with parameters a and b is simply equal
  # pbeta(x, a, b) * beta(a, b)
  p_sequence <- tau^S * (1 - tau)^(nvec - S + 1) /
    (stats::pbeta(1 - tau, nvec - S + 1, S + 1) * beta(nvec - S + 1, S + 1))

  # Then calculate the overall p-value
  p <- min(c(1, p_sequence[!is.nan(p_sequence)]))

  # return the p-value
  p
}
