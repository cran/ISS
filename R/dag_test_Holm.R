
#' dag_test_Holm
#'
#' Given a vector of p-values, each concerning a row in the matrix X0,
#' \code{dag_test_Holm()} first applies Holm's method to the p-values and then also rejects
#' hypotheses corresponding to points coordinate-wise greater or equal to any
#' point whose hypothesis has been rejected.
#'
#' @param X0 a numeric matrix giving points corresponding to hypotheses.
#' @param p a numeric vector taking values in (0, 1] such that \code{length(p) == nrow(X0)}.
#' @param alpha a numeric value in (0, 1] specifying the Type I error rate.
#'
#' @return A boolean vector of the same length as \code{p} with each element being \code{TRUE} if the corresponding hypothesis is rejected and \code{FALSE} otherwise.
#' @export
#'
#' @examples
#' X0 <- rbind(c(0.5, 0.5), c(0.8, 0.9), c(0.4, 0.6))
#' p <- c(0.01, 0.1, 0.05)
#' alpha <- 0.05
#' dag_test_Holm(X0, p, alpha)
#'
dag_test_Holm <- function(X0, p, alpha) {
  # check which p-values are smaller than the curved boundary that Holm's method induces
  n <- length(p)
  d <- ncol(X0)
  p_small <- (sort(p) <= alpha / (n - seq_along(p) + 1))

  # determine the p-values that lead to a rejection according to Holm's method
  if (all(p_small)) {
    return(rep(TRUE, n))
  }
  if (all(!p_small)) {
    return(rep(FALSE, n))
  }

  reject_number <- which.min(p_small) - 1
  reject_ind <- order(p)[1:reject_number]

  # get the hypotheses that are necessarily false if the previous rejections were appropriate
  rej_count <- rep(0, n)

  for (i in 1:length(reject_ind)) {
    x_i <- X0[reject_ind[i], ] # i-th rejected point
    rej_count <- rej_count + (colSums(t(X0) >= x_i) == d)
  }

  # return boolean concerning rejection (TRUE = reject null, FALSE = retain null)
  rej_count > 0
}
