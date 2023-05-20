
#' dag_test_FS
#'
#' Implements the fixed sequence testing procedure of familywise error rate control. The sequence is given through ordering elements of \code{X0} decreasingly.
#'
#' @param X0 a numeric vector or matrix with one column giving points corresponding to hypotheses.
#' @param p a numeric vector taking values in (0, 1] such that \code{length(p) == nrow(X0)} if X0 is a matrix (or \code{length(p) == length(X0)} if X0 is a numeric vector).
#' @param alpha a numeric value in (0, 1] specifying the Type I error rate.
#'
#' @return A boolean vector of the same length as \code{p} with each element being \code{TRUE} if the corresponding hypothesis is rejected and \code{FALSE} otherwise.
#' @export
#'
#' @examples
#' X0 <- c(0.5, 0, 1)
#' p <- c(0.01, 0.1, 0.05)
#' alpha <- 0.05
#' dag_test_FS(X0, p, alpha)
#'
dag_test_FS <- function(X0, p, alpha) {
  if (is.vector(X0) & is.atomic(X0)) {
    x0 <- X0
  } else if (is.atomic(X0) & ncol(X0) == 1) {
    x0 <- as.vector(X0)
  } else {
    stop("Data not obviously one-dimensional.")
  }

  n <- length(x0)

  if (all(p <= alpha)) {
    return(rep(TRUE, n))
  }

  # sort covariates decreasingly
  x0_order <- order(x0, decreasing = TRUE)

  # reject hypotheses corresponding to accordingly sorted p-values until first p-value is above alpha
  first_large_pvalue <- which.min(p[x0_order] <= alpha)

  if (first_large_pvalue == 1) {
    return(rep(FALSE, n))
  } else {
    R <- rep(FALSE, n)
    R[x0_order[1:(first_large_pvalue - 1)]] <- TRUE
    return(R)
  }
}
