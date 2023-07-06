
#' dag_test_FS
#'
#' Implements the fixed sequence testing procedure of familywise error rate control. The sequence is given through ordering elements of \code{p_order} increasingly.
#'
#' @param p_order a numeric vector or matrix with one column whose order determines the sequence of tests.
#' @param p a numeric vector taking values in (0, 1] such that \code{length(p) == nrow(p_order)} if p_order is a matrix (or \code{length(p) == length(p_order)} if p_order is a numeric vector).
#' @param alpha a numeric value in (0, 1] specifying the Type I error rate.
#' @param decreasing a boolean value determining whether the order of p_order should be understood in decreasing order.
#'
#' @return A boolean vector of the same length as \code{p} with each element being \code{TRUE} if the corresponding hypothesis is rejected and \code{FALSE} otherwise.
#' @export
#'
#' @examples
#' p_order <- c(0.5, 0, 1)
#' p <- c(0.01, 0.1, 0.05)
#' alpha <- 0.05
#' dag_test_FS(p_order, p, alpha, decreasing = TRUE)
#'
dag_test_FS <- function(p_order, p, alpha, decreasing = FALSE) {
  if (is.vector(p_order) & is.atomic(p_order)) {
    p_order_vec <- p_order
  } else if (is.atomic(p_order) & ncol(p_order) == 1) {
    p_order_vec <- as.vector(p_order)
  } else {
    stop("Data not obviously one-dimensional.")
  }

  n <- length(p_order_vec)

  if (all(p <= alpha)) {
    return(rep(TRUE, n))
  }

  # sort covariates decreasingly
  p_order_vec_order <- order(p_order_vec, decreasing = decreasing)

  # reject hypotheses corresponding to accordingly sorted p-values until first p-value is above alpha
  first_large_pvalue <- which.min(p[p_order_vec_order] <= alpha)

  if (first_large_pvalue == 1) {
    return(rep(FALSE, n))
  } else {
    R <- rep(FALSE, n)
    R[p_order_vec_order[1:(first_large_pvalue - 1)]] <- TRUE
    return(R)
  }
}
