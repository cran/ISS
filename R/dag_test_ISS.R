
#' dag_test_ISS
#'
#' Implements the DAG testing procedure given in Algorithm 1 by \insertCite{MRCS2023;textual}{ISS}.
#'
#' @param X0 a numeric matrix giving points corresponding to hypotheses.
#' @param p a numeric vector taking values in (0, 1] such that \code{length(p) == nrow(X0)}.
#' @param alpha a numeric value in (0, 1] specifying the Type I error rate.
#'
#' @return A boolean vector of the same length as \code{p} with each element being \code{TRUE} if the corresponding hypothesis is rejected and \code{FALSE} otherwise.
#' @export
#'
#' @references \insertRef{MRCS2023}{ISS}
#'
#' @examples
#' X0 <- rbind(c(0.5, 0.6), c(0.8, 0.9), c(0.9, 0.8))
#' p <- c(0.02, 0.025, 0.1)
#' alpha <- 0.05
#' dag_test_ISS(X0, p, alpha)
#'
dag_test_ISS <- function(X0, p, alpha) {

  # Get  number of hypotheses
  n <- nrow(X0)
  d <- ncol(X0)

  # The procedure reduces to a much simpler procedure if d == 1.
  # The next if-else-branching is thus only for computational efficiency, could
  # just as well only look at the else-branch.

  if (d == 1) {
    dag_test_FS(p_order = X0, p = p, alpha = alpha, decreasing = TRUE)
  } else {
    iG <- get_DAG(X0, sparse = FALSE)
    iF <- get_DAG(X0, sparse = TRUE, twoway = TRUE)

    an_G <- iG$ancestors
    pa_F <- iF$parents
    de_F <- iF$descendants
    topo_F <- iF$topological_ordering
    leaves_F <- iF$leaves

    # Initialise rejected hypotheses (none to start with).
    R <- rep(FALSE, n)

    # Start the iterative procedure.
    # This will take AT MOST n steps (i.e. one hypothesis rejected at each iteration).
    for (j in 1:n) {
      # Setup
      R_rej <- which(R)
      R_unrej <- which(!R)
      # alpha_vec <- rep(0, n) # For budget on each node.
      unrej_leaves <- intersect(leaves_F, R_unrej) # Leaves that have not been rejected yet.
      candidate_bool <- rep(NA, n)
      alpha_vec <- rep(0, n)
      for (i in 1:n) {
        candidate_bool[i] <- (i %in% R_unrej) & all(pa_F[[i]] %in% R_rej)
        if (candidate_bool[i]) {
          alpha_vec[i] <- length(intersect(c(i, de_F[[i]]), unrej_leaves)) * alpha / length(unrej_leaves)
        }
      }

      if (all(p > alpha_vec)) break

      new_rejections <- which(p <= alpha_vec)

      R[new_rejections] <- TRUE
      for (i in new_rejections) R[an_G[[i]]] <- TRUE

      if (all(R)) break
    }

    R
  }
}
