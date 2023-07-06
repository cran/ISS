
#' dag_test_MG
#'
#' Implements the graph-testing procedures proposed by
#' \insertCite{meijer2015multiple;textual}{ISS} for one-way logical relationships.
#' Here implemented for the specific application to isotonic subgroup selection.
#'
#' @param X0 a numeric matrix giving points corresponding to hypotheses.
#' @param p a numeric vector taking values in (0, 1] such that \code{length(p) == nrow(X0)}.
#' @param alpha a numeric value in (0, 1] specifying the Type I error rate.
#' @param version either \code{"all"}
#' for the all-parent version of the procedure or \code{"any"} for the any-parent version of the procedure.
#' @param leaf_weights optional weights for the leaf nodes. Would have to be a numeric vector
#' of the same length as there are leaf nodes in the DAG (resp. polytree, see \code{sparse}) induced by \code{X0}.
#' @param sparse a logical value specifying whether \code{X0} should be used to
#' induce a DAG (\code{FALSE}) or a polytree (\code{TRUE}).
#'
#'
#' @return A boolean vector of the same length as \code{p} with each element being \code{TRUE} if the corresponding hypothesis is rejected and \code{FALSE} otherwise.
#' @export
#'
#' @references \insertRef{meijer2015multiple}{ISS}
#'
#' @examples
#' X0 <- rbind(c(0.5, 0.6), c(0.8, 0.9), c(0.9, 0.8))
#' p <- c(0.02, 0.025, 0.1)
#' alpha <- 0.05
#' dag_test_MG(X0, p, alpha)
#' dag_test_MG(X0, p, alpha, version = "any")
#' dag_test_MG(X0, p, alpha, sparse = TRUE)
#'
dag_test_MG <- function(X0, p, alpha, version = c("all", "any"), leaf_weights, sparse = FALSE) {
  version <- match.arg(version)

  # Get  number of hypotheses
  n <- nrow(X0)
  d <- ncol(X0)

  # Both procedurs by MG reduce to a much simpler procedure if d == 1.
  # The next if-else-branching is thus only for computational efficiency, could
  # just as well only look at the else-branch.

  if (d == 1) {
    dag_test_FS(p_order = X0, p = p, alpha = alpha, decreasing = TRUE)
  } else {

    # Calculate the ancestors and parents of each hypothesis in the induced graph and
    # calculate leaves, roots and topological ordering of the same graph.
    G <- get_DAG(X0, sparse = sparse)
    an <- G$ancestors
    pa <- G$parents
    leaves <- G$leaves
    topo <- G$topological_ordering
    roots <- which(sapply(pa, FUN = length) == 0)

    # Define default leaves_weights if none specified
    if (missing(leaf_weights)) leaf_weights <- rep(1, length(leaves))
    leaves_weights_full <- rep(NA, n)
    leaves_weights_full[leaves] <- leaf_weights

    # Initialise rejected hypotheses (none to start with).
    R <- rep(FALSE, n)

    # Start the iterative procedure.
    # This will take AT MOST n steps (i.e. one hypothesis rejected at each iteration).
    for (j in 1:n) {
      # Setup
      R_rej <- which(R)
      R_unrej <- which(!R)
      w <- rep(0, n) # For weights on each node.
      unrej_leaves <- intersect(leaves, R_unrej) # Leaves that have not been rejected yet.
      v <- leaves_weights_full[unrej_leaves] # Weights for the unrejected leaves.
      w[unrej_leaves] <- v # Initial weight assignment for an iteration.

      # Calculate topological ordering of the DAG
      topo_unrej <- topo[topo %in% R_unrej]

      if (version == "all") {

        # Determine which hypotheses can be considered for rejection next
        no_unrej_parents <- which(sapply(pa, FUN = function(x) all(x %in% R_rej)))
        candidates <- intersect(no_unrej_parents, R_unrej)

        # Update weights
        for (H in topo_unrej) {
          pa_H <- pa[[H]]
          unrej_pa <- R_unrej[R_unrej %in% pa_H]

          if (length(unrej_pa) > 0) {
            for (i in unrej_pa) {
              w[i] <- w[i] + w[H] / length(unrej_pa)
            }
            w[H] <- 0
          }
        }
      } else {

        # Determine which hypotheses can be considered for rejection next
        one_rej_parent <- which(sapply(pa, FUN = function(x) any(x %in% R_rej)))
        one_rej_parent <- union(one_rej_parent, roots)
        R_unrej <- which(!R)
        candidates <- intersect(one_rej_parent, R_unrej)

        for (H in topo_unrej) {
          pa_H <- pa[[H]]
          unrej_pa <- R_unrej[R_unrej %in% pa_H]
          if (length(unrej_pa) > 0) {
            for (i in unrej_pa) {
              w[i] <- w[i] + w[H] / length(pa_H)
            }
            w[H] <- w[H] * (1 - length(unrej_pa) / length(pa_H))
          }
        }
      }


      # Test hypotheses by comparing p-values with corrected thresholds
      alpha_H <- alpha * w / sum(v)
      small_p <- (p <= alpha_H)

      # Stop procedure if no more rejections occur
      if (all(!small_p)) break

      # If rejections occured, update the rejection list R
      R[candidates] <- small_p[candidates]

      # free rejections to make R congruent
      if (version == "any") {
        reject_for_free <- unlist(an[R])
        R[reject_for_free] <- TRUE
      }

      # can stop if everything has been rejected
      if (all(R)) break
    }

    R
  }
}
