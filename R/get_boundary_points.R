
#' get_boundary_points
#'
#' Given a set of points, returns the minimal subset with the same upper hull.
#'
#' @param X a numeric matrix with one point per row.
#'
#' @return A numeric matrix of the same number of columns as \code{X}.
#' @export
#'
#' @examples
#' X <- rbind(c(0, 1), c(1, 0), c(1, 0), c(1, 1))
#' get_boundary_points(X)
#'
get_boundary_points <- function(X) {
  # given a collection of points, get a minimal number of points with the same
  # upper hull as the original collection
  n <- nrow(X)

  if (n == 0) {
    X # returns the input matrix, so also simply empty output
  } else {
    boundary <- rep(NA, n)
    for (i in 1:n) {
      boundary[i] <- !any(apply(X[-i, , drop = FALSE],
        MARGIN = 1,
        FUN = function(x) all(x <= X[i, ]) & any(x != X[i, ])
      ))
    }
    Xb <- X[boundary, , drop = FALSE]

    unique(Xb)
  }
}
