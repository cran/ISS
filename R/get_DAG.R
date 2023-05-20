#' get_DAG
#'
#' This function is used to construct the induced DAG, induced polyforest and
#' reverse topological orderings thereof from a numeric matrix \code{X0}. See
#' Definition 2 in \insertCite{MRCS2023;textual}{ISS}.
#'
#' @param X0 a numeric matrix.
#' @param sparse logical. Either the induced DAG (\code{FALSE}) or the induced
#' polyforest (\code{TRUE}) is constructed.
#' @param twoway logical. If \code{FALSE}, only leaves, parents, ancestors
#' and reverse topological ordering are returned. If TRUE, then roots, children
#' and descendants are also provided.
#'
#'
#' @return A list with named elements giving the leaves, parents, ancestors and
#' reverse topological ordering and additionally, if \code{twoway == TRUE}, the
#' roots, children and descendants, of the constructed graph.
#' @export
#'
#' @references \insertRef{MRCS2023}{ISS}
#'
#' @examples
#' X <- rbind(
#'   c(0.2, 0.8), c(0.2, 0.8), c(0.1, 0.7),
#'   c(0.2, 0.1), c(0.3, 0.5), c(0.3, 0)
#' )
#' get_DAG(X0 = X)
#' get_DAG(X0 = X, sparse = TRUE, twoway = TRUE)
#'

get_DAG <- function(X0, sparse = FALSE, twoway = FALSE) {
  if (is.vector(X0) & is.atomic(X0)) X0 <- matrix(X0, ncol = 1)
  if (!(is.numeric(X0) & is.matrix(X0))) stop("X0 needs to be a numeric matrix.")

  d <- ncol(X0)
  n <- nrow(X0)

  upper_points <- vector(n, mode = "list")
  for (i in 1:n) {
    x_i <- X0[i, ]
    # duplicate values with larger index (subsetting shifts index)
    duplicates_i_incl_i <- which(colSums(t(X0[i:n, ]) == x_i) == d) + (i - 1)
    # strictly larger values
    strictly_larger_i <- which(colSums(t(X0) >= x_i) == d & colSums(t(X0) == x_i) < d)
    upper_points[[i]] <- c(setdiff(duplicates_i_incl_i, i), strictly_larger_i)
  }

  if (sparse == FALSE) {
    ancestors <- upper_points

    if (n >= 2000) {
      no_cores <- parallel::detectCores()
      cl <- parallel::makeCluster(no_cores)
      parents <- parallel::parLapply(
        cl = cl, 1:length(ancestors),
        fun = function(i, ancestors) {
          an_i <- ancestors[[i]]
          an_i[!(an_i %in% unlist(ancestors[an_i]))]
        },
        ancestors
      )
      parallel::stopCluster(cl)
    } else {
      parents <- ancestors
      for (i in 1:n) {
        an_i <- ancestors[[i]]
        parents[[i]] <- an_i[!(an_i %in% unlist(ancestors[an_i]))]
      }
    }
  } else {
    parents <- vector(n, mode = "list") # technically an atomic vector suffices
    children <- vector(n, mode = "list")

    for (i in 1:n) {
      x_i <- X0[i, ]
      upper_i <- upper_points[[i]]
      X_upper <- X0[upper_i, , drop = FALSE]
      parent_ind <- which.min(apply(t(X_upper) - x_i, MARGIN = 2, max))
      pa_i <- upper_i[parent_ind]
      parents[[i]] <- pa_i
      if (length(pa_i) > 0) children[[pa_i]] <- c(children[[pa_i]], i)
    }

    leaves <- which(sapply(children, FUN = length) == 0)
    ancestors <- as.list(data.frame(matrix(integer(0), ncol = n)))
    names(ancestors) <- NULL

    for (j in leaves) {
      chain <- rep(NA, n)
      i <- j
      for (k in 1:n) {
        chain[k] <- i
        pa_i <- parents[[i]]
        if (length(pa_i) == 0) break
        i <- pa_i
      }

      chain <- chain[!is.na(chain)]
      K <- length(chain)
      if (K > 1) {
        for (k in 1:(K - 1)) {
          i <- chain[k]
          ancestors[[i]] <- chain[(k + 1):K]
        }
      }
    }
  }


  leaves <- (1:n)[-unlist(ancestors)]

  topo_list <- vector(n, mode = "list")
  topo_list[[1]] <- leaves
  remaining <- setdiff(1:n, leaves) # those points are not in the ordering yet
  for (k in 2:n) {
    # those points come after points that are not in the ordering yet
    later_in_topo_ordering <- unlist(ancestors[remaining])

    if (length(later_in_topo_ordering) == 0) {
      topo_list[[k]] <- remaining
      break
    }

    # those can already be included in the ordering
    # (contains those already included)
    can_include <- (1:n)[-later_in_topo_ordering]

    # actual new vertices
    include_now <- intersect(can_include, remaining)

    # update list and vector "remaining"
    topo_list[[k]] <- include_now
    remaining <- setdiff(remaining, include_now)
  }
  topological_ordering <- unlist(topo_list)

  res <- list(
    ancestors = ancestors, parents = parents,
    leaves = leaves, topological_ordering = topological_ordering
  )


  if (twoway & !sparse) {
    roots <- which(sapply(parents, length) == 0)

    children <- vector(n, mode = "list")
    descendants <- vector(n, mode = "list")

    for (i in topological_ordering) {
      children[[i]] <- which(sapply(parents, FUN = function(x) i %in% x))
      descendants[[i]] <- integer(0)
      if (length(children[[i]]) > 0) {
        for (k in children[[i]]) {
          descendants[[i]] <- c(descendants[[i]], k, descendants[[k]])
        }
      }
      descendants[[i]] <- sort(descendants[[i]])
    }

    res$roots <- roots
    res$children <- children
    res$descendants <- descendants
  } else if (twoway & sparse) { # have already determined children in this case
    roots <- which(sapply(parents, length) == 0)

    descendants <- vector(n, mode = "list")

    for (i in topological_ordering) {
      descendants[[i]] <- integer(0)
      if (length(children[[i]]) > 0) {
        for (k in children[[i]]) {
          descendants[[i]] <- c(descendants[[i]], k, descendants[[k]])
        }
      }
      descendants[[i]] <- sort(descendants[[i]])
    }

    res$roots <- roots
    res$children <- children
    res$descendants <- descendants
  }

  return(res)
}
