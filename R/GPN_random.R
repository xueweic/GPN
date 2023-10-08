

#' Generate a random bipartite GPN
#'
#' Generates a random bipartite GPN, in which
#' the degrees of all vertices are the same as the original
#' bipartite GPN.
#'
#' @param Net A adjacency matrix with genotype or genes as rows,
#'  phenotype as columns and the association as entries.
#'  The adjacency matrix of the bipartite network which
#'  could be a simple matrix or data frame.
#'
#' @return A list of object type containing: the generated
#'  random bipartite GPN.
#' @export
#' @keywords GPN
#' @examples
#' \dontrun{
#' set.seed(123)
#' Net <- matrix(rnorm(2000)^2, nrow = 500, ncol = 4)
#' out <- GPN_random(Net)
#' }
#'

GPN_random <- function(Net) {

  if (!is.matrix(Net))
    Net <- as.matrix(Net)
  # Get the indices and value of the nonzero elements
  nonzero_indices <- which(Net != 0, arr.ind = TRUE)
  edge.weights <- Net[nonzero_indices]

  # choose the links
  edges <- which(abs(Net) >= 0, arr.ind = TRUE)
  n.e <- length(edge.weights)
  ind <- sample(1:nrow(edges), size = n.e, replace = FALSE)
  Net.R <- matrix(0, nrow = nrow(Net), ncol = ncol(Net))
  Net.R[edges[ind, ]] <- edge.weights
  list(Net.R = Net.R)

}
