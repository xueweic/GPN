

#' Connectance of bipartie network
#'
#' Connectance of bipartite network represent the proportion
#' of realized interactions from the pool of all possible interactions
#' between two kinds of node in a network.
#'
#' @param Net A adjacency matrix with genotype or genes as rows,
#'  phenotype as columns, and the association as entries.
#'  The adjacency matrix of the bipartite network which
#'  could be a simple matrix or data frame.
#'
#' @return A list of object type containing: a value for
#'  bipartite network connectance.
#' @export
#' @keywords GPN
#' @examples
#' \dontrun{
#' set.seed(123)
#' Net <- matrix(rnorm(2000)^2, nrow = 500, ncol = 4)
#' out <- connectance(Net)
#' }


connectance <- function(Net){
  list(connectance=sum(Net!=0)/prod(dim(Net)))
}


