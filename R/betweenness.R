
#' The Approximate Betweenness Centrality of Bipartite Network
#'
#' In the bipartite genetic variant and phenotype network (GPN),
#' the Approximate Betweenness Centrality is used to evaluate
#' the importance of genetic variant. For instance, the genetic variant
#' have a high Approximate Betweenness Centrality which means
#' an important connector between phenotypes.
#'
#' @param Net A adjacency matrix with genotype or genes as rows,
#'  phenotype as columns and the association as entries.
#'  The adjacency matrix of the bipartite network which
#'  could be a simple matrix or data frame.
#'
#' @return A list of object type containing: the first vector is
#'  is the Approximate Betweenness Centrality of genetic variants,
#'  the second vector is standardized Approximate Betweenness
#'  Centralityof genetic variants.
#' @references Latapy, M., Magnien, C., & Del Vecchio, N. (2008).
#'  \emph{Basic notions for the analysis of large two-mode networks.}
#'  Social networks, 30(1), 31-48.
#' @export
#' @keywords GPN
#' @examples
#' Net <- matrix(c(1,1,0,
#'               1,1,1,
#'               1,1,1),nrow = 3, ncol = 3, byrow = TRUE)
#' out <- betweenness(Net)


betweenness <- function(Net) {
  if (!is.matrix(Net))
    Net <- as.matrix(Net)
  Net[Net != 0] <- 1
  K <- ncol(Net)
  M <- nrow(Net)
  tl <- sapply(1:K, function(i) {
    list(which(Net[, i] != 0))
  })
  names(tl) <- colnames(Net)
  B <- rep(0, M)

  for (i in 1:(K - 1)) {
    for (j in (i + 1):K) {
      ind.snp <- intersect(tl[[i]], tl[[j]])
      if (length(ind.snp) > 0) {
        B[ind.snp] <- B[ind.snp] + 1/length(ind.snp)
      }
    }
  }
  names(B) <- rownames(Net)
  list(Betweenness_SNP = B, Betweenness_SNP_sd = B/sum(B))
}
