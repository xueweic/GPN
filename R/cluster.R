
#' The cluster coefficient for bipartite network
#'
#' The cluster coefficient measures the extent to which
#' two SNPs share the same set of associated phenotypes.
#' A high cluster coefficient indicates that SNPs tend
#' to be associated with similar phenotypes, while
#' a low cluster coefficient indicates that SNPs
#' tend to be associated with distinct sets of phenotypes.
#' The cluster coefficient can be used to identify SNP
#' sub-networks that are highly interconnected and
#' potentially functionally related. (It will be time consuming
#' for high dimensional matrix)
#'
#' @param Net A adjacency matrix with genotype or genes as rows,
#'  phenotype as columns and the association as entries.
#'  The adjacency matrix of the bipartite network which
#'  could be a simple matrix or data frame.
#'
#' @return A list of object type containing: first value
#'  which is the cluster coefficient of all SNPs/genes,
#'  the second value is the cluster coefficient of all phenotypes,
#'  the third value is the cluster coefficient of whole network,
#'  the last value is the average cluster coefficient of SNPs and
#'  phenotypes.
#' @export
#' @keywords GPN
#' @examples
#' Net <- matrix(c(1,0,0,1,
#'               1,1,0,0,
#'               0,1,1,0),nrow = 4, ncol = 3)
#' out <- cluster(Net)


cluster <- function(Net) {
  if (!is.matrix(Net))
    Net <- as.matrix(Net)
  M <- nrow(Net)
  K <- ncol(Net)
  Net <- check_sparsity(Net)[[1]]

  sub_fun <- function(Net) {
    M <- nrow(Net)
    N1 <- sapply(1:M, function(i) {
      which(Net[i, ] != 0)
    }, simplify = FALSE, USE.NAMES = TRUE)

    c <- c()
    for (i.node in 1:M) {
      # print(paste0('i.node=', i.node))
      pos <- N1[[i.node]]
      temp <- as.matrix(Net[, pos])
      N2 <- which(rowSums(temp) != 0)
      N2 <- N2[-which(N2 == i.node)]
      cc <- sum(unlist(sapply(N2, function(i) {
        length(intersect(N1[[i.node]], N1[[i]]))/length(union(N1[[i.node]],
                                                              N1[[i]]))
      })))/length(N2)
      c <- c(c, cc)
    }
    c[is.na(c)] <- 0
    names(c) <- rownames(Net)
    return(c)
  }

  cc_SNP <- sub_fun(Net)
  cc_phenotype <- sub_fun(t(Net))
  list(cluster_SNP = sum(cc_SNP)/M,
       cluster_phenotype = sum(cc_phenotype)/K,
       cluster_all = sum(c(cc_SNP, cc_phenotype))/(K + M),
       cluster_SUM = sum(cc_SNP)/M + sum(cc_phenotype)/K)
}
