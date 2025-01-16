


#' Improved Community Detection method
#'
#' @param Net M by K matrix of the signed bipartite genotype and
#'  phenotype network, row denotes the M SNPs and column denotes K phenotypes.
#' @param between_cluster Threshold for hierachical clusters. Default value is 0.8.
#'
#' @return A list of ClusterNumber: the optimal number of clusters (L);
#'  Cluster: K by L matrix of the optimal cluster results;
#'  Q_modularity: modularity information.
#' @export
#'
#' @examples
#' M <- 100
#' K <- 10
#' Sigma <- diag(1, nrow = K)
#' Net <- MASS::mvrnorm(M, rep(0,K), Sigma)
#' res <- CommunityDetection2(Net, between_cluster = 0.8)

CommunityDetection2 <- function(Net, between_cluster = 0.8){
  # - perform hierachical clustering appraoch
  cormat <- cor(Net)
  hc <- hclust(as.dist(1-cormat))
  # - get the optimal number of cluster
  opt_cluster <- get_n_cluster(hc, cormat, between_cluster = between_cluster)
  n_cluster <- opt_cluster$n_cluster
  Q_modularity <- opt_cluster$Qmodularity
  # - obtain the final clusters
  index <- cutree(hc,n_cluster)
  B <- sapply(1:n_cluster, function(t) as.numeric(index==t))
  B <- as.matrix(B)
  L <- ncol(B)
  return(list("ClusterNumber"= L,
              "cluster" = B,
              "Q_modularity" = Q_modularity))
}





get_n_cluster <- function(hc, Sigma, m=ncol(Sigma), between_cluster = 0.8){
  if (min(Sigma) > between_cluster){
    IND = 1
    Q = 1
  } else {
    Q <- c()
    if (ncol(Sigma) < 10){m = ncol(Sigma)}
    for(i in 1:m)
    {
      index=cutree(hc,i)
      B=sapply(1:i, function(t) as.numeric(index==t))
      Q[i] <- get_modularity(Sigma, B)
    }

    IND=which(Q==max(Q))
    L=length(IND)
    if (L>1) IND=IND[1]
  }
  return(list("n_cluster" = IND,
              "Qmodularity" = Q))
}



get_modularity <- function(Weight, B){
  # ------ Calculate modularity ----------
  if (dim(Weight)[1] == 1){
    Q <- 0
  } else {
    W_pos <- Weight * (Weight > 0)
    W_neg <- Weight * (Weight < 0)
    N <- dim(Weight)[1]
    K_pos <- colSums(W_pos)
    K_neg <- colSums(W_neg)
    m_pos <- sum(K_pos)
    m_neg <- sum(K_neg)
    m <- m_pos + m_neg
    cate <- B %*% t(B)
    if (m_pos == 0 & m_neg == 0){
      Q <- 0
    } else {
      if (m_pos == 0){
        Q_positive <- 0
        Q_negative <- sum((W_neg - K_neg %*% t(K_neg) / m_neg) * cate) / m_neg
      } else if (m_neg == 0){
        Q_positive <- sum((W_pos - K_pos %*% t(K_pos) / m_pos) * cate) / m_pos
        Q_negative <- 0
      } else {
        Q_positive <- sum((W_pos - K_pos %*% t(K_pos) / m_pos) * cate) / m_pos
        Q_negative <- sum((W_neg - K_neg %*% t(K_neg) / m_neg) * cate) / m_neg
      }
    }
    Q <- m_pos / m * Q_positive - m_neg / m * Q_negative
  }
}

