
#' Community Detection function
#'
#' Community detection method: Partition K phenotypes into disjoint network
#' modules using the hierarchical clustering method based on the constructed
#' signed bipartite network, GPN - Genotype and Phenotype Network.
#' The number of network modules and the final cluster results are determined
#' by the perturbation procedure.
#'
#' @param Net M by K matrix of the signed bipartite genotype and
#'  phenotype network, row denotes the M SNPs and column denotes K phenotypes.
#' @param npert The number of perturbation times. default: 100
#' @param mMax The maximum searching number of clusters (mMax < K). default: K
#' @param cluster.method The agglomeration method to be used in the
#'  hierarchical clustering method. This should be one of "ward.D",
#'  "ward.D2","single","complete","average","mcquitty","median","centroid".
#'  default: "ward.D2"
#' @param quiet decide if show the perturbation times. default:
#'  FALSE (show procedure)
#'
#' @return A list of ClusterNumber: the optimal number of clusters (C);
#'  Clusters: K by C matrix of the optimal cluster results;
#'  AF: the area under the curve of CDF calculated from the perturbation results.
#'
#' @importFrom zoo fortify.zoo
#' @importFrom methods hasArg
#' @import dplyr
#' @export
#' @examples
#' M <- 100
#' K <- 10
#' Sigma <- diag(1, nrow = K)
#' Net <- MASS::mvrnorm(M, rep(0,K), Sigma)
#' res <- CommunityDetection(Net, npert = 100, quiet = TRUE)
#'

CommunityDetection <- function(Net, npert = 100, mMax = "None", cluster.method = "ward.D2",
                               quiet = FALSE) {
  requireNamespace("dplyr", quietly = TRUE)
  # Check the Bipartite Network
  if (!hasArg(Net)) {
    stop("Error: Please input the Network matrix! You can obtain from ConstructBiNet function.")
  } else {
    if (sum(!is.finite(Net)) != 0) {
      temp <- fortify.zoo(Net) %>%
        mutate_all(function(x) ifelse(!is.finite(x), 0, x))
      Net <- as.matrix(temp)[, -1]
    } else {
      Net <- as.matrix(Net)
    }
  }

  # Set the maximum searching number of clusters
  if (mMax == "None") {
    mMax <- ncol(Net)
  } else {
    if (mMax > ncol(Net)) {
      stop("Error: The maximum searching number of clusters must smaller
           than the number of phenotypes!")
    }
  }

  # Check agglomeration method in the hierarchical clustering
  # method
  cMethod <- c("ward.D", "ward.D2", "single", "complete", "average",
               "mcquitty", "median", "centroid")
  if (!cluster.method %in% cMethod) {
    stop(paste("Error: the cluster.method must be one of", paste(cMethod,
                                                                 collapse = ", "), "!"))
  }

  # main function
  n <- nrow(Net)
  m <- ncol(Net)
  S <- Sp <- D <- B <- list()
  for (i in 2:(mMax - 1)) S[[i - 1]] <- Sp[[i - 1]] <- matrix(0, m, m)
  for (i in 2:(mMax - 1)) B[[i - 1]] <- matrix(0, m, i)
  sigma2 <- sqrt(median(apply(Net, 1, var))) * 1

  AF <- rep(0, mMax - 2)
  Sigma <- cor(Net)
  dist <- 1 - Sigma
  hc <- hclust(as.dist(dist), method = cluster.method)

  for (i in 2:(mMax - 1)) {
    index <- cutree(hc, i)
    for (j in 1:i) {
      a <- ((index >= j) & (index < j + 1))
      S[[i - 1]] <- S[[i - 1]] + a %*% t(a)
      B[[i - 1]][, j] <- a
    }
  }
  for (k in 1:npert) {
    ypert <- Net + matrix(sigma2 * rnorm(n * m), n, m)
    Sigma <- cor(ypert)
    dist <- 1 - Sigma
    hc = hclust(as.dist(dist), method = cluster.method)
    for (i in 2:(mMax - 1)) {
      index <- cutree(hc, i)
      for (j in 1:i) {
        a <- ((index >= j) & (index < j + 1))
        Sp[[i - 1]] <- Sp[[i - 1]] + a %*% t(a)
      }
    }
    if (!quiet) {
      print(paste("Perturbation times:", k, "among", npert, "!"))
    }
  }
  for (i in 2:(m - 1)) {
    Sp[[i - 1]] <- Sp[[i - 1]]/npert
    D[[i - 1]] <- abs(S[[i - 1]] - Sp[[i - 1]])
    for (k in 1:100) AF[i - 1] <- AF[i - 1] + sum((D[[i - 1]] <= (k -
                                                                    0.5)/100))
  }
  L.max <- length(AF)
  IND <- which.max(abs(AF[-1] - AF[-L.max])) + 1
  L <- length(IND)
  if (L > 1)
    IND <- IND[L]
  khat <- IND
  Bb <- B[[khat - 1]]
  return(list(ClusterNumber = khat, Clusters = Bb, AF = AF))
}
