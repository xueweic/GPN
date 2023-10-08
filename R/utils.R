
### internal function

#' check the adjacency matrix since some nodes have no coonection
#'
#' @param Net A adjacency matrix with genotype or genes as rows,
#'  phenotype as columns and the association as entries.
#'  The adjacency matrix of the bipartite network which
#'  could be a simple matrix or data frame.
#'
#' @return
#'
#' @examples
#' Net <- matrix(rnorm(2000)^2, nrow = 500, ncol = 4)
#' out <- check_sparsity(Net)
#'
check_sparsity <- function(Net) {
  if (!is.matrix(Net))
    Net <- as.matrix(Net)
  ind.SNP <- which(rowSums(Net) != 0)
  ind.phe <- which(colSums(Net) != 0)
  if (length(ind.phe) > 0) {
    Net <- Net[, ind.phe]
  }
  if (length(ind.SNP) > 0) {
    Net <- Net[ind.SNP, ]
  }
  list(Net = Net, index_SNP = ind.SNP, index_phe = ind.phe)
}


######
# internal function: extract the well defined network
adjust_adjacency <- function(Net, p, tau, method = "qvalue", ...) {
  Net <- as.matrix(Net)
  temp_p <- adjust_p(p, method, ...)
  ad.p <- as.matrix(temp_p$adjust_p)
  Net[ad.p > tau] <- 0
  list(Sparsity.Net = Net)
}




##########
#' Calculate R inverse
#'
#' Calculate the inverse of the squared matrix by either ginv function
#' in MASS package or the proposed method, (Yan et al. 2022, Cao et al. 2022).
#'
#' @param R The squared matrix
#' @param ratio The threshold that the smaller eigenvalues to be removed
#'   in the calculation of inverse matrix (0 < ratio <=1). default: ratio = 1,
#'   that is, the exact inverse matrix; otherwise, ratio < 1 is
#'   the adjusted inverse matrix.
#' @return The inverse of matrix R
#' @import MASS ginv
#' @examples
#' y <- replicate(5, rnorm(100))
#' R <- cor(y)
#' ginv_ratio(R)

ginv_ratio <- function(R, ratio = 1){
  R <- as.matrix(R)
  # Check the correlation matrix
  if (ncol(R) != nrow(R)){
    stop("Error: matrix must be squared!")
  }
  if (ncol(R) == 1){
    Rinv <- 1/R
  } else {
    # Check if MASS::ginv can be used
    temp <- try(Rinv <- ginv(R), silent = TRUE)
    if ("try-error" %in% class(temp)){
      Rstar <- R %*% t(R)
      eig <- eigen(Rstar)
      e.value <- eig$values
      e.vec <- eig$vectors
      L1 <- sum(cumsum(e.value)/sum(e.value) < ratio) + 1
      m <- e.vec[,1:L1]
      v <- as.matrix(diag(1/sqrt(e.value[1:L1])))
      Rinv <- m%*%v%*%t(m)
    }
  }
  return(Rinv)
}



########
#' Adjusted test statistics calculation
#'
#' Calculate the adjusted test statistics for the binary phenotypes
#' with the extremely unbalanced case-control ratios.
#'
#' @param x The individual-level genotypes with dimension n individuals by M SNPs.
#'   It could be only one SNPs with a n by 1 genotype vector.
#' @param y The individual-level binary phenotypes with the extremely unbalanced
#'   case-control ratios with dimension n individuals by K phenotypes. It could be
#'   only one phenotype with a n by 1 phenotype vector.
#' @param cov The covariate matrix with dimesion n individuals by C covariates.
#'   default: there is no covariates.
#'
#' @return M by K matrix of adjusted test statistics
#' @keywords Test
#'
#' @examples
#' N <- 10000
#' M <- 10
#' y <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.999, 0.001)))
#' x <- replicate(M, rbinom(N,2,0.3))
#' Adjusted_Tstat(x,y)
#'
Adjusted_Tstat <- function(x,y,cov = NULL){

  x <- as.matrix(x)
  y <- as.matrix(y)
  # Check dimension of the individual-level data: x and y
  if (nrow(x) != nrow(y)){
    stop("Error: Please check the sample size of x and y. They should be the same!")
  }
  # Check if binary phenotype
  if (!all(y %in% c(0,1))){
    stop("Error: Only the test statistics for testing the assocition between BINARY
         phenotypes with extremely unbalanced case-control ratios and SNPs need to
         be adjusted by the saddlepoint approximation!")
  }
  # Check the existing of the covariance matrix
  if (!is.null(cov)){
    cov <- as.matrix(cov)
    if (nrow(cov) != nrow(x)){
      stop("Error: Please check the sample size of covariance matrix, which
           should be the same as the sample size of  individual-level genotyps
           and phenotypes!")
    }
  }

  BiNet <- GPN(x, y, cov, output.T = FALSE, well.Net = FALSE)
  AdjTstat <- sign(BiNet) * sqrt(abs(BiNet))
  return(AdjTstat)
}

