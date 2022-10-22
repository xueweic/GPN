

#' Adjusted test statistics calculation
#'
#' Calculate the adjusted test statistics for the binary phenotypes with the extremely unbalanced case-control ratios.
#'
#' @param x The individual-level genotyes with dimension n individuals by M SNPs. It could be only one SNPs with a n by 1 genotype vector.
#' @param y The individual-level binary phentypes with the extremely unbalanced case-control ratios with dimesion n individuals by K phenotypes. It could be only one phenotype with a n by 1 phenotype vector.
#' @param cov The covariate matrix with dimesion n individuals by C covariates. defalut: there is no covariates.
#'
#' @return M by K matrix of adjusted test statistics
#' @export
#'
#' @examples
#' N <- 10000
#' M <- 10
#' y <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.999, 0.001))) # two binary phenotypes with the extremely unbalanced case-contorl ratios.
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
    stop("Error: Only the test statistics for testing the assocition between BINARY phenotypes with extremely unbalanced case-control ratios and SNPs need to be adjusted by the saddlepoint approximation!")
  }
  # Check the existing of the covariance matrix
  if (!is.null(cov)){
    cov <- as.matrix(cov)
    if (nrow(cov) != nrow(x)){
      stop("Error: Please check the sample size of covariance matrix, which should be the same as the sample size of  individual-level genotyps and phenotypes!")
    }
  }

  BiNet <- GPN(x, y, cov, output.T = FALSE)
  AdjTstat <- sign(BiNet) * sqrt(abs(BiNet))
  return(AdjTstat)
}
