

#' Construct a signed bipartite genotype and phenotype network
#'
#' @param x The individual-level genotyes with dimension n individuals by M SNPs.
#' @param y The individual-level phentypes with dimesion n individuals by K phenotypes.
#' @param cov The covariate matrix with dimesion n individuals by C covariates. defalut: there is no covariates.
#' @param output.T Determine if output the test statistics. default: output.T = FALSE
#'
#' @return M by K matrix of the signed bipartite genotype and phenotype network (GPN)
#' @export
#'
#' @examples
#' N <- 10000
#' M <- 10
#' y1 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.999, 0.001))) # two binary phenotypes with the extremely unbalanced case-contorl ratios.
#' y2 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.5, 0.5))) # two binary phenotypes with the balanced case-contorl ratios.
#' y3 <- replicate(2, rnorm(N)) # two quantitative (continuous) phenotypes
#' y <- cbind(y1, y2, y3)
#' x <- replicate(M, rbinom(N,2,0.3))
#' res <- GPN(x,y)
#'
#'
GPN <- function(x, y, cov = NULL, output.T = FALSE){
  x <- as.matrix(x)
  y <- as.matrix(y)
  # Check dimension of the individual-level data: x and y
  if (nrow(x) != nrow(y)){
    stop("Error: Please check the sample size of x and y. They should be the same!")
  }
  # Check the existing of the covariance matrix
  if (!is.null(cov)){
    cov <- as.matrix(cov)
    if (nrow(cov) != nrow(x)){
      stop("Error: Please check the sample size of covariance matrix, which should be the same as the sample size of  individual-level genotyps and phenotypes!")
    }
  }

  M <- ncol(x)
  K <- ncol(y)
  pv <- matrix(1, ncol = M, nrow = K)
  for (i.phen in 1:K){
    y.temp <- y[,i.phen]

    # For quantitative phenotype
    if (!all(y.temp %in% c(0,1))){
      pv[i.phen, ] <- SCORE(x, y.temp, cov)$pvalue
    } else {

      # For binary phenotype
      ccratio <- sum(y.temp)/length(y.temp)
      if (ccratio >= 0.4 & ccratio <= 0.6){
        pv[i.phen, ] <- SCORE(x, y.temp, cov)$pvalue
      } else {
        pv[i.phen, ] <- SPA_SCORE(x, y.temp, cov, output.T = output.T)$pvalue
      }

    }
  }

  y_mean <- colMeans(y)
  S <- (t(y) - y_mean) %*% x
  BiNet <- t(sign(S) * qchisq(pv, 1, lower.tail=FALSE))
  return(BiNet)
}





