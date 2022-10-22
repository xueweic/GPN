
#' SPA_SCORE function
#'
#' Score test by saddlepoint approximation for testing the association between only binary phenotypes and SNPs. Note that SPA_SCORE function can use to binary phenotypes with the unbalanced or extremely unbalanced cast-control ratios.
#'
#' @param x The individual-level genotyes with dimension n individuals by M SNPs. It could be only one SNPs with a n by 1 genotype vector.
#' @param y The individual-level binary phentypes with dimesion n individuals by K phenotypes. It could be only one phenotype with a n by 1 phenotype vector.
#' @param cov The covariate matrix with dimesion n individuals by C covariates. defalut: there is no covariates.
#' @param output.T Determine if output the test statistics. default: output.T = FALSE
#'
#' @return A list of pvalue: K by M matrix of p-values. If output.T = TRUE, Tstat: K by M matrix of the test statistics.
#' @export
#'
#' @examples
#' N <- 10000
#' M <- 5
#' K <- 10
#' x <- replicate(M, rbinom(N,2,0.3))
#' y <- replicate(K, sample(c(0,1), N, replace = TRUE, prob = c(0.999, 0.001)))
#' res <- SPA_SCORE(x,y)
#'
SPA_SCORE <- function(x, y, cov = NULL, output.T = FALSE)
{
  requireNamespace("SPAtest", quietly = TRUE)

  x <- as.matrix(x)
  y <- as.matrix(y)
  # Check dimension of the individual-level data: x and y
  if (nrow(x) != nrow(y)){
    stop("Error: Please check the sample size of x and y. They should be the same!")
  }
  # Check if binary phenotype
  if (!all(y %in% c(0,1))){
    stop("Error: SPA_SCORE only can be applied to the binary phenotypes")
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
  if (!output.T){
    for (i.phen in 1:K){
      res <- ScoreTest_SPA(genos = t(x), pheno = as.numeric(y[,i.phen]), cov = cov, method = "fastSPA")
      pv[i.phen, ] <- res$p.value
    }
    return(list("pvalue" = pv))
  } else if (output.T){
    Tstat <- matrix(1, ncol = ncol(x), nrow = ncol(y))
    for (i.phen in 1:K){
      res <- ScoreTest_SPA(genos = t(x), pheno = as.numeric(y[,i.phen]), cov = cov, method = "fastSPA",
                           beta.out=TRUE, beta.Cutoff=1)
      pv[i.phen, ] <- res$p.value
      Tstat[i.phen, ] <- res$beta / res$SEbeta
    }
    return(list("pvalue" = pv,
                "Tstat" = Tstat))
  }
}
