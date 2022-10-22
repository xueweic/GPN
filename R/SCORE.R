




#' SCORE function
#'
#' Score test for testing the association between phenotypes and SNPs. Note that the phenotypes could be the binary phenotypes with the balanced case-control ratios or quantitative phenotypes. If you consider the phenotypes with the unbalanced or extremely unbalanced cast-control ratios, you need to use SPA_SCORE function.
#'
#' @param x The individual-level genotyes with dimension n individuals by M SNPs. It could be only one SNPs with a n by 1 genotype vector.
#' @param y The individual-level phentypes with dimesion n individuals by K phenotypes. It could be only one phenotype with a n by 1 phenotype vector.
#' @param cov The covariate matrix with dimesion n individuals by C covariates. defalut: there is no covariates.
#'
#' @return A list of pvalue: K by M matrix of p-values; Tstat: K by M matrix of the test statistics.
#' @export
#'
#' @examples
#' N <- 100
#' M <- 5
#' K <- 10
#' x <- replicate(M, rbinom(N,2,0.3))
#' y <- replicate(K, rnorm(N))
#' res <- SCORE(x,y)
#'
SCORE <- function(x,y,cov=NULL)
{
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
    x1 <- apply(x, 2, function(t) return(residuals(lm(t ~ 1+cov))))
    y1 <- apply(y, 2, function(t){
      if (all(t %in% c(0,1))){
        return(residuals(glm(t ~ 1+cov, family = binomial(link = "logit"))))
      } else {
        return(residuals(lm(t ~ 1+cov)))
      }
    })
    x <- as.matrix(x1)
    y <- as.matrix(y1)
  }

  n <- nrow(x)
  Tstat <- sqrt(n)*cor(y,x)
  pv <- pchisq(Tstat^2,1,lower.tail = FALSE)
  return(list("pvalue" = pv,
              "Tstat" = Tstat))
}

