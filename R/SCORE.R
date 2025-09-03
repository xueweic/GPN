

#' SCORE function
#'
#' Score test for testing the association between phenotypes and SNPs.
#' Note that the phenotypes could be the binary phenotypes or
#' quantitative phenotypes.
#'
#' \code{SCORE} Default method;
#'
#' \code{SPA_SCORE} Score test by saddlepoint approximation for testing the
#'  association between only binary phenotypes and SNPs. Note that SPA_SCORE
#'  function can use to binary phenotypes with the unbalanced or extremely
#'  unbalanced cast-control ratios.
#'
#' @aliases SCORE SPA_SCORE
#' @param x The individual-level genotyes with dimension n individuals by M SNPs.
#'  It could be only one SNPs with a n by 1 genotype vector.
#' @param y The individual-level phentypes with dimesion n individuals by K
#'  phenotypes. It could be only one phenotype with a n by 1 phenotype vector.
#' @param cov The covariate matrix with dimesion n individuals by C covariates.
#'  defalut: there is no covariates.
#' @param output.T Determine if output the test statistics. default:
#'  output.T = FALSE
#' @import stats
#' @return A list of pvalue: K by M matrix of p-values; Tstat: K by M matrix
#'  of the test statistics.
#' @export
#' @examples
#' N <- 100
#' M <- 200
#' K <- 10
#' x <- replicate(M, rbinom(N,2,0.3))
#' y1 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.9, 0.1)))
#' y2 <- replicate(2, rnorm(N))
#' res <- SCORE(x, y2)
#' res1 <- SPA_SCORE(x, y1)
#'


SCORE <- function(x, y, cov = NULL) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  # Check dimension of the individual-level data: x and y
  if (nrow(x) != nrow(y)) {
    stop("Error: Please check the sample size of x and y.
         They should be the same!")
  }

  # Check the existing of the covariance matrix
  if (!is.null(cov)) {
    cov <- as.matrix(cov)
    if (nrow(cov) != nrow(x)) {
      stop("Error: Please check the sample size of covariance matrix,
           which should be the same as the sample size of  individual-level
           genotyps and phenotypes!")
    }
    x1 <- apply(x, 2, function(t) return(residuals(lm(t ~ 1 + cov))))
    y1 <- apply(y, 2, function(t) {
      if (all(t %in% c(0, 1))) {
        return(residuals(glm(t ~ 1 + cov, family = binomial(link = "logit"))))
      } else {
        return(residuals(lm(t ~ 1 + cov)))
      }
    })
    x <- as.matrix(x1)
    y <- as.matrix(y1)
  }

  n <- nrow(x)
  Tstat <- sqrt(n) * cor(y, x)
  pv <- pchisq(Tstat^2, 1, lower.tail = FALSE)
  return(list(pvalue = pv, Tstat = Tstat))
}




#' @rdname SCORE
#' @importFrom SPAtest ScoreTest_SPA
#' @export
SPA_SCORE <- function(x, y, cov = NULL, output.T = FALSE) {
  requireNamespace("SPAtest", quietly = TRUE)

  x <- as.matrix(x)
  y <- as.matrix(y)
  # Check dimension of the individual-level data: x and y
  if (nrow(x) != nrow(y)) {
    stop("Error: Please check the sample size of x and y. They should be the same!")
  }
  # Check if binary phenotype
  if (!all(y %in% c(0, 1))) {
    stop("Error: SPA_SCORE only can be applied to the binary phenotypes")
  }
  # Check the existing of the covariance matrix
  if (!is.null(cov)) {
    cov <- as.matrix(cov)
    if (nrow(cov) != nrow(x)) {
      stop("Error: Please check the sample size of covariance matrix,
           which should be the same as the sample size of  individual-level
           genotyps and phenotypes!")
    }
  }

  M <- ncol(x)
  K <- ncol(y)
  pv <- matrix(1, ncol = M, nrow = K)
  if (!output.T) {
    for (i.phen in 1:K) {
      res <- ScoreTest_SPA(genos = t(x), pheno = as.numeric(y[, i.phen]),
                           cov = cov, method = "fastSPA")
      pv[i.phen, ] <- res$p.value
    }
    return(list(pvalue = pv))
  } else if (output.T) {
    Tstat <- Betahat <- seBetahat <- matrix(1, ncol = ncol(x), nrow = ncol(y))
    for (i.phen in 1:K) {
      res <- ScoreTest_SPA(genos = t(x), pheno = as.numeric(y[, i.phen]),
                           cov = cov, method = "fastSPA", beta.out = TRUE, beta.Cutoff = 1)
      pv[i.phen, ] <- res$p.value
      Tstat[i.phen, ] <- res$beta/res$SEbeta
      Betahat[i.phen, ] <- res$beta
      seBetahat[i.phen, ] <- res$SEbeta
    }
    return(list(pvalue = pv, Tstat = Tstat, Betahat = Betahat, seBetahat = seBetahat))
  }
}





