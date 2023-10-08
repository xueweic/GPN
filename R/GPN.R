

#' Construct a signed bipartite genotype and phenotype network based on
#' individual data
#'
#' @param x The individual-level genotyes with dimension n individuals by M SNPs.
#' @param y The individual-level phentypes with dimesion n individuals by K phenotypes.
#' @param \dots Arguments passed to the internal function \code{gpn_train}.
#'  Argument "replicates": the number of replicates is used to choose optimal
#'  parameter \code{tau}, default is 50; Argument "weight": Logical obejct,
#'  if TRUE, bipartite projection will include edge weights. Defaults to TRUE.
#'  Argument "metric": a string specifying which entropy to use.
#'  Possible values are: \code{"kl_divergence_snp"},
#'  \code{"kl_divergence_phe"},\code{"kl_divergence_all"},
#'  \code{"cross_entropy_snp"},\code{"cross_entropy_phe"},
#'  \code{"cross_entropy_all"}. The default is
#'  \code{"kl_divergence_all"}. Argument "tuneGrid": the vector with possible
#'  tuning values \code{tau}. It's should between 0 to 1, default value
#'  from 0.05 to 0.95 with a step size of 0.05.
#' @param cov The covariate matrix with dimesion n individuals by C covariates.
#'  defalut: there is no covariates.
#' @param well.Net Logical object. If TRUE, return a well-defined GPN; If FALSE,
#'  return a denser GPN. Default is TRUE.
#'
#' @return M by K matrix of the signed bipartite genotype and phenotype network (GPN)
#' @export
#'
#' @examples
#' N <- 100
#' M <- 200
#' y1 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.9, 0.1)))
#' y2 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.5, 0.5)))
#' y3 <- replicate(2, rnorm(N))
#' y <- cbind(y1, y2, y3)
#' x <- replicate(M, rbinom(N,2,0.3))
#' res <- GPN(x,y)


GPN <- function(x, y, cov = NULL, well.Net = TRUE, ...) {
  x <- as.matrix(x)
  y <- as.matrix(y)
  # Check dimension of the individual-level data: x and y
  if (nrow(x) != nrow(y)) {
    stop("Error: Please check the sample size of x and y. They should be the same!")
  }
  # Check the existing of the covariance matrix
  if (!is.null(cov)) {
    cov <- as.matrix(cov)
    if (nrow(cov) != nrow(x)) {
      stop("Error: Please check the sample size of covariance matrix,
           which should be the same as the sample size of individual-level
           genotyps and phenotypes!")
    }
  }

  M <- ncol(x)
  K <- ncol(y)
  pv <- matrix(1, ncol = M, nrow = K)
  for (i.phen in 1:K) {
    y.temp <- y[, i.phen]
    # For quantitative phenotype
    if (!all(y.temp %in% c(0, 1))) {
      pv[i.phen, ] <- SCORE(x, y.temp, cov)$pvalue
    } else {
      # For binary phenotype
      ccratio <- sum(y.temp)/length(y.temp)
      if (ccratio >= 0.4 & ccratio <= 0.6) {
        pv[i.phen, ] <- SCORE(x, y.temp, cov)$pvalue
      } else {
        pv[i.phen, ] <- SPA_SCORE(x, y.temp, cov, ...)$pvalue
      }

    }
  }
  y_mean <- colMeans(y)
  S <- (t(y) - y_mean) %*% x
  Net <- t(sign(S) * qchisq(pv, 1, lower.tail = FALSE))
  if (well.Net == TRUE) {
    out <- gpn_train(Net, ...)
  } else {
    out <- Net
  }
  return(out)
}









