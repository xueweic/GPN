
#' Modified Computational Efficient Clustering linear combination (ceCLC) test
#'
#' ceCLC test for testing the association between K phenotypes with a SNP. The phenotypes can be either qualitative or binary, espectially the binary phenotypes with the extremely case-control ratio (the test statistics has been adjusted by the saddlepoint approximation)
#'
#' @param x The genotyes for only one SNP with dimension n by 1.
#' @param y The phentypes with dimension n individuals by K phenotypes.
#' @param cov The covariate matrix with dimension n individuals by C covariates. defalut: there is no covariates.
#' @param cluster.method The agglomeration method to be used in the hierarchical clustering method. This should be one of "ward.D","ward.D2","single","complete","average","mcquitty","median","centroid". default: "ward.D2"
#'
#' @return p-value of the ceCLC test statistic
#' @export
#'
#' @examples
#' N <- 10000
#' y1 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.999, 0.001)))
#' y2 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.5, 0.5)))
#' y3 <- replicate(2, rnorm(N))
#' y <- cbind(y1, y2, y3)
#' x <- rbinom(N,2,0.3)
#' ceCLC(x, y, cov = NULL)
#'
ceCLC <- function(x, y, cov = NULL, cluster.method = "ward.D2")
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
  }
  # Check agglomeration method in the hierarchical clustering method
  cMethod <- c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")
  if (!cluster.method %in% cMethod){
    stop(paste("Error: the cluster.method must be one of", paste(cMethod, collapse = ", "), "!"))
  }

  # score test statistic
  Tstat <- apply(y, 2, function(ytemp){
    if (!all(ytemp %in% c(0,1))){
      stat <- SCORE(x, ytemp, cov)$Tstat
    } else {
      ccratio <- sum(ytemp)/length(ytemp)
      stat <- ifelse((ccratio >= 0.4 & ccratio <= 0.6),
                     SCORE(x, ytemp, cov)$Tstat,
                     Adjusted_Tstat(x,ytemp,cov = NULL))
    }
    return(stat)
  })
  Tstat <- as.matrix(Tstat)


  # Check if adjust for the covariates
  if (!is.null(cov)){
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


  L0 <- ncol(y)
  if (L0 == 1){
    CLC <- Tstat^2
    pv0 <- pchisq(CLC,L0,lower.tail = FALSE)
  } else {
    Sigma <- cor(y)
    dist <- 1-Sigma
    hc <- hclust(as.dist(dist), method = cluster.method)
    pv0 <- rep(9999,L0)
    U <- list()
    for (L in c(1:L0))
    {
      index1 <- cutree(hc,L)
      B <- sapply(1:L, function(t) as.numeric(index1==t))
      W <- t(B)%*%ginv_ratio(Sigma)
      U[[L]] <- t(W)%*%ginv_ratio(W%*%Sigma%*%t(W))%*%W
      CLC <- t(Tstat)%*%U[[L]]%*%Tstat
      pv0[L] <- pchisq(CLC,L,lower.tail = FALSE)
    }
  }
  is.small <- (pv0<1e-15)
  pv0[!is.small] <- tan((0.5-pv0[!is.small])*pi)
  pv0[is.small] <- 1/pv0[is.small]/pi
  ACAT <- mean(pv0)
  pvalue <- 0.5-atan(ACAT)/pi
  return(pvalue)
}
