

#' Modified clustering linear combination method based on hierarchical clustering (HCLC)
#'
#' HCLC test for testing the association between K phenotypes with a SNP. The phenotypes can be either qualitative or binary, espectially the binary phenotypes with the extremely case-control ratio (the test statistics has been adjusted by the saddlepoint approximation)
#'
#' Reference:
#' Li, X., Zhang, S., & Sha, Q. (2020). Joint analysis of multiple phenotypes using a clustering linear combination method based on hierarchical clustering. Genetic epidemiology, 44(1), 67-78.
#' Liang, X., Cao, X., Sha, Q., & Zhang, S. (2022). HCLC-FC: a novel statistical method for phenome-wide association studies. bioRxiv.
#'
#' @param x The genotyes for only one SNP with dimension n by 1.
#' @param y The phentypes with dimension n individuals by K phenotypes.
#' @param cov The covariate matrix with dimension n individuals by C covariates. defalut: there is no covariates.
#' @param cluster.method The agglomeration method to be used in the hierarchical clustering method. This should be one of "ward.D","ward.D2","single","complete","average","mcquitty","median","centroid". default: "average"
#'
#' @return p-value of the HCLC test statistic
#' @export
#'
#' @examples
#' N <- 10000
#' y1 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.999, 0.001)))
#' y2 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.5, 0.5)))
#' y3 <- replicate(2, rnorm(N))
#' y <- cbind(y1, y2, y3)
#' x <- rbinom(N,2,0.3)
#' HCLC(x, y)
#'
HCLC <- function(x, y, cov = NULL, cluster.method = "average"){
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

  # Test
  K <- ncol(y)
  if (K == 1){
    pv <- pchisq(Tstat^2, 1, lower.tail = FALSE)
  } else if ( K == 2){
    Sigma <- cor(y)
    W <- ginv_ratio(Sigma)
    U <- t(W)%*%ginv_ratio(W%*%Sigma%*%t(W))%*%W
    CLC <- t(Tstat)%*%U%*%Tstat
    pv <- pchisq(CLC,2,lower.tail = FALSE)
  } else {
    # Hierarchical clustering based on the number of clusters with max height difference
    Sigma <- as.matrix(cor(y))
    dist <- 1 - Sigma
    hc <- hclust(as.dist(dist),cluster.method)
    H <- which.max(diff(hc$height))+1
    index <- cutree(hc,H)
    L <- max(index)
    B <- sapply(1:L, function(t) as.numeric(index==t))
    # CLC test
    W <- ginv_ratio( t(B)%*%ginv_ratio(Sigma)%*%B ) %*% t(B) %*% ginv_ratio(Sigma)
    CLC <- t( W%*%Tstat ) %*% ginv_ratio( W%*%Sigma%*%t(W) ) %*% ( W%*%Tstat )
    pv <- pchisq(CLC,L,lower.tail = FALSE)
  }
  return(pv)
}
