


#' Multiple phenotypes assocaition test
#'
#' Hypothesis test for testing the association between K phenotypes with a SNP.
#' The phenotypes can be either qualitative or binary, espectially the
#' binary phenotypes with the extremely case-control ratio (the test statistics
#' has been adjusted by the saddlepoint approximation)
#'
#' \code{ceCLC()} Modified Computational Efficient Clustering linear combination
#'  (ceCLC) test
#'
#' \code{CLC()} Modified Clustering linear combination (CLC) test.
#'
#' \code{HCLC()} Modified clustering linear combination method based on hierarchical
#'  clustering (HCLC)
#'
#' \code{OBrien()} Modified O'Brien test
#'
#' \code{omnibus()} Modified omnibus (Chi-square) test
#'
#' \code{TATES()} Modified TATES test
#'
#' \code{acat.test()} Cauchy combination test
#'
#' \code{hmp.test()} Harmonic mean p-value
#'
#' \code{MultiPhen()} MultiPhen test
#'
#' @aliases ceCLC CLC HCLC OBrien omnibus TATES acat.test hmp.test MultiPhen
#' @param x The genotyes for only one SNP with dimension n by 1.
#' @param y The phentypes with dimension n individuals by K phenotypes.
#' @param cov The covariate matrix with dimension n individuals by C covariates.
#'  defalut: there is no covariates.
#' @param cluster.method The agglomeration method to be used in the hierarchical
#'  clustering method. This should be one of "ward.D","ward.D2","single",
#'  "complete","average","mcquitty","median","centroid". default: "ward.D2"
#'
#' @return p-value of the test statistic
#'
#' @references Wang, M., Zhang, S., & Sha, Q. (2022). A computationally efficient
#'  clustering linear combination approach to jointly analyze multiple phenotypes
#'  for GWAS. PloS one, 17(4), e0260911.
#'
#' @references Sha, Q., Wang, Z., Zhang, X., & Zhang, S. (2019). A clustering
#'  linear combination approach to jointly analyze multiple phenotypes for GWAS.
#'  Bioinformatics, 35(8), 1373-1379.
#'
#' @references Li, X., Zhang, S., & Sha, Q. (2020). Joint analysis of multiple
#'  phenotypes using a clustering linear combination method based on hierarchical
#'  clustering. Genetic epidemiology, 44(1), 67-78.
#' @references Liang, X., Cao, X., Sha, Q., & Zhang, S. (2022). HCLC-FC:
#'  a novel statistical method for phenome-wide association studies. bioRxiv.
#'
#' @references O'Brien, P. C. (1984). Procedures for comparing samples with
#'  multiple endpoints. Biometrics, 1079-1087.
#'
#' @references Van der Sluis, S., Posthuma, D., & Dolan, C. V. (2013).
#'  TATES: efficient multivariate genotype-phenotype analysis for genome-wide
#'  association studies. PLoS genetics, 9(1), e1003235.
#'
#' @references Liu, Y., Chen, S., Li, Z., Morrison, A. C., Boerwinkle, E., & Lin,
#'  X. (2019). ACAT: a fast and powerful p value combination method for
#'  rare-variant analysis in sequencing studies. The American Journal of
#'  Human Genetics, 104(3), 410-421.
#'
#' @references Wilson, D. J. (2019). The harmonic mean p-value for combining
#'  dependent tests. Proceedings of the National Academy of Sciences, 116(4),
#'  1195-1200.
#'
#'
#' @references O’Reilly, P. F., Hoggart, C. J., Pomyen, Y., Calboli, F. C.,
#'  Elliott, P., Jarvelin, M. R., & Coin, L. J. (2012). MultiPhen: joint model
#'  of multiple phenotypes can increase discovery in GWAS. PloS one, 7(5), e34861.
#'
#' @examples
#' N <- 10000
#' y1 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.999, 0.001)))
#' y2 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.5, 0.5)))
#' y3 <- replicate(2, rnorm(N))
#' y <- cbind(y1, y2, y3)
#' x <- rbinom(N,2,0.3)
#' ceCLC(x, y)
#' CLC(x, y)
#' HCLC(x, y)
#' OBrien(x, y)
#' omnibus(x, y)
#' TATES(x, y)
#' acat.test(x, y)
#' hmp.test(x, y)
#' MultiPhen(x, y)




###### ceCLC #######
#' @rdname Test
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



###### CLC #######
#' @rdname Test
CLC <- function(x, y, cov = NULL, cluster.method = "average"){
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
    # Hierarchical clustering based on the number of clusters with minimum height
    Sigma <- as.matrix(cor(y))
    dist <- 1 - Sigma
    hc <- hclust(as.dist(dist),cluster.method)
    H <- which.min(diff(hc$height))
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



###### HCLC #######
#' @rdname Test
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


###### OBrien #######
#' @rdname Test
OBrien <- function(x, y, cov=NULL)
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
  n <- nrow(y)
  if (K == 1){
    pv <- pchisq(Tstat^2, 1, lower.tail = FALSE)
  } else {
    Sigma <- as.matrix(cor(y))
    U <- rep(1,K) %*% ginv_ratio(Sigma) %*% Tstat
    U1 <- (U/sqrt(sum(ginv_ratio(Sigma))))^2
    pv <- pchisq(U1,1,lower.tail = FALSE)
  }
  return(pv)
}


###### omnibus #######
#' @rdname Test
omnibus <- function(x, y, cov=NULL)
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

  n <- nrow(y)
  K <- ncol(y)
  pv <- pchisq(t(Tstat)%*%ginv_ratio(cor(y))%*%Tstat, K, lower.tail = FALSE)
  return(pv)
}


###### TATES #######
#' @rdname Test
TATES <- function(x, y, cov=NULL)
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


  K <- ncol(y)
  n <- nrow(y)
  if (K == 1)
  {
    pv <- pchisq(Tstat^2, 1, lower.tail = FALSE)
  }
  if (K > 1)
  {
    pv0 <- pchisq(Tstat^2, 1, lower.tail = FALSE)
    pv_j <- sort(pv0, index.return = TRUE)
    idx <- pv_j$ix
    y <- y[,idx]
    r <- cor(y)
    rho <- -0.0008-0.0023*r+0.6226*r^2+0.0149*r^3+0.1095*r^4-0.0219*r^5+0.2179*r^6
    me <- numeric(K)
    for(i in 1:K){
      rhoi <- rho[1:i,1:i]
      eigeni <- eigen(rhoi)$values
      me[i] <- i-sum(ifelse(eigeni>1,eigeni,1)-1)
    }
    pv <- min(me[K]*pv_j$x/me)
  }
  return(pv)
}



###### ACAT #######
#' @rdname Test
#' @importFrom ACAT ACAT
acat.test <- function(x, y, cov=NULL)
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
  pv <- as.vector(pchisq(Tstat^2, 1, lower.tail = FALSE))
  pv <- ACAT(pv)
  return(pv)
}



###### HMP #######
#' @rdname Test
#' @importFrom harmonicmeanp p.hmp
hmp.test <- function(x, y, cov=NULL)
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
  pv <- as.vector(pchisq(Tstat^2, 1, lower.tail = FALSE))
  L <- length(pv)
  pv <- p.hmp(pv, L = L)
  return(pv)
}







###### MultiPhen #######
#' @rdname Test
MultiPhen <- function(x, y, cov=NULL)
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

  # Check if adjust for the covariates
  if (!is.null(cov)){
    y1 <- apply(y, 2, function(t){
      if (all(t %in% c(0,1))){
        return(residuals(glm(t ~ 1+cov, family = binomial(link = "logit"))))
      } else {
        return(residuals(lm(t ~ 1+cov)))
      }
    })
    y <- as.matrix(y1)
  }

  x <- factor(x)
  fit1 <- polr_fixed(x~1)
  fit2 <- polr_fixed(x~y)
  res <- anova(fit1,fit2)
  pv <- res$"Pr(Chi)"[2]
  return(pv)
}



# ---- Fixed the polr function for MultiPhen.

# file MASS/R/polr.R
# copyright (C) 1994-2008 W. N. Venables and B. D. Ripley
# Use of transformed intercepts contributed by David Firth
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 or 3 of the License
#  (at your option).
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/
#
#
polr_fixed <- function(formula, data, weights, start, ..., subset,
                       na.action, contrasts = NULL, Hess = FALSE,
                       model = TRUE,
                       method = c("logistic", "probit", "cloglog", "cauchit"))
{
  logit <- function(p) log(p/(1 - p))

  fmin <- function(beta) {
    theta <- beta[pc + 1L:q]
    gamm <- c(-Inf, cumsum(c(theta[1L], exp(theta[-1L]))), Inf)
    eta <- offset
    if (pc > 0)
      eta <- eta + drop(x %*% beta[1L:pc])
    pr <- pfun(gamm[y + 1] - eta) - pfun(gamm[y] - eta)
    if (all(pr > 0))
      -sum(wt * log(pr))
    else Inf
  }

  gmin <- function(beta)
  {
    jacobian <- function(theta) { ## dgamma by dtheta matrix
      k <- length(theta)
      etheta <- exp(theta)
      mat <- matrix(0 , k, k)
      mat[, 1] <- rep(1, k)
      for (i in 2:k) mat[i:k, i] <- etheta[i]
      mat
    }
    theta <- beta[pc + 1L:q]
    gamm <- c(-Inf, cumsum(c(theta[1L], exp(theta[-1L]))), Inf)
    eta <- offset
    if(pc > 0) eta <- eta + drop(x %*% beta[1L:pc])
    pr <- pfun(gamm[y+1] - eta) - pfun(gamm[y] - eta)
    p1 <- dfun(gamm[y+1] - eta)
    p2 <- dfun(gamm[y] - eta)
    g1 <- if(pc > 0) t(x) %*% (wt*(p1 - p2)/pr) else numeric(0)
    xx <- .polrY1*p1 - .polrY2*p2
    g2 <- - t(xx) %*% (wt/pr)
    g2 <- t(g2) %*% jacobian(theta)
    if(all(pr > 0)) c(g1, g2) else rep(NA, pc+q)
  }

  m <- match.call(expand.dots = FALSE)
  method <- match.arg(method)
  pfun <- switch(method, logistic = plogis, probit = pnorm,
                 cloglog = pgumbel, cauchit = pcauchy)
  dfun <- switch(method, logistic = dlogis, probit = dnorm,
                 cloglog = dgumbel, cauchit = dcauchy)
  if(is.matrix(eval.parent(m$data)))
    m$data <- as.data.frame(data)
  m$start <- m$Hess <- m$method <- m$model <- m$... <- NULL
  m[[1L]] <- as.name("model.frame")
  m <- eval.parent(m)
  Terms <- attr(m, "terms")
  x <- model.matrix(Terms, m, contrasts)
  xint <- match("(Intercept)", colnames(x), nomatch=0L)
  n <- nrow(x)
  pc <- ncol(x)
  cons <- attr(x, "contrasts") # will get dropped by subsetting
  if(xint > 0) {
    x <- x[, -xint, drop=FALSE]
    pc <- pc - 1
  } else warning("an intercept is needed and assumed")
  wt <- model.weights(m)
  if(!length(wt)) wt <- rep(1, n)
  offset <- model.offset(m)
  if(length(offset) <= 1) offset <- rep(0, n)
  y <- model.response(m)
  if(!is.factor(y)) stop("response must be a factor")
  lev <- levels(y)
  if(length(lev) <= 2) stop("response must have 3 or more levels")
  y <- unclass(y)
  q <- length(lev) - 1
  Y <- matrix(0, n, q)
  .polrY1 <- col(Y) == y
  .polrY2 <- col(Y) == y - 1
  if(missing(start)) {
    # try something that should always work -tjb
    u <- as.integer(table(y))
    u <- (cumsum(u)/sum(u))[1:q]
    zetas <-
      switch(method,
             "logistic"= qlogis(u),
             "probit"=   qnorm(u),
             "cauchit"=  qcauchy(u),
             "cloglog"=  -log(-log(u)) )
    s0 <- c(rep(0,pc),zetas[1],log(diff(zetas)))
  } else if(length(start) != pc + q)
    stop("'start' is not of the correct length")
  else {
    s0 <- if(pc > 0) c(start[seq_len(pc+1)], log(diff(start[-seq_len(pc)])))
    else c(start[1L], log(diff(start)))
  }
  res <- optim(s0, fmin, gmin, method="BFGS", hessian = Hess, ...)
  beta <- res$par[seq_len(pc)]
  theta <- res$par[pc + 1L:q]
  zeta <- cumsum(c(theta[1L],exp(theta[-1L])))
  deviance <- 2 * res$value
  niter <- c(f.evals=res$counts[1L], g.evals=res$counts[2L])
  names(zeta) <- paste(lev[-length(lev)], lev[-1L], sep="|")
  if(pc > 0) {
    names(beta) <- colnames(x)
    eta <- offset + drop(x %*% beta)
  } else eta <- offset + rep(0, n)

  cumpr <- matrix(pfun(matrix(zeta, n, q, byrow=TRUE) - eta), , q)
  fitted <- t(apply(cumpr, 1L, function(x) diff(c(0, x, 1))))
  dimnames(fitted) <- list(row.names(m), lev)
  fit <- list(coefficients = beta, zeta = zeta, deviance = deviance,
              fitted.values = fitted, lev = lev, terms = Terms,
              df.residual = sum(wt) - pc - q, edf = pc + q, n = sum(wt),
              nobs = sum(wt),
              call = match.call(), method = method,
              convergence = res$convergence, niter = niter, lp = eta)
  if(Hess) {
    dn <- c(names(beta), names(zeta))
    H <- res$hessian
    dimnames(H) <- list(dn, dn)
    fit$Hessian <- H
  }
  if(model) fit$model <- m
  fit$na.action <- attr(m, "na.action")
  fit$contrasts <- cons
  fit$xlevels <- .getXlevels(Terms, m)
  class(fit) <- "polr"
  fit
}

print.polr <- function(x, ...)
{
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl, control=NULL)
  }
  if(length(coef(x))) {
    cat("\nCoefficients:\n")
    print(coef(x), ...)
  } else {
    cat("\nNo coefficients\n")
  }
  cat("\nIntercepts:\n")
  print(x$zeta, ...)
  cat("\nResidual Deviance:", format(x$deviance, nsmall=2), "\n")
  cat("AIC:", format(x$deviance + 2*x$edf, nsmall=2), "\n")
  if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
  if(x$convergence > 0)
    cat("Warning: did not converge as iteration limit reached\n")
  invisible(x)
}

vcov.polr <- function(object, ...)
{
  jacobian <- function(theta) { ## dgamma by dtheta matrix
    k <- length(theta)
    etheta <- exp(theta)
    mat <- matrix(0 , k, k)
    mat[, 1] <- rep(1, k)
    for (i in 2:k) mat[i:k, i] <- etheta[i]
    mat
  }

  if(is.null(object$Hessian)) {
    message("\nRe-fitting to get Hessian\n")
    utils::flush.console()
    object <- update(object, Hess=TRUE,
                     start=c(object$coefficients, object$zeta))
  }
  vc <- ginv(object$Hessian)
  pc <- length(coef(object))
  gamma <- object$zeta
  z.ind <- pc + seq_along(gamma)
  theta <- c(gamma[1L], log(diff(gamma)))
  J <- jacobian(theta)
  A <- diag(pc + length(gamma))
  A[z.ind, z.ind] <- J
  V <- A %*% vc %*% t(A)
  structure(V,  dimnames = dimnames(object$Hessian))
}

summary.polr <- function(object, digits = max(3, .Options$digits - 3),
                         correlation = FALSE, ...)
{
  cc <- c(coef(object), object$zeta)
  pc <- length(coef(object))
  q <- length(object$zeta)
  coef <- matrix(0, pc+q, 3, dimnames=list(names(cc),
                                           c("Value", "Std. Error", "t value")))
  coef[, 1] <- cc
  vc <- vcov(object)
  coef[, 2] <- sd <- sqrt(diag(vc))
  coef[, 3] <- coef[, 1]/coef[, 2]
  object$coefficients <- coef
  object$pc <- pc
  object$digits <- digits
  if(correlation)
    object$correlation <- (vc/sd)/rep(sd, rep(pc+q, pc+q))
  class(object) <- "summary.polr"
  object
}

print.summary.polr <- function(x, digits = x$digits, ...)
{
  if(!is.null(cl <- x$call)) {
    cat("Call:\n")
    dput(cl, control=NULL)
  }
  coef <- format(round(x$coefficients, digits=digits))
  pc <- x$pc
  if(pc > 0) {
    cat("\nCoefficients:\n")
    print(x$coefficients[seq_len(pc), , drop=FALSE], quote = FALSE,
          digits = digits, ...)
  } else {
    cat("\nNo coefficients\n")
  }
  cat("\nIntercepts:\n")
  print(coef[(pc+1):nrow(coef), , drop=FALSE], quote = FALSE,
        digits = digits, ...)
  cat("\nResidual Deviance:", format(x$deviance, nsmall=2), "\n")
  cat("AIC:", format(x$deviance + 2*x$edf, nsmall=2), "\n")
  if(nzchar(mess <- naprint(x$na.action))) cat("(", mess, ")\n", sep="")
  if(!is.null(correl <- x$correlation)) {
    cat("\nCorrelation of Coefficients:\n")
    ll <- lower.tri(correl)
    correl[ll] <- format(round(correl[ll], digits))
    correl[!ll] <- ""
    print(correl[-1, -ncol(correl)], quote = FALSE, ...)
  }
  invisible(x)
}

predict.polr <- function(object, newdata, type=c("class","probs"), ...)
{
  if(!inherits(object, "polr")) stop("not a \"polr\" object")
  type <- match.arg(type)
  if(missing(newdata)) Y <- object$fitted
  else {
    newdata <- as.data.frame(newdata)
    Terms <- delete.response(object$terms)
    m <- model.frame(Terms, newdata, na.action = function(x) x,
                     xlev = object$xlevels)
    if (!is.null(cl <- attr(Terms, "dataClasses")))
      .checkMFClasses(cl, m)
    X <- model.matrix(Terms, m, contrasts = object$contrasts)
    xint <- match("(Intercept)", colnames(X), nomatch=0L)
    if(xint > 0) X <- X[, -xint, drop=FALSE]
    n <- nrow(X)
    q <- length(object$zeta)
    eta <- drop(X %*% object$coefficients)
    pfun <- switch(object$method, logistic = plogis, probit = pnorm,
                   cloglog = pgumbel, cauchit = pcauchy)
    cumpr <- matrix(pfun(matrix(object$zeta, n, q, byrow=TRUE) - eta), , q)
    Y <- t(apply(cumpr, 1L, function(x) diff(c(0, x, 1))))
    dimnames(Y) <- list(rownames(X), object$lev)
  }
  if(missing(newdata) && !is.null(object$na.action))
    Y <- napredict(object$na.action, Y)
  switch(type, class = {
    Y <- factor(max.col(Y), levels=seq_along(object$lev),
                labels=object$lev)
  }, probs = {})
  drop(Y)
}

extractAIC.polr <- function(fit, scale = 0, k = 2, ...)
{
  edf <- fit$edf
  c(edf, deviance(fit) + k * edf)
}

model.frame.polr <- function(formula, ...)
{
  dots <- list(...)
  nargs <- dots[match(c("data", "na.action", "subset"), names(dots), 0)]
  if(length(nargs) || is.null(formula$model)) {
    m <- formula$call
    m$start <- m$Hess <- m$... <- NULL
    m[[1L]] <- as.name("model.frame")
    m[names(nargs)] <- nargs
    if (is.null(env <- environment(formula$terms))) env <- parent.frame()
    data <- eval(m, env)
    if(!is.null(mw <- m$weights)) {
      nm <- names(data)
      nm[match("(weights)", nm)] <- as.character(mw)
      names(data) <- nm
    }
    data
  } else formula$model
}

pgumbel <- function(q, loc = 0, scale = 1, lower.tail = TRUE)
{
  q <- (q - loc)/scale
  p <- exp(-exp(-q))
  if (!lower.tail) 1 - p else p
}

dgumbel <- function (x, loc = 0, scale = 1, log = FALSE)
{
  x <- (x - loc)/scale
  d <- log(1/scale) - x - exp(-x)
  d[is.nan(d)] <- -Inf                # -tjb
  if (!log) exp(d) else d
}

anova.polr <- function (object, ..., test = c("Chisq", "none"))
{
  test <- match.arg(test)
  dots <- list(...)
  if (length(dots) == 0L)
    stop('anova is not implemented for a single "polr" object')
  mlist <- list(object, ...)
  nt <- length(mlist)
  dflis <- sapply(mlist, function(x) x$df.residual)
  s <- order(dflis, decreasing = TRUE)
  mlist <- mlist[s]
  if (any(!sapply(mlist, inherits, "polr")))
    stop('not all objects are of class "polr"')
  ns <- sapply(mlist, function(x) length(x$fitted.values))
  if(any(ns != ns[1L]))
    stop("models were not all fitted to the same size of dataset")
  rsp <- unique(sapply(mlist, function(x) paste(formula(x)[2L])))
  mds <- sapply(mlist, function(x) paste(formula(x)[3L]))
  dfs <- dflis[s]
  lls <- sapply(mlist, function(x) deviance(x))
  tss <- c("", paste(1L:(nt - 1), 2:nt, sep = " vs "))
  df <- c(NA, -diff(dfs))
  x2 <- c(NA, -diff(lls))
  pr <- c(NA, 1 - pchisq(x2[-1L], df[-1L]))
  out <- data.frame(Model = mds, Resid.df = dfs, Deviance = lls,
                    Test = tss, Df = df, LRtest = x2, Prob = pr)
  names(out) <- c("Model", "Resid. df", "Resid. Dev", "Test",
                  "   Df", "LR stat.", "Pr(Chi)")
  if (test == "none") out <- out[, 1L:6]
  class(out) <- c("Anova", "data.frame")
  attr(out, "heading") <-
    c("Likelihood ratio tests of ordinal regression models\n",
      paste("Response:", rsp))
  out
}

polr.fit <- function(x, y, wt, start, offset, method)
{
  logit <- function(p) log(p/(1 - p))

  fmin <- function(beta) {
    theta <- beta[pc + 1L:q]
    gamm <- c(-Inf, cumsum(c(theta[1L], exp(theta[-1L]))), Inf)
    eta <- offset
    if (pc > 0)
      eta <- eta + drop(x %*% beta[1L:pc])
    pr <- pfun(gamm[y + 1] - eta) - pfun(gamm[y] - eta)
    if (all(pr > 0))
      -sum(wt * log(pr))
    else Inf
  }

  gmin <- function(beta)
  {
    jacobian <- function(theta) { ## dgamma by dtheta matrix
      k <- length(theta)
      etheta <- exp(theta)
      mat <- matrix(0 , k, k)
      mat[, 1] <- rep(1, k)
      for (i in 2:k) mat[i:k, i] <- etheta[i]
      mat
    }
    theta <- beta[pc + 1L:q]
    gamm <- c(-Inf, cumsum(c(theta[1L], exp(theta[-1L]))), Inf)
    eta <- offset
    if(pc > 0) eta <- eta + drop(x %*% beta[1L:pc])
    pr <- pfun(gamm[y+1] - eta) - pfun(gamm[y] - eta)
    p1 <- dfun(gamm[y+1] - eta)
    p2 <- dfun(gamm[y] - eta)
    g1 <- if(pc > 0) t(x) %*% (wt*(p1 - p2)/pr) else numeric(0)
    xx <- .polrY1*p1 - .polrY2*p2
    g2 <- - t(xx) %*% (wt/pr)
    g2 <- t(g2) %*% jacobian(theta)
    if(all(pr > 0)) c(g1, g2) else rep(NA, pc+q)
  }

  pfun <- switch(method, logistic = plogis, probit = pnorm,
                 cloglog = pgumbel, cauchit = pcauchy)
  dfun <- switch(method, logistic = dlogis, probit = dnorm,
                 cloglog = dgumbel, cauchit = dcauchy)
  n <- nrow(x)
  pc <- ncol(x)
  lev <- levels(y)
  if(length(lev) <= 2L) stop("response must have 3 or more levels")
  y <- unclass(y)
  q <- length(lev) - 1L
  Y <- matrix(0, n, q)
  .polrY1 <- col(Y) == y
  .polrY2 <- col(Y) == y - 1L
  # pc could be 0.
  s0 <- if(pc > 0) c(start[seq_len(pc+1)], diff(start[-seq_len(pc)]))
  else c(start[1L], diff(start))
  res <- optim(s0, fmin, gmin, method="BFGS")
  beta <- res$par[seq_len(pc)]
  theta <- res$par[pc + 1L:q]
  zeta <- cumsum(c(theta[1L],exp(theta[-1L])))
  deviance <- 2 * res$value
  names(zeta) <- paste(lev[-length(lev)], lev[-1L], sep="|")
  if(pc > 0) {
    names(beta) <- colnames(x)
    eta <- drop(x %*% beta)
  } else {
    eta <- rep(0, n)
  }
  list(coefficients = beta, zeta = zeta, deviance = deviance)
}

profile.polr <- function(fitted, which = 1L:p, alpha = 0.01,
                         maxsteps = 10, del = zmax/5, trace = FALSE, ...)
{
  Pnames <- names(B0 <- coefficients(fitted))
  pv0 <- t(as.matrix(B0))
  p <- length(Pnames)
  if(is.character(which)) which <- match(which, Pnames)
  summ <- summary(fitted)
  std.err <- summ$coefficients[, "Std. Error"]
  mf <- model.frame(fitted)
  n <- length(Y <- model.response(mf))
  O <- model.offset(mf)
  if(!length(O)) O <- rep(0, n)
  W <- model.weights(mf)
  if(length(W) == 0L) W <- rep(1, n)
  OriginalDeviance <- deviance(fitted)
  X <- model.matrix(fitted)[, -1L, drop=FALSE] # drop intercept
  zmax <- sqrt(qchisq(1 - alpha, 1))
  profName <- "z"
  prof <- vector("list", length=length(which))
  names(prof) <- Pnames[which]
  start <- c(fitted$coefficients, fitted$zeta)
  for(i in which) {
    zi <- 0
    pvi <- pv0
    Xi <- X[,  - i, drop = FALSE]
    pi <- Pnames[i]
    for(sgn in c(-1, 1)) {
      if(trace) {
        message("\nParameter:", pi, c("down", "up")[(sgn + 1)/2 + 1])
        utils::flush.console()
      }
      step <- 0
      z <- 0
      ## LP is the linear predictor including offset.
      ## LP <- X %*% fitted$coef + O
      while((step <- step + 1) < maxsteps && abs(z) < zmax) {
        bi <- B0[i] + sgn * step * del * std.err[i]
        o <- O + X[, i] * bi
        fm <- polr.fit(x = Xi, y = Y, wt = W, start = start[-i],
                       offset = o, method = fitted$method)
        ri <- pv0
        ri[, names(coef(fm))] <- coef(fm)
        ri[, pi] <- bi
        pvi <- rbind(pvi, ri)
        zz <- fm$deviance - OriginalDeviance
        if(zz > - 1e-3) zz <- max(zz, 0)
        else stop("profiling has found a better solution, so original fit had not converged")
        z <- sgn * sqrt(zz)
        zi <- c(zi, z)
      }
    }
    si <- order(zi)
    prof[[pi]] <- structure(data.frame(zi[si]), names = profName)
    prof[[pi]]$par.vals <- pvi[si, ]
  }
  val <- structure(prof, original.fit = fitted, summary = summ)
  class(val) <- c("profile.polr", "profile")
  val
}

confint.polr <- function(object, parm, level = 0.95, trace = FALSE, ...)
{
  pnames <- names(coef(object))
  if(missing(parm)) parm <- seq_along(pnames)
  else if(is.character(parm))  parm <- match(parm, pnames, nomatch = 0L)
  message("Waiting for profiling to be done...")
  utils::flush.console()
  object <- profile(object, which = parm, alpha = (1. - level)/4.,
                    trace = trace)
  confint(object, parm=parm, level=level, trace=trace, ...)
}

confint.profile.polr <-
  function(object, parm = seq_along(pnames), level = 0.95, ...)
  {
    of <- attr(object, "original.fit")
    pnames <- names(coef(of))
    if(is.character(parm))  parm <- match(parm, pnames, nomatch = 0L)
    a <- (1-level)/2
    a <- c(a, 1-a)
    pct <- paste(round(100*a, 1), "%")
    ci <- array(NA, dim = c(length(parm), 2L),
                dimnames = list(pnames[parm], pct))
    cutoff <- qnorm(a)
    for(pm in parm) {
      pro <- object[[ pnames[pm] ]]
      if(length(pnames) > 1L)
        sp <- spline(x = pro[, "par.vals"][, pm], y = pro[, 1])
      else sp <- spline(x = pro[, "par.vals"], y = pro[, 1])
      ci[pnames[pm], ] <- approx(sp$y, sp$x, xout = cutoff)$y
    }
    drop(ci)
  }

logLik.polr <- function(object, ...)
  structure(-0.5 * object$deviance, df = object$edf, class = "logLik")

simulate.polr <- function(object, nsim = 1, seed = NULL, ...)
{
  if(!is.null(object$model) && any(model.weights(object$model) != 1))
    stop("weighted fits are not supported")

  rgumbel <- function(n, loc = 0, scale = 1) loc - scale*log(rexp(n))

  ## start the same way as simulate.lm
  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)                     # initialize the RNG if necessary
  if(is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  rfun <- switch(object$method, logistic = rlogis, probit = rnorm,
                 cloglog = rgumbel, cauchit = rcauchy)
  eta <- object$lp
  n <- length(eta)
  res <- cut(rfun(n*nsim, eta),
             c(-Inf, object$zeta, Inf),
             labels = colnames(fitted(object)),
             ordered_result = TRUE)
  val <- split(res, rep(seq_len(nsim), each=n))
  names(val) <- paste("sim", seq_len(nsim), sep="_")
  val <- as.data.frame(val)
  if (!is.null(nm <- rownames(fitted(object)))) row.names(val) <- nm
  attr(val, "seed") <- RNGstate
  val
}

