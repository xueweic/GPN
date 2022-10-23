

#' Modified O'Brien test
#'
#' O'Brien test for testing the association between K phenotypes with a SNP. The phenotypes can be either qualitative or binary, espectially the binary phenotypes with the extremely case-control ratio (the test statistics has been adjusted by the saddlepoint approximation)
#'
#' Reference:
#' O'Brien, P. C. (1984). Procedures for comparing samples with multiple endpoints. Biometrics, 1079-1087.
#'
#' @param x The genotyes for only one SNP with dimension n by 1.
#' @param y The phentypes with dimension n individuals by K phenotypes.
#' @param cov The covariate matrix with dimension n individuals by C covariates. defalut: there is no covariates.
#'
#' @return p-value of the OBrien test statistic
#' @export
#'
#' @examples
#' N <- 10000
#' y1 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.999, 0.001)))
#' y2 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.5, 0.5)))
#' y3 <- replicate(2, rnorm(N))
#' y <- cbind(y1, y2, y3)
#' x <- rbinom(N,2,0.3)
#' OBrien(x, y)
#'
#'
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
