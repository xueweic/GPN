

#' Calculate R inverse
#'
#' Calculate the inverse of the squared matrix by either ginv function in MASS package or the proposed method, (Yan et al. 2022, Cao et al. 2022).
#'
#' @param R The squared matrix
#' @param ratio The threshold that the smaller eigenvalues to be removed in the calculation of inverse matrix (0 < ratio <=1). default: ratio = 1, that is, the exact inverse matrix; otherwise, ratio < 1 is the adjusted inverse matrix.
#'
#' @return The inverse of matrix R
#' @export
#'
#' @examples
#'
#' y <- replicate(5, rnorm(100))
#' R <- cor(y)
#' ginv_ratio(R)
#'
ginv_ratio <- function(R, ratio = 1){
  R <- as.matrix(R)
  # Check the correlation matrix
  if (ncol(R) != nrow(R)){
    stop("Error: matrix must be squared!")
  }
  if (ncol(R) == 1){
    Rinv <- 1/R
  } else {
    # Check if MASS::ginv can be used
    temp <- try(Rinv <- ginv(R), silent = TRUE)
    if ("try-error" %in% class(temp)){
      Rstar <- R %*% t(R)
      eig <- eigen(Rstar)
      e.value <- eig$values
      e.vec <- eig$vectors
      L1 <- sum(cumsum(e.value)/sum(e.value) < ratio) + 1
      m <- e.vec[,1:L1]
      v <- as.matrix(diag(1/sqrt(e.value[1:L1])))
      Rinv <- m%*%v%*%t(m)
    }
  }
  return(Rinv)
}
