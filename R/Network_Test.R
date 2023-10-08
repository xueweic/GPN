#'
#'
#' Multiple phenotype association test based on the network modules
#' derived from genotype and phenotype network (GPN)
#'
#' Multiple phenotype association test for testing the association between
#' K phenotypes with a SNP. The phenotypes can be either qualitative or binary,
#' espectially the binary phenotypes with the extremely case-control ratio
#' (the test statistics has been adjusted by the saddlepoint approximation)
#'
#'
#' @param x The genotyes for only one SNP with dimension n by 1.
#' @param y The phentypes with dimension n individuals by K phenotypes.
#' @param B The clustering results of K phenotypes with dimension K by C.cate
#'   categories.
#' @param test.method The multiple phenotype association test. This should be
#'   one of "ceCLC","CLC","HCLC","OBrien","omnibus","TATES","acat.test","HMPhmp.test",
#'   "MultiPhen". default: "acat.test"
#' @param cov The covariate matrix with dimension n individuals by C covariates.
#'   defalut: there is no covariates.
#' @param cluster.method The agglomeration method to be used in the hierarchical
#'   clustering method. This should be one of "ward.D","ward.D2","single",
#'   "complete","average","mcquitty","median","centroid". default: "ward.D2"
#'
#' @return A list of pv.Modules and pv.NET. pv.Modules is a C.cate-dimension
#'   vector of p-values, where pvalues are the association results for testing
#'   the SNP and phenotypes in each of network modules. pv.NET is the pvalue
#'   of modified ceCLC with network modules.
#' @export
#'
#' @examples
#' N <- 10000
#' y1 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.999, 0.001)))
#' y2 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.5, 0.5)))
#' y3 <- replicate(2, rnorm(N))
#' y <- cbind(y1, y2, y3)
#' x <- rbinom(N,2,0.3)
#' B <- matrix(c(1,1,0,0,0,0,0,0,1,1,0,0,0,0,0,0,1,1), nrow = 6)
#' out <- Network_Test(x, y, B, test.method = "acat.test")
#'
#'
Network_Test <- function(x, y, B, test.method = "acat.test", cov = NULL, cluster.method = "ward.D2"){
  x <- as.matrix(x)
  y <- as.matrix(y)
  B <- as.matrix(B)
  # Check dimension of the individual-level data: x and y
  if (nrow(x) != nrow(y)){
    stop("Error: Please check the sample size of x and y. They should be the same!")
  }
  # Check the existing of the covariance matrix
  if (!is.null(cov)){
    cov <- as.matrix(cov)
    if (nrow(cov) != nrow(x)){
      stop("Error: Please check the sample size of covariance matrix,
           which should be the same as the sample size of  individual-level
           genotyps and phenotypes!")
    }
  }
  # Check the cluster
  if (nrow(B) != ncol(y)){
    stop("Error: Check the cluster results!")
  }

  # Check agglomeration method in the hierarchical clustering method
  cMethod <- c("ward.D","ward.D2","single","complete","average","mcquitty","median","centroid")
  if (!cluster.method %in% cMethod){
    stop(paste("Error: the cluster.method must be one of", paste(cMethod, collapse = ", "), "!"))
  }

  # test method
  tMethod <- c("ceCLC","CLC","HCLC","OBrien","omnibus","TATES","acat.test","hmp.test","MultiPhen")
  if (!test.method %in% tMethod){
    stop(paste("Error: the test.method must be one of", paste(tMethod, collapse = ", "), "!"))
  }

  n.C <- ncol(B)
  num.y <- apply(B, 2, sum)
  pvalue <- c()
  for (i.c in 1:n.C)
  {
    pos <- which(B[, i.c] == 1)
    y.cate <- as.matrix(y[,pos])
    if (test.method == "ceCLC"){
      pvalue[i.c] <- ceCLC(x,y.cate,cov,cluster.method)
    } else if (test.method == "CLC"){
      pvalue[i.c] <- CLC(x,y.cate,cov,cluster.method)
    } else if (test.method == "HCLC"){
      pvalue[i.c] <- HCLC(x,y.cate,cov,cluster.method)
    } else if (test.method == "OBrien"){
      pvalue[i.c] <- OBrien(x,y.cate,cov)
    } else if (test.method == "omnibus"){
      pvalue[i.c] <- omnibus(x,y.cate,cov)
    } else if (test.method == "TATES"){
      pvalue[i.c] <- TATES(x,y.cate,cov)
    } else if (test.method == "acat.test") {
      pvalue[i.c] <- acat.test(x,y.cate,cov)
    } else if (test.method == "hmp.test") {
      pvalue[i.c] <- hmp.test(x,y.cate,cov)
    } else if (test.method == "MultiPhen"){
      pvalue[i.c] <- OBrien(x,y.cate,cov)
    }
  }
  padj <- length(pvalue) * min(pvalue)

  return(list("pv.Modules" = pvalue,
              "pv.NET" = padj))
}
