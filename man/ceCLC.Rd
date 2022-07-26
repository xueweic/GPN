% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ceCLC.R
\name{ceCLC}
\alias{ceCLC}
\title{Modified Computational Efficient Clustering linear combination (ceCLC) test}
\usage{
ceCLC(x, y, cov = NULL, cluster.method = "ward.D2")
}
\arguments{
\item{x}{The genotyes for only one SNP with dimension n by 1.}

\item{y}{The phentypes with dimension n individuals by K phenotypes.}

\item{cov}{The covariate matrix with dimension n individuals by C covariates. defalut: there is no covariates.}

\item{cluster.method}{The agglomeration method to be used in the hierarchical clustering method. This should be one of "ward.D","ward.D2","single","complete","average","mcquitty","median","centroid". default: "ward.D2"}
}
\value{
p-value of the ceCLC test statistic
}
\description{
ceCLC test for testing the association between K phenotypes with a SNP. The phenotypes can be either qualitative or binary, espectially the binary phenotypes with the extremely case-control ratio (the test statistics has been adjusted by the saddlepoint approximation)
}
\details{
Reference:
Wang, M., Zhang, S., & Sha, Q. (2022). A computationally efficient clustering linear combination approach to jointly analyze multiple phenotypes for GWAS. PloS one, 17(4), e0260911.
}
\examples{
N <- 10000
y1 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.999, 0.001)))
y2 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.5, 0.5)))
y3 <- replicate(2, rnorm(N))
y <- cbind(y1, y2, y3)
x <- rbinom(N,2,0.3)
ceCLC(x, y)

}
