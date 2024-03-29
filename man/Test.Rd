% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/test_utils.R
\name{ceCLC}
\alias{ceCLC}
\alias{CLC}
\alias{HCLC}
\alias{OBrien}
\alias{omnibus}
\alias{TATES}
\alias{acat.test}
\alias{hmp.test}
\alias{MultiPhen}
\title{Multiple phenotypes assocaition test}
\usage{
ceCLC(x, y, cov = NULL, cluster.method = "ward.D2")

CLC(x, y, cov = NULL, cluster.method = "average")

HCLC(x, y, cov = NULL, cluster.method = "average")

OBrien(x, y, cov = NULL)

omnibus(x, y, cov = NULL)

TATES(x, y, cov = NULL)

acat.test(x, y, cov = NULL)

hmp.test(x, y, cov = NULL)

MultiPhen(x, y, cov = NULL)
}
\arguments{
\item{x}{The genotyes for only one SNP with dimension n by 1.}

\item{y}{The phentypes with dimension n individuals by K phenotypes.}

\item{cov}{The covariate matrix with dimension n individuals by C covariates.
defalut: there is no covariates.}

\item{cluster.method}{The agglomeration method to be used in the hierarchical
clustering method. This should be one of "ward.D","ward.D2","single",
"complete","average","mcquitty","median","centroid". default: "ward.D2"}
}
\value{
p-value of the test statistic
}
\description{
Hypothesis test for testing the association between K phenotypes with a SNP.
The phenotypes can be either qualitative or binary, espectially the
binary phenotypes with the extremely case-control ratio (the test statistics
has been adjusted by the saddlepoint approximation)
}
\details{
\code{ceCLC()} Modified Computational Efficient Clustering linear combination
 (ceCLC) test

\code{CLC()} Modified Clustering linear combination (CLC) test.

\code{HCLC()} Modified clustering linear combination method based on hierarchical
 clustering (HCLC)

\code{OBrien()} Modified O'Brien test

\code{omnibus()} Modified omnibus (Chi-square) test

\code{TATES()} Modified TATES test

\code{acat.test()} Cauchy combination test

\code{hmp.test()} Harmonic mean p-value

\code{MultiPhen()} MultiPhen test
}
\examples{
N <- 10000
y1 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.999, 0.001)))
y2 <- replicate(2, sample(c(0,1), N, replace = TRUE, prob = c(0.5, 0.5)))
y3 <- replicate(2, rnorm(N))
y <- cbind(y1, y2, y3)
x <- rbinom(N,2,0.3)
ceCLC(x, y)
CLC(x, y)
HCLC(x, y)
OBrien(x, y)
omnibus(x, y)
TATES(x, y)
acat.test(x, y)
hmp.test(x, y)
MultiPhen(x, y)
}
\references{
Wang, M., Zhang, S., & Sha, Q. (2022). A computationally efficient
 clustering linear combination approach to jointly analyze multiple phenotypes
 for GWAS. PloS one, 17(4), e0260911.

Sha, Q., Wang, Z., Zhang, X., & Zhang, S. (2019). A clustering
 linear combination approach to jointly analyze multiple phenotypes for GWAS.
 Bioinformatics, 35(8), 1373-1379.

Li, X., Zhang, S., & Sha, Q. (2020). Joint analysis of multiple
 phenotypes using a clustering linear combination method based on hierarchical
 clustering. Genetic epidemiology, 44(1), 67-78.

Liang, X., Cao, X., Sha, Q., & Zhang, S. (2022). HCLC-FC:
 a novel statistical method for phenome-wide association studies. bioRxiv.

O'Brien, P. C. (1984). Procedures for comparing samples with
 multiple endpoints. Biometrics, 1079-1087.

Van der Sluis, S., Posthuma, D., & Dolan, C. V. (2013).
 TATES: efficient multivariate genotype-phenotype analysis for genome-wide
 association studies. PLoS genetics, 9(1), e1003235.

Liu, Y., Chen, S., Li, Z., Morrison, A. C., Boerwinkle, E., & Lin,
 X. (2019). ACAT: a fast and powerful p value combination method for
 rare-variant analysis in sequencing studies. The American Journal of
 Human Genetics, 104(3), 410-421.

Wilson, D. J. (2019). The harmonic mean p-value for combining
 dependent tests. Proceedings of the National Academy of Sciences, 116(4),
 1195-1200.

O’Reilly, P. F., Hoggart, C. J., Pomyen, Y., Calboli, F. C.,
 Elliott, P., Jarvelin, M. R., & Coin, L. J. (2012). MultiPhen: joint model
 of multiple phenotypes can increase discovery in GWAS. PloS one, 7(5), e34861.
}
