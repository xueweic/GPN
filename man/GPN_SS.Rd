% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/GPN_SS.R
\name{GPN_SS}
\alias{GPN_SS}
\title{Construct a signed bipartite genotype and phenotype network
based on GWAS summary statistics}
\usage{
GPN_SS(z, gwas, well.Net = TRUE, ...)
}
\arguments{
\item{z}{A matrix of z-scores with sign, where positive sign is the protective
effect, negative sign is the risk effect. It defined as z=beta/se(beta).}

\item{gwas}{Alternative summary data giving the estimated effects
(a vector of length p). This, together with \code{shat}, may be
provided instead of \code{z}.}

\item{well.Net}{Logical object. If TRUE, return a well-defined GPN; If FALSE,
return a denser GPN. Default is TRUE.}

\item{\dots}{Arguments passed to the internal function \code{gpn_train}.
Argument "replicates": the number of replicates is used to choose optimal
parameter \code{tau}, default is 50; Argument "weight": Logical obejct,
if TRUE, bipartite projection will include edge weights. Defaults to TRUE.
Argument "metric": a string specifying which entropy to use.
Possible values are: \code{"kl_divergence_snp"},
\code{"kl_divergence_phe"},\code{"kl_divergence_all"},
\code{"cross_entropy_snp"},\code{"cross_entropy_phe"},
\code{"cross_entropy_all"}. The default is
\code{"kl_divergence_all"}. Argument "tuneGrid": the vector with possible
tuning values \code{tau}. It's should between 0 to 1, default value
from 0.05 to 0.95 with a step size of 0.05.}
}
\value{
A list of object type containing: the generated
 random bipartite GPN.
}
\description{
Construct a signed bipartite genotype and phenotype network
based on GWAS summary statistics
}
\examples{
data(gwas_z)
out1 <- GPN_SS(z=gwas_z)
data(gwas_data_all)
out2 <- GPN_SS(gwas=gwas_data_all)
}
