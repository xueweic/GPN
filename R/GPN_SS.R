

#' Construct a signed bipartite genotype and phenotype network
#' based on GWAS summary statistics
#'
#' @param z A matrix of z-scores with sign, where positive sign is the protective
#'  effect, negative sign is the risk effect. It defined as z=beta/se(beta).
#'
#' @param gwas Alternative summary data giving the estimated effects
#'   (a vector of length p). This, together with \code{shat}, may be
#'   provided instead of \code{z}.
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
#'
#' @param well.Net Logical object. If TRUE, return a well-defined GPN; If FALSE,
#'  return a denser GPN. Default is TRUE.
#'
#' @return
#' @import stats
#' @export
#'
#' @examples
#' data(gwas_z)
#' out1 <- GPN_SS(z=gwas_z)
#' data(gwas_data_all)
#' out2 <- GPN_SS(gwas=gwas_data_all)



GPN_SS <- function(z, gwas, well.Net = TRUE, ...) {
  # Check inputs z, or gwas
  if (sum(c(missing(z), missing(gwas))) != 1)
    stop("Please provide either z or gwas, but not both")
  if (missing(gwas) & !missing(gwas)) {
    z[is.na(z)] = 0
  }
  if (missing(z) & !missing(gwas)) {
    snp.set <- as.vector(sapply(1:length(gwas), function(x) {
      gwas[[x]]$SNP
    }))
    snp.set <- unique(snp.set)
    M <- length(snp.set)
    K <- length(gwas)
    z <- matrix(NA, nrow = M, ncol = K)
    for (i in 1:length(gwas)) {
      temp <- sign(gwas[[i]]$Beta) * qnorm(gwas[[i]]$P/2, lower.tail = FALSE)
      # temp <- gwas[[i]]$Beta/gwas[[i]]$SE
      ind <- match(gwas[[i]]$SNP, snp.set)
      z[ind, i] <- temp
    }
    z[is.na(z)] = 0
  }
  Net <- sign(z) * (z^2)
  if (well.Net == TRUE) {
    out <- gpn_train(Net, ...)
  } else {
    out <- Net
  }
  return(out)
}
