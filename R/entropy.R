

#' Bipartite Network Entropy
#'
#' Various functions for calculating the Bipartite network entropy
#'
#' These functions can be used to calculate the entropy of the GPN.
#'
#' \code{degree_entropy} is one of the Shannon entropy. The Degree entropy
#' of GPN can be used to measure the diversity of associations
#' between genetic variants and phenotypes.
#'
#' \code{cross_entropy} It is used to determine the diversity
#' between a bipartite GPN and a random network.
#'
#' \code{kl_divergence} It is used to determine the diversity
#' between a bipartite GPN and a random network.
#'
#' @aliases degree_entropy cross_entropy kl_divergence
#' @param Net A adjacency matrix with genotype or genes as rows,
#'  phenotype as columns and the association as entries.
#'  The adjacency matrix of the bipartite network which
#'  could be a simple matrix or data frame.
#' @param Net.r It is a adjacency matrix of the random GPN in which
#'  the degrees of all vertices are the same as the original
#'  bipartite GPN. The adjacency matrix could be
#'  a simple matrix or data frame.
#' @param weight Logical; if TRUE, bipartite projection
#'  will include edge weights. Defaults to TRUE.
#' @param \dots Arguments passed to the internal function
#'  \code{degree_shrink} to avoid no definition of entropy.
#'  The default is \code{ratio = 1e-10}.
#' @keywords GPN
#' @return A list of object type containing: first value
#'  represents the degree entropy of genetic variants/genes;
#'  second value represents the degree entropy of phenotype;
#'  the third value represents the degree entropy of whole
#'  bipartite network.
#' @references Kullback, S., & Leibler, R. A. (1951).
#'  \emph{On information and sufficiency}. The annals of mathematical
#'  statistics, 22(1), 79-86.
#' @references Murphy, K. P. (2012).\emph{ Machine learning:
#'  a probabilistic perspective}. MIT press.
#' @examples
#' set.seed(123)
#' Net <- matrix(rnorm(2000)^2, nrow = 500, ncol = 4)
#' Net.r <- GPN_random(Net)$Net.R
#' out <- degree_entropy(Net, weight = TRUE)
#' out2 <- cross_entropy(Net, Net.r)
#' out3 <- kl_divergence(Net, Net.r)
#' @export degree_entropy

degree_entropy <- function(Net, weight = TRUE, ...) {
  if (!is.logical(weight)) {
    stop(paste("Weight", weight, "should be logical value"), call. = FALSE)
  }
  if (!is.matrix(Net))
    Net <- as.matrix(Net)

  d <- degree_shrink(Net, weight, ...)
  d.snp <- d$degree_SNP_sd
  d.phe <- d$degree_phenotype_sd
  H.snp <- -d.snp %*% log(d.snp)
  H.phe <- -d.phe %*% log(d.phe)
  H <- H.snp + H.phe
  list(degree_entropy_snp = as.numeric(H.snp),
       degree_entropy_phe = as.numeric(H.phe),
       degree_entropy_all = as.numeric(H))
}


#' @rdname entropy
#' @export
cross_entropy <- function(Net, Net.r, weight = TRUE, ...) {
  if (!is.logical(weight)) {
    stop(paste("Weight", weight, "should be logical value"), call. = FALSE)
  }

  if (!is.matrix(Net))
    Net <- as.matrix(Net)

  d <- degree_shrink(Net, weight, ...)
  dr <- degree_shrink(Net.r, weight, ...)

  d.snp <- d$degree_SNP_sd_shrink
  d.phe <- d$degree_phenotype_sd_shrink
  dr.snp <- dr$degree_SNP_sd_shrink
  dr.phe <- dr$degree_phenotype_sd_shrink

  H.snp <- -sum(d.snp * log(dr.snp))
  H.phe <- -sum(d.phe * log(dr.phe))
  H <- H.snp + H.phe
  list(cross_entropy_snp = as.numeric(H.snp),
       cross_entropy_phe = as.numeric(H.phe),
       cross_entropy_all = as.numeric(H))
}


#' @rdname entropy
#' @export
kl_divergence <- function(Net, Net.r, weight = TRUE, ...) {
  if (!is.logical(weight)) {
    stop(paste("Weight", weight, "should be logical value"), call. = FALSE)
  }

  if (!is.matrix(Net))
    Net <- as.matrix(Net)

  d <- degree_shrink(Net, weight, ...)
  dr <- degree_shrink(Net.r, weight, ...)
  d.snp <- d$degree_SNP_shrink
  d.phe <- d$degree_phenotype_shrink
  dr.snp <- dr$degree_SNP_shrink
  dr.phe <- dr$degree_phenotype_shrink

  r1 <- d.snp/dr.snp
  r2 <- d.phe/dr.phe
  KL.SNP <- d.snp %*% log(r1)
  KL.phe <- d.phe %*% log(r2)
  KL.all <- KL.SNP + KL.phe

  list(kl_divergence_snp = as.numeric(KL.SNP),
       kl_divergence_phe = as.numeric(KL.phe),
       kl_divergence_all = as.numeric(KL.all))
}
