

#' Degree of the bipartite network
#'
#' @aliases degree
#' @param Net A adjacency matrix with genotype or genes as rows,
#'  phenotype as columns and the association as entries.
#'  The adjacency matrix of the bipartite network which
#'  could be a simple matrix or data frame.
#' @param weight Logical; if TRUE, bipartite projection
#'  will include edge weights. Defaults to TRUE.
#'
#' @return A list of object type containing: node degree
#'  and standardized node degree
#' @export
#' @keywords GPN
#' @examples
#' set.seed(123)
#' Net <- matrix(rnorm(2000)^2, nrow = 500, ncol = 4)
#' out <- degree(Net)
#'



degree <- function(Net, weight = TRUE) {
  if (!is.logical(weight)) {
    stop(paste("Weight", weight, "should be logical value"), call. = FALSE)
  }
  if (!is.matrix(Net))
    Net <- as.matrix(Net)
  Net <- abs(Net)
  if (weight) {
    d_snp <- rowSums(Net)
    d_phe <- colSums(Net)
  } else {
    d_snp <- rowSums(Net != 0)
    d_phe <- colSums(Net != 0)
  }
  # standardized
  d_snp_sd <- (d_snp - min(d_snp))/(max(d_snp) - min(d_snp))
  d_phe_sd <- (d_phe - min(d_phe))/(max(d_phe) - min(d_phe))
  list(degree_SNP = d_snp, degree_phenotype = d_phe, degree_SNP_sd = d_snp_sd,
       degree_phenotype_sd = d_phe_sd)
}


### internal function
## \dots Arguments passed to the internal function, where default ratio <- 1e-10

degree_shrink <- function(Net, weight = TRUE, ...) {
  # shrink_method=c('all','partial')
  if (!is.logical(weight)) {
    stop(paste("Weight", weight, "should be logical value"), call. = FALSE)
  }

  if (!is.matrix(Net))
    Net <- as.matrix(Net)
  Net <- abs(Net)

  additional_args <- list(...)
  if ("ratio" %in% names(additional_args)) {
    ratio <- additional_args$ratio
  } else {
    ratio <- 1e-10
  }
  # shrink, using the ratio to eliminate 0 standard degree
  if (weight) {
    d_snp <- rowSums(Net) + ratio
    d_phe <- colSums(Net) + ratio
  } else {
    d_snp <- rowSums(Net != 0) + ratio
    d_phe <- colSums(Net != 0) + ratio
  }
  # standardized
  d_snp_sd <- (d_snp - min(d_snp))/(max(d_snp) - min(d_snp) + ratio)
  d_phe_sd <- (d_phe - min(d_phe))/(max(d_phe) - min(d_phe) + ratio)
  d_snp_sd[which(d_snp_sd >= 1)] <- 1 - ratio
  d_snp_sd[which(d_snp_sd == 0)] <- ratio
  d_phe_sd[which(d_phe_sd >= 1)] <- 1 - ratio
  d_phe_sd[which(d_phe_sd == 0)] <- ratio

  list(degree_SNP_shrink = d_snp, degree_phenotype_shrink = d_phe, degree_SNP_sd_shrink = d_snp_sd,
       degree_phenotype_sd_shrink = d_phe_sd)

}

