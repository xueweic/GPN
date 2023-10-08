

#' Adjust the Pvalue Matrix for Bipartite Network
#'
#' `adjust_p()` : Given a set of p-values, returns p-values adjusted
#'  using one of several methods.
#'
#' @param p The p-value matrix is an object where SNPs are in rows and
#'  phentoypes are in columns. The entries represent the strength of
#'  association between SNP and phenotype. This could be a simple matrix
#'  or data frame.
#' @param method A string that specifies what methods will
#'  be used to adjust p-value, which can be "qvalue","lfdr (local FDR)",
#'  or "BH (Benjamini & Hochberg (1995))". The "BH" method is time consuming
#'  for the huge matrix. The default method is "qvalue".
#' @param \dots Arguments passed to the \code{qvalue()}.
#' @return A matrix of adjusted p-value
#' @importFrom qvalue qvalue
#' @importFrom stats p.adjust.methods
#' @importFrom stats p.adjust
#' @export
#' @keywords GPN
#' @examples
#' set.seed(123)
#' Net <- matrix(rnorm(2000)^2, nrow = 500, ncol = 4)
#' p <- get_pvalue(Net)$p
#' out <- adjust_p(p, method = "qvalue")




adjust_p <- function(p, method = "qvalue", ...) {
  allMethods <- c("qvalue", "lfdr", "BH")
  if (!(method %in% allMethods))
    stop(paste("p-value adjust methods are limited to:", paste(allMethods,
                                                        collapse = ", ")), call. = FALSE)

  p <- as.matrix(p)
  re <- p
  for (ind in 1:ncol(p)) {
    ind.K <- which(p[, ind] != 1)
    if (method %in% c("qvalue") || is.null(method)) {
      re[ind.K, ind] <- qvalue(p[ind.K, ind], ...)$qvalues
    } else if (method == "lfdr") {
      re[ind.K, ind] <- qvalue(p[ind.K, ind], ...)$lfdr
    } else if (method %in% p.adjust.methods) {
      #'BH'
      re <- p.adjust(p, method = method)
    } else {
      stop("Please specific the adjust p-value method!")
    }
  }
  out <- list(adjust_p = re)
  out
}




#' @param Net A adjacency matrix with genotype or genes as rows,
#'  phenotype as columns and the association as entries.
#'  The adjacency matrix of the bipartite network which
#'  could be a simple matrix or data frame.
#' @return a list object
#'
#' @importFrom stats pchisq
get_pvalue <- function(Net) {
  list(p = pchisq(Net, df = 1, lower.tail = FALSE))
}
