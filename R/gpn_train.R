

#' Bipartite GPN over Different Tuning Parameters
#'
#' @param Net A adjacency matrix with genotype or genes as rows,
#'  phenotype as columns and the association as entries.
#'  The adjacency matrix of the bipartite network which
#'  could be a simple matrix or data frame.
#' @param replicates The number of replicates is used to
#'  choose optimal parameter \code{tau}, default is 50.
#' @param weight Logical; if TRUE, bipartite projection
#'  will include edge weights. Defaults to TRUE.
#' @param metric A string specifying which entropy to use.
#'  Possible values are: \code{"kl_divergence_snp"},
#'  \code{"kl_divergence_phe"},\code{"kl_divergence_all"},
#'  \code{"cross_entropy_snp"},\code{"cross_entropy_phe"},
#'  \code{"cross_entropy_all"}. The default is
#'  \code{"kl_divergence_all"}.
#' @param tuneGrid The vector with possible tuning values \code{tau}.
#'  it's should between 0 to 1, default argument from 0.05 to 0.95
#'  with a step size of 0.05
#' @param \dots Arguments passed to the internal function
#'  \code{degree_shrink} to avoid no definition of entropy.
#'  The default is \code{ratio = 1e-10}.
#' @return A list of object type containing:
#'  \item{modelIonfo }{The original Bipartite GPN and
#'  predefined parameters.} \item{parameter}{A vector
#'  contains all tuning parameters.} \item{metric}{the method
#'  is used to choose the well-defined bipartite GPN.}
#'  \item{entropy}{A list consists of degree entropy,
#'  cross-entropy, and kl_divergence over different
#'  tuning parameters.} \item{perfEntropy}{The performance
#'  of different metric.} \item{bestTune}{The best tuning parameter
#'  corresponds to the well-define bipartite GPN.}
#'  \item{finalModel}{The final output includes the denser and
#'  well-defined adjacency matrix of GPN, and the network index,
#'  such as degree, betweenness, cluster,etc. }
#' @import foreach
#' @keywords GPN
#'
#' @examples
#' set.seed(123)
#' Net <- matrix(rnorm(2000)^2, nrow = 500, ncol = 4)
#' out <- gpn_train(Net,replicates=10,metric = "kl_divergence_snp",tuneGrid=seq(0.1,0.9,0.1))
#'


gpn_train <- function(Net, replicates = 50, weight = TRUE,
                     metric = "kl_divergence_snp",
                     tuneGrid = seq(0.05,0.95,0.05), ...) {

  startTime <- proc.time()

  if (!is.matrix(Net)) {
    Net <- as.matrix(Net)
  }


  if (!is.logical(weight)) {
    stop(paste("Weight", weight, "should be logical value"), call. = FALSE)
  }

  if (!is.numeric(replicates)) {
    stop(paste("replicates", replicates, "should be positive interger"),
         call. = FALSE)
  }

  metric.all <- c("kl_divergence_snp", "kl_divergence_phe", "kl_divergence_all",
                  "cross_entropy_snp", "cross_entropy_phe", "cross_entropy_all")
  if (!metric %in% metric.all) {
    stop(paste("Metric", metric, "not applicable for this models"),
         call. = FALSE)
  }


  if (is.null(tuneGrid) || !is.vector(tuneGrid) || max(tuneGrid) >= 1 ||
      min(tuneGrid) <= 0) {
    stop(paste("tuneGrid", tuneGrid, "not applicable for this models"),
         call. = FALSE)
  }

  ## Remove duplicates from grid
  tuneGrid <- tuneGrid[!duplicated(tuneGrid), drop = FALSE]
  tuneGrid <- sort(tuneGrid, decreasing = FALSE)

  # save the result
  res.DE <- data.frame(degree_entropy_snp = rep(NA, 0),
                       degree_entropy_phe = rep(NA,0),
                       degree_entropy_all = rep(NA, 0))
  res.CE <- data.frame(cross_entropy_snp = rep(NA, 0),
                       cross_entropy_phe = rep(NA, 0),
                       cross_entropy_all = rep(NA, 0))
  res.KL <- data.frame(kl_divergence_snp = rep(NA, 0),
                       kl_divergence_phe = rep(NA, 0),
                       kl_divergence_all = rep(NA, 0))
  # calculate the index for original GPN
  p <- get_pvalue(Net)$p
  modelInfo <- list(Net = Net, pvalue = p, metric = metric, replicates = replicates,
                    weight = weight, tuneGrid = tuneGrid)
  # the first row is the original DE, the remaining rows are random DE
  tmp.DE <- as.data.frame(degree_entropy(Net, weight, ...))
  res.DE <- rbind(res.DE, tmp.DE)

  # replicates 1000 times, generate the random network, find optimal tau
  init_result <- lapply(1:length(tuneGrid), function(i) {
    Net.tau <- adjust_adjacency(Net, p, tau = tuneGrid[i], method = "qvalue", ...)$Sparsity.Net
    foreach(rep = 1:replicates, .combine = "c", .verbose = FALSE, .errorhandling = "stop") %do%
      {
        Net.r <- GPN_random(Net.tau)[[1]]
        tmp.DE <- as.data.frame(degree_entropy(Net.tau, weight, ...))
        res.DE <- rbind(res.DE, tmp.DE)
        tmp.CE <- as.data.frame(cross_entropy(Net.tau, Net.r, weight, ...))
        res.CE <- rbind(res.CE, tmp.CE)
        tmp.KL <- as.data.frame(kl_divergence(Net.tau, Net.r, weight, ...))
        res.KL <- rbind(res.KL, tmp.KL)
      }
    list(degree_entropy = res.DE, cross_entropy = res.CE, kl_divergence = res.KL)
  })
  names(init_result) <- sapply(1:length(tuneGrid), function(i) {
    paste("tau_", tuneGrid[i])
  })

  # find the optimal parameter
  perfEntropy <- sapply(1:length(tuneGrid), function(i) {
    c(colMeans(init_result[[i]]$kl_divergence, dims = 1),
      colMeans(init_result[[i]]$cross_entropy, dims = 1))
  })
  perfEntropy <- as.data.frame(perfEntropy)
  colnames(perfEntropy) <- sapply(1:length(tuneGrid), function(i) {
    paste("tau_", tuneGrid[i])
  })

  bestTune <- lapply(1:nrow(perfEntropy), function(i) {
    tuneGrid[which.max(perfEntropy[i, ])]
  })

  names(bestTune) <- c("kl_divergence_snp", "kl_divergence_phe", "kl_divergence_all",
                       "cross_entropy_snp", "cross_entropy_phe", "cross_entropy_all")

  best.tau <- bestTune[[which(metric.all == metric)]]
  Net.well <- adjust_adjacency(Net, p, tau = best.tau, method = "qvalue",
                             ...)$Sparsity.Net

  # calculate the index for original network and well-defined network
  # finalModel <- list(metric = metric, tau = best.tau, Net = Net, pvalue = p,
  #                    Net.well = Net.well, denser = list(connectance = connectance(Net),
  #                                                       betweenness = betweenness(Net),
  #                                                       cluster = cluster(Net),
  #                                                       degree = degree(Net, weight)),
  #                    well.defined = list(connectance = connectance(Net.well),
  #                                        betweenness = betweenness(Net.well),
  #                                        cluster = cluster(Net.well),
  #                                        degree = degree(Net.well, weight)))
  finalModel <- list(metric = metric, tau = best.tau, Net = Net, pvalue = p,
                     Net.well = Net.well, denser = list(connectance = connectance(Net),
                                                        betweenness = betweenness(Net),
                                                        degree = degree(Net, weight)),
                     well.defined = list(connectance = connectance(Net.well),
                                         betweenness = betweenness(Net.well),
                                         degree = degree(Net.well, weight)))

  endTime <- proc.time()
  times <- endTime - startTime
  out <- structure(list(metric = metric, tau = best.tau,
                        Net.well = Net.well,
                        Net.topology = list(connectance = connectance(Net.well),
                                            betweenness = betweenness(Net.well),
                                            degree = degree(Net.well, weight)),
                        times = times))
  out
}
