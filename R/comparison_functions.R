#' Compare correlations of two sets of spatial maps
#'
#' @details
#' `compare_2d_map_correlations` takes two sets of spatial maps as arguments and
#' returns the correlation between them as well as information about sorting
#' them to be as aligned as possible.
#'
#' @param map1 A Q x V matrix, with each row corresponding to a component
#' @param map2 A Q x V matrix, with each row corresponding to a component
#' to the number of independent components that you plan to include in the model.
#' @param Q Optional parameter giving the number of components This will be used
#' to fixed the orientation of map1 and map2 if needed
#' @param cormethod technique used to calculate the correlation. Default is "pearson"
#' @return
#' \itemize{
#'   \item correlation_unsorted - The raw correlation matrix between the two sets
#' of maps
#'   \item sorted_order - A vector that can be used to sort map1 to be most aligned with map2
#'   \item correlation_sorted - The correlation between the two sets of maps after map1 has been sorted.
#' For this result, the largest correlation should be along the diagonal.
#' }
#'
#' @export
compare_2d_map_correlations <- function(map1, map2, Q = NULL, cormethod = "pearson"){

  if (is.null(Q)){
    Q = nrow(map1)
  }

  if (dim(map1)[2] != Q){
    map1 <- t(map1)
  }

  if (dim(map2)[2] != Q){
    map2 <- t(map2)
  }

  cor_mat <- cor(map1, map2, method = cormethod)

  # Get the sorted correlations
  abs_cor <- abs(cor_mat)
  sorted_order <- rep(0, size)
  for (i in 1:size){
    maxcor = which(abs_cor == max(abs_cor), arr.ind = T)[1,]
    sorted_order[maxcor[2]] =  maxcor[1]
    abs_cor[,maxcor[2]] = 0
    abs_cor[maxcor[1],] = 0
  }

  cor_mat_sorted <- cor(map1[,sorted_order], map2, method = cormethod)

  return(list(
    correlation_unsorted = cor_mat,
    sorted_order = sorted_order,
    correlation_sorted = cor_mat_sorted
  ))

}


#' Get TPR, FDR, etc
#'
#' @param pvalues p-values for the method of interest
#' @param true_effect_indicators Matrix with 1 where truth is non-zero and 0 otherwise
#' @return A list containing standard evaluation quantities
#'
# true effect indicators should be all 0 (not sig) or 1 (sig)
evaluate_effect_discovery <- function(pvalues, true_effect_indicators){

  # By specific effect/IC combination
  nThr <- 101
  threshold <- seq(from = 0.0, to = 1.0, length.out = nThr)
  TP_all <- NULL
  FP_all <- NULL
  FN_all <- NULL
  TN_all <- NULL

  for (thr_level in threshold){
    sig <- pvalues < thr_level
    if (thr_level == 1.0){sig <- pvalues <= thr_level }
    TPi = sum((sig == 1) * (true_effect_indicators == 1))
    FPi = sum((sig == 1) * (true_effect_indicators == 0))
    FNi = sum((sig == 0) * (true_effect_indicators == 1))
    TNi = sum((sig == 0) * (true_effect_indicators == 0))

    TP_all <- c(TP_all, TPi)
    FP_all <- c(FP_all, FPi)
    FN_all <- c(FN_all, FNi)
    TN_all <- c(TN_all, TNi)
  }

  precision <- TP_all / (TP_all + FP_all)
  precision[is.nan(precision)] = 1.0 # no discoveries
  recall    <- TP_all / (TP_all + FN_all)
  TPR <- recall
  FPR <- FP_all / (FP_all + TN_all)
  FDR <- FP_all / (FP_all + TP_all)
  FDR[is.nan(FDR)] = 0.0

  results <- as_tibble(cbind(threshold, TP_all, FP_all, FN_all, TN_all, precision, recall, TPR, FPR, FDR))

  return(results)

}


