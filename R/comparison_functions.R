compare_2d_map_correlations <- function(map1, map2, size = NULL, cormethod = "pearson"){

  if (is.null(size)){
    size = nrow(map1)
  }

  if (dim(map1)[2] != size){
    map1 <- t(map1)
  }

  if (dim(map2)[2] != size){
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
    correleation_unsorted = cor_mat,
    sorted_order = sorted_order,
    correlation_sorted = cor_mat_sorted
  ))

}





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




evaluate_power_by_effect_magnitude <- function(pvalues, true_effect, thr = 0.05){

  unique_values <- sort(unique(c(true_effect)))
  unique_values <- unique_values[unique_values != 0]

  sig <- pvalues < thr

  results <- tibble()

  for (val in unique_values){

    TPv    <- sum((sig == 1) * (true_effect == val))
    FNv    <- sum((sig == 0) * (true_effect == val))
    Powerv <- TPv / (TPv + FNv)
    results <- rbind(results, tibble(EffectSize = val, TPR = Powerv))

  }

  return(results)

}


