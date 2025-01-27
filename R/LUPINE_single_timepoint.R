#' LUPINE for single time point
#'
#' @param data A 2D array of counts or transformed data with dimensions samples x variables
#' @param lib_size A matrix of library sizes for each sample and time point, optional. Currently if used regression models will be run on log+1 and libsize accounted for
#' @param method The method to use for dimensionality reduction. Options are "pca", "ica", "rpca"
#' @param ncomp The number of components to use for dimensionality reduction
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating the type of parallelisation
#
#' @return A list of correlations and p-values
#' @export
#'
LUPINE_single_timepoint <- function(data_timepoint,
                                    lib_size = NULL,
                                    method = "pca",
                                    ncomp = 1,
                                    BPPARAM = SerialParam()) {

  # Extract variable names
  var_names <- colnames(data_timepoint)
  nVar <- length(var_names)

  # Generate pairwise variable combinations
  len <- nVar * (nVar - 1) / 2
  var1 <- unlist(lapply(1:nVar, function(i) rep(i, (nVar - i))))
  var2 <- unlist(lapply(2:nVar, function(i) seq(i, nVar, 1)))

  # Initialisation
  pcor <- matrix(NA, nrow = nVar, ncol = nVar)
  pcor.pval <- matrix(NA, nrow = nVar, ncol = nVar)
  colnames(pcor) <- colnames(pcor.pval) <- var_names
  rownames(pcor) <- rownames(pcor.pval) <- var_names

  # Run 1 dimensional approximation
  if (method == "pca") {
    loadings_m <- PCA_approx(data_timepoint, ncomp = ncomp)
  } else if (method == "ica") {
    loadings_m <- RPCA_approx(data_timepoint, ncomp = ncomp)
  } else if (method == "rpca") {
    loadings_m <- ICA_approx(data_timepoint, ncomp = ncomp)
  } else {
    stop("Method not supported\n. Use one of pca, ica, rpca.")
  }

  # Apply two indepdent log linear regression models for each combination of variables
  if (is.null(lib_size)) {
    print("No library size detected, continuing without accounting for library size...")
    pcor_list <- bplapply(1:len, function(i) { # apply through each pairwise combination
      loading_tmp <- loadings_m
      loading_tmp[c(var1[i], var2[i]), ] <- 0 # loading vectors with the key variables removed
      u1 <- scale(data_timepoint, center = TRUE, scale = TRUE) %*% loading_tmp
      r_i <- lm(data_timepoint[, var1[i]] ~ u1)
      r_j <- lm(data_timepoint[, var2[i]] ~ u1)
      # partial correlation calculation
      list(
        var1 = var1[i],
        var2 = var2[i],
        pcor_estimate = cor.test(r_i$residuals, r_j$residuals, method = "pearson")$estimate,
        pcor_pval = cor.test(r_i$residuals, r_j$residuals, method = "pearson")$p.value
      )
    }, BPPARAM = BPPARAM)
    for (result in pcor_list) {
      pcor[result$var1, result$var2] <- pcor[result$var2, result$var1] <- result$pcor_estimate
      pcor.pval[result$var1, result$var2] <- pcor.pval[result$var2, result$var1] <- result$pcor_pval
    }
  } else {
    print("Library size detected, accounting for library size...")
    pcor_list <- bplapply(1:len, function(i) { # apply through each pairwise combination
      loading_tmp <- loadings_m
      loading_tmp[c(var1[i], var2[i]), ] <- 0 # loading vectors with the key variables removed
      u1 <- scale(data_timepoint, center = TRUE, scale = TRUE) %*% loading_tmp
      r_i <- lm(log(data_timepoint[, var1[i]] + 1) ~ u1, offset = log(lib_size))
      r_j <- lm(log(data_timepoint[, var2[i]] + 1) ~ u1, offset = log(lib_size))
      # partial correlation calculation
      list(
        var1 = var1[i],
        var2 = var2[i],
        pcor_estimate = cor.test(r_i$residuals, r_j$residuals, method = "pearson")$estimate,
        pcor_pval = cor.test(r_i$residuals, r_j$residuals, method = "pearson")$p.value
      )
    }, BPPARAM = BPPARAM)
    for (result in pcor_list) {
      pcor[result$var1, result$var2] <- pcor[result$var2, result$var1] <- result$pcor_estimate
      pcor.pval[result$var1, result$var2] <- pcor.pval[result$var2, result$var1] <- result$pcor_pval
    }
  }


  pcor.full <- matrix(NA, nrow = nVar, ncol = nVar)
  pcor.pval.full <- matrix(NA, nrow = nVar, ncol = nVar)
  colnames(pcor.full) <- colnames(pcor.pval.full) <- colnames(data_timepoint)
  rownames(pcor.full) <- rownames(pcor.pval.full) <- colnames(data_timepoint)

  common_rows_cols <- intersect(rownames(pcor.full), rownames(pcor))

  pcor.full[common_rows_cols, common_rows_cols] <- pcor[common_rows_cols, common_rows_cols]
  pcor.pval.full[common_rows_cols, common_rows_cols] <- pcor.pval[common_rows_cols, common_rows_cols]

  res <- list(Estimate = pcor.full, pvalue = pcor.pval.full)

  return(res)

}
