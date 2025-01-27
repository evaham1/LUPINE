#' LUPINE for single time point
#'
#' @param data A 2D array of counts or transformed data with dimensions samples x variables
#' @param lib_size A matrix of library sizes for each sample and time point, optional. 
#'        If used, regression models will be run on log-transformed data with library size as an offset.
#' @param method The method to use for dimensionality reduction. Options are "pca", "ica", "rpca"
#' @param ncomp The number of components to use for dimensionality reduction
#' @param BPPARAM A \linkS4class{BiocParallelParam} object indicating the type of parallelisation
#'
#' @return A list of correlations and p-values
#' @export
#'
LUPINE_single_timepoint <- function(data_timepoint,
                                    lib_size = NULL,
                                    method = "pca",
                                    ncomp = 1,
                                    BPPARAM = SerialParam()) {

  # Extract variable names and calculate the number of variables
  var_names <- colnames(data_timepoint)
  nVar <- length(var_names)

  # Generate all pairwise combinations of variables for correlation analysis
  len <- nVar * (nVar - 1) / 2
  var1 <- unlist(lapply(1:nVar, function(i) rep(i, (nVar - i))))
  var2 <- unlist(lapply(2:nVar, function(i) seq(i, nVar, 1)))

  # Initialize matrices to store partial correlation coefficients and p-values
  pcor <- matrix(NA, nrow = nVar, ncol = nVar)
  pcor.pval <- matrix(NA, nrow = nVar, ncol = nVar)
  colnames(pcor) <- colnames(pcor.pval) <- var_names
  rownames(pcor) <- rownames(pcor.pval) <- var_names

  # Perform dimensionality reduction using the specified method
  if (method == "pca") {
    loadings_m <- PCA_approx(data_timepoint, ncomp = ncomp)
  } else if (method == "ica") {
    loadings_m <- RPCA_approx(data_timepoint, ncomp = ncomp)
  } else if (method == "rpca") {
    loadings_m <- ICA_approx(data_timepoint, ncomp = ncomp)
  } else {
    stop("Method not supported. Use one of 'pca', 'ica', 'rpca'.")
  }

  # Define a function to compute partial correlations for a single pair
  compute_partial_corr <- function(i) {
    loading_tmp <- loadings_m
    loading_tmp[c(var1[i], var2[i]), ] <- 0  # Remove loading vectors of the current variable pair
    u1 <- scale(data_timepoint, center = TRUE, scale = TRUE) %*% loading_tmp  # Dimensionality-reduced data

    # Fit linear models for both variables - may account for library size / log normalise
    r_i <- if (is.null(lib_size)) {
      lm(data_timepoint[, var1[i]] ~ u1)
    } else {
      lm(log(data_timepoint[, var1[i]] + 1) ~ u1, offset = log(lib_size))
    }
    r_j <- if (is.null(lib_size)) {
      lm(data_timepoint[, var2[i]] ~ u1)
    } else {
      lm(log(data_timepoint[, var2[i]] + 1) ~ u1, offset = log(lib_size))
    }

    # Calculate partial correlation and p-value
    list(
      var1 = var1[i],
      var2 = var2[i],
      pcor_estimate = cor.test(r_i$residuals, r_j$residuals, method = "pearson")$estimate,
      pcor_pval = cor.test(r_i$residuals, r_j$residuals, method = "pearson")$p.value
    )
  }

  # Apply the function in parallel to all variable pairs
  pcor_list <- bplapply(1:len, compute_partial_corr, BPPARAM = BPPARAM)

  # Populate the partial correlation and p-value matrices
  for (result in pcor_list) {
    pcor[result$var1, result$var2] <- pcor[result$var2, result$var1] <- result$pcor_estimate
    pcor.pval[result$var1, result$var2] <- pcor.pval[result$var2, result$var1] <- result$pcor_pval
  }

  # Create full matrices for partial correlations and p-values
  pcor.full <- matrix(NA, nrow = nVar, ncol = nVar)
  pcor.pval.full <- matrix(NA, nrow = nVar, ncol = nVar)
  colnames(pcor.full) <- colnames(pcor.pval.full) <- colnames(data_timepoint)
  rownames(pcor.full) <- rownames(pcor.pval.full) <- colnames(data_timepoint)

  # Populate the full matrices with the calculated partial correlations
  common_rows_cols <- intersect(rownames(pcor.full), rownames(pcor))
  pcor.full[common_rows_cols, common_rows_cols] <- pcor[common_rows_cols, common_rows_cols]
  pcor.pval.full[common_rows_cols, common_rows_cols] <- pcor.pval[common_rows_cols, common_rows_cols]

  # Return the results as a list
  res <- list(Estimate = pcor.full, pvalue = pcor.pval.full)
  return(res)
}
