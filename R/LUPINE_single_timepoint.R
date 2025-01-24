#' LUPINE for single time point
#'
#' @param data A 2D array of counts or transformed data with dimensions samples x taxa
#' @param lib_size A matrix of library sizes for each sample and time point, optional. Currently if used regression models will be run on log+1 and libsize accounted for
#' @param method The method to use for dimensionality reduction. Options are "pca", "ica", "rpca"
#' @param ncomp The number of components to use for dimensionality reduction
#'
#' @return A list of correlations and p-values
#' @export
#'
LUPINE_single_timepoint <- function(data_timepoint,
                                    lib_size = NULL,
                                    method = "pca",
                                    ncomp = 1) {

  # Extract variable names
  taxa_names <- colnames(data_timepoint)
  nVar <- length(taxa_names)

  # Generate pairwise variable combinations
  len <- nVar * (nVar - 1) / 2
  taxa1 <- unlist(lapply(1:nVar, function(i) rep(i, (nVar - i))))
  taxa2 <- unlist(lapply(2:nVar, function(i) seq(i, nVar, 1)))

  # Initialisation
  pcor <- matrix(NA, nrow = nVar, ncol = nVar)
  pcor.pval <- matrix(NA, nrow = nVar, ncol = nVar)
  colnames(pcor) <- colnames(pcor.pval) <- taxa_names
  rownames(pcor) <- rownames(pcor.pval) <- taxa_names

  # This should be calculated outside the LUPINE function
  # Creating library size by summing taxa counts per sample and time point
  # if (is.null(lib_size) & !is.transformed) {
  #   lib_size <- apply(data, c(1, 3), sum)
  # }

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
  for (i in 1:len) { # loop through each pairwise combination
    # print(i)
    loading_tmp <- loadings_m
    loading_tmp[c(taxa1[i], taxa2[i]), ] <- 0 # loading vectors with the key variables removed
    u1 <- scale(data_timepoint, center = TRUE, scale = TRUE) %*% loading_tmp
    if (is.null(lib_size)) {
      # principal component regression directly on counts
      print("No library size detected, continuing without accounting for library size...")
      r_i <- lm(data_timepoint[, taxa1[i]] ~ u1)
      r_j <- lm(data_timepoint[, taxa2[i]] ~ u1)
    } else {
      # principal component regression on log counts+1 with an offset for library size
      print("Library size detected, accounting for library size...")
      r_i <- lm(log(data_timepoint[, taxa1[i]] + 1) ~ u1, offset = log(lib_size))
      r_j <- lm(log(data_timepoint[, taxa2[i]] + 1) ~ u1, offset = log(lib_size))
    }

    # partial correlation calculation
    pcor[taxa1[i], taxa2[i]] <- pcor[taxa2[i], taxa1[i]] <- cor.test(r_i$residuals,
                                                                     r_j$residuals,
                                                                     method = "pearson"
    )$estimate
    pcor.pval[taxa1[i], taxa2[i]] <- pcor.pval[taxa2[i], taxa1[i]] <- cor.test(r_i$residuals,
                                                                               r_j$residuals,
                                                                               method = "pearson"
    )$p.value
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
