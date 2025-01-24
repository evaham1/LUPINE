#' LUPINE for single time point
#'
#' @param data A 3D array of counts or transformed data with dimensions samples x taxa x time points
#' @param timepoint The time point for which the correlations are calculated
#' @param excluded_taxa A list of taxa to be excluded at each time point
#' @param is.transformed A logical indicating whether the data is transformed or not
#' @param lib_size A matrix of library sizes for each sample and time point
#' @param method The method to use for dimensionality reduction. Options are "pca", "ica", "rpca"
#' @param ncomp The number of components to use for dimensionality reduction
#'
#' @return A list of correlations and p-values
#' @export
#'
LUPINE_single <- function(data,
                          timepoint,
                          excluded_taxa = NULL,
                          is.transformed = FALSE,
                          lib_size = NULL,
                          method = "pca",
                          ncomp = 1) {

  # Extract total number of variables (p)
  nVar_total <- dim(data)[2]

  # Extract matrix for the current time point
  data_timepoint <- data[, , timepoint]

  # Extract variable names after excluded variables
  taxa_names <- colnames(data_timepoint)[!(colnames(data_timepoint) %in% excluded_taxa[[timepoint]])]

  # pairwise variable combinations
  nVar <- length(taxa_names)
  len <- nVar * (nVar - 1) / 2
  taxa1 <- unlist(lapply(1:nVar, function(i) rep(i, (nVar - i))))
  taxa2 <- unlist(lapply(2:nVar, function(i) seq(i, nVar, 1)))

  # Initialisation
  pcor <- matrix(NA, nrow = nVar, ncol = nVar)
  pcor.pval <- matrix(NA, nrow = nVar, ncol = nVar)
  colnames(pcor) <- colnames(pcor.pval) <- taxa_names
  rownames(pcor) <- rownames(pcor.pval) <- taxa_names

  if (is.null(lib_size) & !is.transformed) {
    # Creating library size by summing taxa counts per sample and time point
    lib_size <- apply(data, c(1, 3), sum)
  }
  # Extract count array after excluding taxa
  data_filt <- data[, taxa_names, ]
  if(length(dim(data_filt))> 2) {
    data_timepoint_f <- data_filt[, , timepoint]
  } else {
    data_timepoint_f <- data_filt
  }

  if (method == "pca") {
    loadings_m <- PCA_approx(data_timepoint_f, ncomp = ncomp)
  } else if (method == "ica") {
    loadings_m <- RPCA_approx(data_timepoint_f, ncomp = ncomp)
  } else if (method == "rpca") {
    loadings_m <- ICA_approx(data_timepoint_f, ncomp = ncomp)
  } else {
    stop("Method not supported\n. Use one of pca, ica, rpca.")
  }

  for (i in 1:len) {
    # print(i)
    loading_tmp <- loadings_m
    loading_tmp[c(taxa1[i], taxa2[i]), ] <- 0
    u1 <- scale(data_timepoint_f, center = TRUE, scale = TRUE) %*% loading_tmp
    if (!is.transformed) {
      ##** LUPINE_single with counts**##
      # principal component regression on log counts+1 with an offset for library size
      r_i <- lm(log(data_timepoint_f[, taxa1[i]] + 1) ~ u1, offset = log(lib_size[, timepoint]))
      r_j <- lm(log(data_timepoint_f[, taxa2[i]] + 1) ~ u1, offset = log(lib_size[, timepoint]))
    } else {
      ##** LUPINE_single with clr**##
      # principal component regression on clr
      r_i <- lm(data_timepoint_f[, taxa1[i]] ~ u1)
      r_j <- lm(data_timepoint_f[, taxa2[i]] ~ u1)
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

  pcor.full <- matrix(NA, nrow = nVar_total, ncol = nVar_total)
  pcor.pval.full <- matrix(NA, nrow = nVar_total, ncol = nVar_total)
  colnames(pcor.full) <- colnames(pcor.pval.full) <- colnames(data_timepoint)
  rownames(pcor.full) <- rownames(pcor.pval.full) <- colnames(data_timepoint)

  common_rows_cols <- intersect(rownames(pcor.full), rownames(pcor))

  pcor.full[common_rows_cols, common_rows_cols] <- pcor[common_rows_cols, common_rows_cols]
  pcor.pval.full[common_rows_cols, common_rows_cols] <- pcor.pval[common_rows_cols, common_rows_cols]

  res <- list(Estimate = pcor.full, pvalue = pcor.pval.full)

  return(res)
}
