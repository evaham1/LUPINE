#' LUPINE for single time point
#'
#' @param data A 3D array of counts or transformed data with dimensions samples x taxa x time points
#' @param day The time point for which the correlations are calculated
#' @param excluded_taxa A list of taxa to be excluded at each time point
#' @param is.transformed A logical indicating whether the data is transformed or not
#' @param lib_size A matrix of library sizes for each sample and time point
#' @param method The method to use for dimensionality reduction. Options are "pca", "ica", "rpca"
#' @param ncomp The number of components to use for dimensionality reduction
#'
#' @return A list of correlations and p-values
#' @export
#'
LUPINE_single <- function(data, day, excluded_taxa = NULL, is.transformed = FALSE, lib_size = NULL,
                          method = "pca", ncomp = 1) {
  # Total number of taxa (p)
  nOTU_total <- dim(data)[2]
  # Extract matrix for current time point
  data_day <- data[, , day]
  # Extract taxa names after excluded taxa
  taxa_names <- colnames(data_day)[!(colnames(data_day) %in% excluded_taxa[[day]])]
  # pairwise taxa combinations
  nOTU <- length(taxa_names)
  len <- nOTU * (nOTU - 1) / 2
  taxa1 <- unlist(lapply(1:nOTU, function(i) rep(i, (nOTU - i))))
  taxa2 <- unlist(lapply(2:nOTU, function(i) seq(i, nOTU, 1)))

  # Initialisation
  pcor <- matrix(NA, nrow = nOTU, ncol = nOTU)
  pcor.pval <- matrix(NA, nrow = nOTU, ncol = nOTU)
  colnames(pcor) <- colnames(pcor.pval) <- taxa_names
  rownames(pcor) <- rownames(pcor.pval) <- taxa_names

  if (is.null(lib_size) & !is.transformed) {
    # Creating library size by summing taxa counts per sample and time point
    lib_size <- apply(data, c(1, 3), sum)
  }
  # Extract count array after excluding taxa
  data_filt <- data[, taxa_names, ]
  data_day_f <- data_filt[, , day]

  if (method == "pca") {
    loadings_m <- PCA_approx(data_day_f, ncomp = ncomp)
  } else if (method == "ica") {
    loadings_m <- RPCA_approx(data_day_f, ncomp = ncomp)
  } else if (method == "rpca") {
    loadings_m <- ICA_approx(data_day_f, ncomp = ncomp)
  } else {
    stop("Method not supported\n. Use one of pca, ica, rpca.")
  }

  for (i in 1:len) {
    # print(i)
    loading_tmp <- loadings_m
    loading_tmp[c(taxa1[i], taxa2[i]), ] <- 0
    u1 <- scale(data_day_f, center = TRUE, scale = TRUE) %*% loading_tmp
    if (!is.transformed) {
      ##** LUPINE_single with counts**##
      # principal component regression on log counts+1 with an offset for library size
      r_i <- lm(log(data_day_f[, taxa1[i]] + 1) ~ u1, offset = log(lib_size[, day]))
      r_j <- lm(log(data_day_f[, taxa2[i]] + 1) ~ u1, offset = log(lib_size[, day]))
    } else {
      ##** LUPINE_single with clr**##
      # principal component regression on clr
      r_i <- lm(data_day_f[, taxa1[i]] ~ u1)
      r_j <- lm(data_day_f[, taxa2[i]] ~ u1)
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

  pcor.full <- matrix(NA, nrow = nOTU_total, ncol = nOTU_total)
  pcor.pval.full <- matrix(NA, nrow = nOTU_total, ncol = nOTU_total)
  colnames(pcor.full) <- colnames(pcor.pval.full) <- colnames(data_day)
  rownames(pcor.full) <- rownames(pcor.pval.full) <- colnames(data_day)

  common_rows_cols <- intersect(rownames(pcor.full), rownames(pcor))

  pcor.full[common_rows_cols, common_rows_cols] <- pcor[common_rows_cols, common_rows_cols]
  pcor.pval.full[common_rows_cols, common_rows_cols] <- pcor.pval[common_rows_cols, common_rows_cols]

  res <- list(Estimate = pcor.full, pvalue = pcor.pval.full)

  return(res)
}
