#' LUPINE for longitudinal data
#'
#' @param data A 3D array of counts or transformed data with dimensions samples x taxa x time points
#' @param day_index The index of the time point for which the correlations are calculated
#' @param num_lags The number of previous time points to be used for the calculation of the correlation
#' @param excluded_taxa A list of taxa to be excluded at each time point
#' @param is.transformed A logical indicating whether the data is transformed or not
#' @param lib_size A matrix of library sizes for each sample and time point
#' @param ncomp The number of components to use for dimensionality reduction
#'
#' @return A list of correlations and p-values
#' @export
#'
LUPINE_t <- function(data, day_index, num_lags=999, excluded_taxa = NULL, is.transformed = FALSE, lib_size = NULL, ncomp = 1) {

  if (num_lags==0) {
    res <- LUPINE_single(data, day_index, excluded_taxa, is.transformed, lib_size,
                         method = "pca", ncomp)
  } else if(num_lags==999) {
    res <- LUPINE_longitudinal(data, day_index, excluded_taxa, is.transformed, lib_size,
                               ncomp)
  } else {
    data_new<- data[,,(day_index-num_lags):day_index]
    day <- num_lags+1
    # Total number of taxa (p)
    nOTU_total <- dim(data_new)[2]
    # Extract matrix for current time point
    count_day <- data_new[, , day]
    # Extract taxa names after excluded taxa
    taxa_names <- colnames(count_day)[!(colnames(count_day) %in% excluded_taxa[[day_index]])]
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
      lib_size <- apply(data_new, c(1, 3), sum)
    }
    # Extract count array after excluding taxa for current day
    data_day_f <- data_new[, taxa_names, day]

    loadings_m <- PLS_approx(data_new, day, ncomp = ncomp, taxa_names = taxa_names)


    for (i in 1:len) {
      # print(i)
      loading_tmp <- loadings_m
      loading_tmp[c(taxa1[i], taxa2[i]), ] <- 0
      u1 <- scale(data_day_f, center = TRUE, scale = TRUE) %*% loading_tmp
      if (!is.transformed) {
        ##** LUPINE with counts**##
        # principal component regression on log counts+1 with an offset for library size
        r_i <- lm(log(data_day_f[, taxa1[i]] + 1) ~ u1, offset = log(lib_size[, day_index]))
        r_j <- lm(log(data_day_f[, taxa2[i]] + 1) ~ u1, offset = log(lib_size[, day_index]))
      } else {
        ##** LUPINE with clr**##
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
    colnames(pcor.full) <- colnames(pcor.pval.full) <- colnames(data)
    rownames(pcor.full) <- rownames(pcor.pval.full) <- colnames(data)

    common_rows_cols <- intersect(rownames(pcor.full), rownames(pcor))

    pcor.full[common_rows_cols, common_rows_cols] <- pcor[common_rows_cols, common_rows_cols]
    pcor.pval.full[common_rows_cols, common_rows_cols] <- pcor.pval[common_rows_cols, common_rows_cols]

    res <- list(Estimate = pcor.full, pvalue = pcor.pval.full)

  }

  return(res)
}
