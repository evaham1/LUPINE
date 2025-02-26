#' LUPINE for longitudinal data
#'
#' @param data A 3D array of counts or transformed data with dimensions samples x taxa x time points. The last timepoint in this data will be considered the 'current' one
#' for which the correlations are calculated
#' @param is.transformed A logical indicating whether the data is transformed or not
#' @param lib_size A matrix of library sizes for each sample and time point
#' @param ncomp The number of components to use for dimensionality reduction
#'
#' @return A list of correlations and p-values
#' @export
#'
LUPINE_multiple_timepoint <- function(data,
                                      is.transformed = FALSE,
                                      lib_size = NULL,
                                      ncomp = 1) {

  # Extract the data from the final timepoint
  data_final_timepoint <- data[, , dim(data)[3]] # equivalent to count_day <- data[, , day]

  # Extract variable names and calculate the number of variables
  var_names <- colnames(data_final_timepoint)
  nVar <- dim(data)[2]

  # pairwise taxa combinations
  nOTU <- length(var_names)
  len <- nOTU * (nOTU - 1) / 2
  taxa1 <- unlist(lapply(1:nOTU, function(i) rep(i, (nOTU - i))))
  taxa2 <- unlist(lapply(2:nOTU, function(i) seq(i, nOTU, 1)))

  # Initialisation
  pcor <- matrix(NA, nrow = nOTU, ncol = nOTU)
  pcor.pval <- matrix(NA, nrow = nOTU, ncol = nOTU)
  colnames(pcor) <- colnames(pcor.pval) <- var_names
  rownames(pcor) <- rownames(pcor.pval) <- var_names

  if (is.null(lib_size) & !is.transformed) {
    # Creating library size by summing taxa counts per sample and time point
    lib_size <- apply(data, c(1, 3), sum)
  }
  # Extract count array after excluding taxa for current day
  data_day_f <- data[, var_names, dim(data)[3]]

  loadings_m <- PLS_approx(data, dim(data)[3], ncomp = ncomp, taxa_names = var_names)


  for (i in 1:len) {
    # print(i)
    loading_tmp <- loadings_m
    loading_tmp[c(taxa1[i], taxa2[i]), ] <- 0
    u1 <- scale(data_day_f, center = TRUE, scale = TRUE) %*% loading_tmp
    if (!is.transformed) {
      ##** LUPINE with counts**##
      # principal component regression on log counts+1 with an offset for library size
      r_i <- lm(log(data_day_f[, taxa1[i]] + 1) ~ u1, offset = log(lib_size[, dim(data)[3]]))
      r_j <- lm(log(data_day_f[, taxa2[i]] + 1) ~ u1, offset = log(lib_size[, dim(data)[3]]))
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

  pcor.full <- matrix(NA, nrow = nVar, ncol = nVar)
  pcor.pval.full <- matrix(NA, nrow = nVar, ncol = nVar)
  colnames(pcor.full) <- colnames(pcor.pval.full) <- colnames(data)
  rownames(pcor.full) <- rownames(pcor.pval.full) <- colnames(data)

  common_rows_cols <- intersect(rownames(pcor.full), rownames(pcor))

  pcor.full[common_rows_cols, common_rows_cols] <- pcor[common_rows_cols, common_rows_cols]
  pcor.pval.full[common_rows_cols, common_rows_cols] <- pcor.pval[common_rows_cols, common_rows_cols]

  res <- list(Estimate = pcor.full, pvalue = pcor.pval.full)

  return(res)
}
