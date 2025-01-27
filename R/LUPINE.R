#' LUPINE function
#'
#' @param data The data matrix (rows are samples, columns are features/variables, and slices are time points)
#' @param is.transformed A logical indicating whether the data is transformed or not
#' @param lib_size A matrix of library sizes for each sample and time point
#' @param ncomp The number of components to use for dimensionality reduction
#' @param single A logical indicating whether to use LUPINE for a single time point or longitudinal data
#' @param singleMethod The method to use for dimensionality reduction for single time point
#' @param cutoff The cutoff value for the p-value to determine significance of the correlation
#' @param excluded_var A list of taxa to be excluded at each time point
#'
#'
#'
#' @return A list of correlations and p-values for all days
#' @export
#'
LUPINE <- function(data, is.transformed = FALSE, lib_size = NULL, ncomp = 1,
                   single = FALSE, singleMethod = "pca", excluded_var = NULL, cutoff = 0.05) {

  # Checks
  if(length(dim(data)) != 3) {
    stop("Data should be a 3D array with dimensions samples x taxa x time points")
  }
  if(is.null(colnames(data))) {
  stop("Data should have variable names in columns")
  }

  # extract number of timepoints
  nTimepoints <- dim(data)[3]

  # call LUPINE_single or LUPINE_longitudinal functions
  if (single) {
    res <- sapply(1:nTimepoints, function(d) {
      net <- LUPINE_single(data,
                           timepoint = d,
                           excluded_var,
                           is.transformed,
                           lib_size,
                           method = singleMethod,
                           ncomp
      )$pvalue
      net <- apply(net<cutoff, c(1, 2), function(x) {
        ifelse(is.na(x), 0, x)
      })
      return(net)
    }, simplify = FALSE)
  } else {
    res <- sapply(2:nTimepoints, function(d) {
      net <- LUPINE_longitudinal(data,
        day = d, excluded_var, is.transformed,
        lib_size, ncomp
      )$pvalue
      net <- apply(net<cutoff, c(1, 2), function(x) {
        ifelse(is.na(x), 0, x)
      })
      return(net)
    }, simplify = FALSE)
  }

  return(res)
}
