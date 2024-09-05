#' PLS low dimensional approximaton for longitudinal data
#'
#' @param X_array The data array
#' @param day The day to approximate
#' @param ncomp The number of components to use for PLS
#' @param taxa_names The taxa names to use for PLS
#'
#' @return The low dimensional data matrix using PLS
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom mixOmics pls
#' @importFrom mixOmics block.pls
#' @importFrom rlist list.append
PLS_approx <- function(X_array, day, ncomp = 1, taxa_names = NULL) {
  # Extract array after excluding taxa
  X_f <- X_array[, taxa_names, ]
  # Extract day
  X_f.day <- X_f[, , day]

  ## checking which taxa has zero variance to avoid errors in pls##
  zero_var_index <- which(X_f.day %>% apply(., 2, var) == 0)

  ## Creating d1 to include all past data while making sure no zero variance taxa##
  day_i <- day
  zero_var_tmp <- which(X_f[, , (day_i - 1)] %>% apply(., 2, var) == 0)
  d1 <- list(X1 = if (length(zero_var_tmp) != 0) {
    X_f[, -zero_var_tmp, (day_i - 1)]
  } else {
    X_f[, , (day_i - 1)]
  })
  day_i <- day_i - 1
  index <- 3
  while (day_i != 1) {
    zero_var_tmp <- which(X_f[, , (day_i - 1)] %>% apply(., 2, var) == 0)
    d_tmp <- if (length(zero_var_tmp) != 0) {
      X_f[, -zero_var_tmp, (day_i - 1)]
    } else {
      X_f[, , (day_i - 1)]
    }
    d1 <- list.append(d1, d_tmp)
    names(d1)[index - 1] <- paste("X", index - 1, sep = "")
    index <- index + 1
    day_i <- day_i - 1
  }

  if (length(zero_var_index) > 0) {
    if (day == 2) {
      res_netPLS1 <- pls(d1$X1, X_f.day[, -zero_var_index],
        ncomp = ncomp, max.iter = 1000
      )
    } else {
      res_netPLS1 <- block.pls(d1, X_f.day[, -zero_var_index],
        ncomp = ncomp, max.iter = 1000
      )
    }

    # loadings matrix
    loadings_m <- res_netPLS1$loadings$Y

    # Including a loading of zero to taxa with zero variance
    for (i in 1:(length(zero_var_index))) {
      if (zero_var_index[i] == 1) {
        loadings_m <- rbind(rep(0, ncomp), loadings_m)
      } else {
        loadings_m <- rbind(
          loadings_m[1:(zero_var_index[i] - 1), ],
          rep(0, ncomp),
          loadings_m[zero_var_index[i]:dim(loadings_m)[1], ]
        )
      }
    }
  } else {
    if (day == 2) {
      res_netPLS1 <- pls(d1$X1, X_f.day,
        ncomp = ncomp, max.iter = 1000
      )
    } else {
      res_netPLS1 <- block.pls(d1, X_f.day,
        ncomp = ncomp, max.iter = 1000
      )
    }
    # loading matrix
    loadings_m <- res_netPLS1$loadings$Y
  }

  return(loadings_m)
}
