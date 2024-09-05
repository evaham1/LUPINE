#' LUPINE_bootsrap function
#'
#' @param data The data matrix (rows are samples, columns are features/variables, and slices are time points)
#' @param is.transformed A logical indicating whether the data is transformed or not
#' @param lib_size A matrix of library sizes for each sample and time point
#' @param ncomp The number of components to use for dimensionality reduction
#' @param single A logical indicating whether to use LUPINE for a single time point or longitudinal data
#' @param singleMethod The method to use for dimensionality reduction for single time point
#' @param day_range day_range The range of days to perform LUPINE bootstrap on
#' @param cutoff cutoff The cutoff value for the p-value to determine significance of the correlation
#' @param nboot The number of bootstrap iterations to perform
#' @param excluded_taxa A list of taxa to be excluded at each time point
#'
#' @return A list of correlations and p-values for all days
#' @export
#' @importFrom pbapply pblapply
#' @importFrom abind abind
LUPINE_bootsrap <- function(data, day_range= NULL, is.transformed = FALSE, lib_size = NULL, ncomp = 1,
                            single = FALSE, singleMethod = "pca", excluded_taxa = NULL, cutoff = 0.05, nboot = 1000) {
  nDays <- dim(data)[3]

  bootstrap_data<- lapply(1:nboot,function(i){
    boot_data<- data[sample(1:dim(data)[1], replace = T),,]
    rownames(boot_data)<-rownames(data)
    return(boot_data)
  })

  # Initialize an empty list to store median ,lower, and upper pvalues
  boot_res <-list()

  if(is.null(day_range)){
    day_range<-1:nDays
  }

  for(i in length(day_range)){
    if(single){
      res <- pblapply(
        1:nboot,
        #Library size is set to NULL so that it is calculated for each bootstrap iteration
        completed_counter(function(iter){net<-LUPINE_single(bootstrap_data[[iter]], day = day_range[i],
                                          excluded_taxa, is.transformed, lib_size=NULL,
                                          method = singleMethod, ncomp)$pvalue
        net <- apply(net, c(1, 2), function(x) {ifelse(is.na(x),0,x)})
        return(net)
        }, nboot))
      pvalues<-res %>%
        abind(., along=3)
      median_mt<-apply(pvalues, c(1,2), median)
      lower_mt<-apply(pvalues, c(1,2), function(x) quantile(x,probs=0.025))
      upper_mt<-apply(pvalues, c(1,2), function(x) quantile(x,probs=0.975))
      boot_res[[i]]<-list(median_mt=median_mt, lower_mt=lower_mt, upper_mt=upper_mt)
      names(boot_res)[i] <- paste0("Day_",  day_range[i])
    }else if(!single & day_range[i]>1){
      res <- pblapply(
        1:nboot,
        #Library size is set to NULL so that it is calculated for each bootstrap iteration
        completed_counter(function(iter){net<-LUPINE_longitudinal(bootstrap_data[[iter]], day = day_range[i],
                                                excluded_taxa, is.transformed, lib_size=NULL, ncomp)$pvalue
        net <- apply(net, c(1, 2), function(x) {ifelse(is.na(x),0,x)})
        return(net)
        }, nboot))
      pvalues<-res %>%
        abind(., along=3)
      median_mt<-apply(pvalues, c(1,2), median)
      lower_mt<-apply(pvalues, c(1,2), function(x) quantile(x,probs=0.025))
      upper_mt<-apply(pvalues, c(1,2), function(x) quantile(x,probs=0.975))
      boot_res[[i]]<-list(median_mt=median_mt, lower_mt=lower_mt, upper_mt=upper_mt)
      names(boot_res)[i] <- paste0("Day_",  day_range[i])
    }
  }

  return(boot_res)
}


