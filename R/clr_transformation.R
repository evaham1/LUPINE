#' Centered log ratio transform
#'
#' Compute the centered log ratio transform of a data matrix
#'
#' @param data The data matrix (rows are samples and columns are features/variables)
#'
#' @return The transformed data matrix
#' @export
clr_trans<-function(data){
  trans_data<-t(SpiecEasi::clr(t(data)))
  return(trans_data)
}
