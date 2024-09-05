#' Distance matrix
#'
#' @param network_ls network list
#'
#' @return data frame with IVI values
#' @export
#'
#' @importFrom NetworkDistance nd.gdd
distance_matrix<- function(network_ls){

  num_net<-length(network_ls)

  g_ls<-lapply(1:num_net,function(k){network_ls[[k]]})

  dst<-nd.gdd(g_ls, out.dist = FALSE)$D
  return(dst)
}
