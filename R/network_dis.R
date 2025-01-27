#' Distance matrix
#'
#' @param network_ls network list
#'
#' @return data frame with IVI values
#' @export
#'
#' @importFrom NetworkDistance nd.gdd
distance_matrix<- function(network_ls){

  # calculate number of networks in network list
  num_net<-length(network_ls)
  # creates a new list where each element is a network
  g_ls<-lapply(1:num_net,function(k){network_ls[[k]]})
  # run graph diffusion distance calculation on networks
  dst<-nd.gdd(g_ls, out.dist = FALSE)$D
  return(dst)
}
