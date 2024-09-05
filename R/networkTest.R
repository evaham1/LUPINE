#' Network Mantel test
#'
#' @param network_ls network list
#'
#' @return matrix of p-values
#' @export
#'
#' @importFrom ade4 mantel.rtest
#' @importFrom e1071 hamming.distance
MantelTest_matrix<- function(network_ls){
  num_net<-length(network_ls)

  mantel.pvalue_matrix <- matrix(0, num_net, num_net)

  # Calculate pairwise mantel.pvalue
  for (i in 1:(num_net-1)) {
    for (j in (i+1):num_net) {

      lapl1<-as.dist(hamming.distance(network_ls[[i]]))

      lapl2<-as.dist(hamming.distance(network_ls[[j]]))

      mantel.pvalue_matrix[i,j]<- mantel.pvalue_matrix[j,i]<-mantel.rtest(lapl1,lapl2)$pvalue
    }
  }
  return(mantel.pvalue_matrix)
}
