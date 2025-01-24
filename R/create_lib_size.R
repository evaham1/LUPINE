#' Create library size for a dataset at one timepoint or multiple timepoints
#'
#' @param data A 3D or 2D array of counts or transformed data with dimensions samples x taxa
#' @return Library size across samples and optionally timepoints
#' @export
#'
create_lib_size <- function(data){

  if (length(dim(data)) == 3){
    print("Data has multiple timepoints, calculating library size across each timepoint...")
    lib_size <- apply(data, c(1, 3), sum)
  }

  else if (length(dim(data)) == 2){
    print("Data has only one timepoint, calculating library size...")
    lib_size <- apply(data, 1, sum)
  }

  else {stop("Data has must have either 2 or 3 dimensions!") }

  return(lib_size)

}

