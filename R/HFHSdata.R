#' @title Example Dataset: HFHS data
#'
#' @description The example data is a list containing 8 elements:
#' \enumerate{
#'   \item OTU data for HFHS diet: an array with dimensions n=23, p=212, t=4.
#'   \item OTU data for Normal diet: an array with dimensions n=23, p=212, t=4.
#'   \item Library size for HFHS diet: a matrix with dimensions n=23, t=4.
#'   \item Library size for Normal diet: a matrix with dimensions n=23, t=4.
#'   \item Filtered sample information: a data frame with dimensions.
#'   \item Filtered taxonomy data: a data frame with dimension.
#'   \item Low abundant taxa for HFHS diet: a list of 4 elements.
#'   \item Low abundant taxa for Normal diet: a list of 4 elements.
#' }
#'
#' @format A list containing 8 elements, including arrays, matrices, data frames, and lists, as described above.

#' @examples
#' data(HFHSdata)
#'
#' @references Kodikara, S., & Le Cao, K. A. (2024). Microbial network inference for longitudinal microbiome studies with LUPINE. bioRxiv, 2024-05.
#'
#' @keywords datasets
"HFHSdata"
