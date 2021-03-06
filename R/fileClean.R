#' cleanR
#'
#' Function to clean cnv files
#' @param cnv_file File to clean [Required]
#' @param region More stringent filter for larger windows as defualt
#' @keywords clean
#' @import dplyr
#' @export


cleanR <- function(cnv_file, region=F) {
  options(scipen=1000000)
  clean_file <- dplyr::filter(cnv_file, is.finite(log2) )

  if(region==F){
    cat("Filtering on number of reads per window (each window must have more than 100 reads in tum & cont)", "\n")
    clean_file <- dplyr::filter(clean_file, test > 100)
  }
  clean_file <- droplevels(clean_file)
  return(clean_file)
}
