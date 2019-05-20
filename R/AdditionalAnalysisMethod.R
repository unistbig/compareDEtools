#' List available *.createRmd functions
#'
#' Print a list of all \code{*.createRmd} functions that are available in the search path. These functions can be used together with the \code{\link{runDiffExp}} function to perform differential expression analysis. Consult the help pages for the respective functions for more information.
#' @export
#' @author Charlotte Soneson
#' @examples
#' listcreateRmd()
listcreateRmd <- function() {
  s <- unlist(sapply(search(), ls, all.names = TRUE))
  print(unname(s[grep("\\.createRmd$", s)]))
}









