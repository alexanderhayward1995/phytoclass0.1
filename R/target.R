#' Title
#'
#' @param x   xx
#'
#' @return
#' @export
#'
#' @examples
target = function(x){
  return(ifelse(x<0,0,exp(-x)))
}
