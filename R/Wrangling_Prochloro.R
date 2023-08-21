#' Converts data-types and selects data for randomisation in 
#' the simulated annealing algorithm with prochlorococcus
#'
#' @param Fl   xx
#' @param min.val xx
#' @param max.val  xx
#'
#' @return
#' @export
#'
#' @examples
Wrangling <- function(Fl, min.val, max.val){
  Fd <- Fl
  Fmin <- as.matrix(Fd)   #### Set up Fmin matrix
  Fmin <- Fmin[,1:(ncol(Fmin)-2)]
  Fmin[Fmin>0] <- min.val  # set all non-zero elements to the minimum values (imported from csv )
  Fmax <- as.matrix(Fd)
  Fmax <- Fmax[,1:(ncol(Fmax)-2)]
  Fmax[Fmax>0] <- max.val

  chlv <- Fd[,ncol(Fd)]  ##### The chlorophyll values once weighted to rowsums for initial F matrix
  pChl <- Fd[nrow(Fd),ncol(Fd)-1]   ##### The chlorophyll prochloro... 


  Fmin <- Fmin * chlv #### multiply the minimum value by weighted chlorophyll to updated ratios
  Fmin[nrow(Fmin),1:ncol(Fmin)-1] <-Fmin[nrow(Fmin),1:ncol(Fmin)-1] *pChl
  Fmin[1:nrow(Fmin)-1,ncol(Fmin)] <- chlv
  Fmin <- cbind(Fmin,chlv)  #### Reassign correct initial chl values
  Fmax[1:nrow(Fmin)-1,1:ncol(Fmax)] <- Fmax[1:nrow(Fmax)-1,1:ncol(Fmax)] * chlv  #### same for max values
  
  Fmax[nrow(Fmin),ncol(Fmax)-1] <-Fmax[nrow(Fmax),ncol(Fmax)-1] *pChl
  Fmax[1:nrow(Fmin)-1,ncol(Fmax)] <- chlv
  Fmin[nrow(Fmin),ncol(Fmin)-1] <- pChl
  Fmax[nrow(Fmax),ncol(Fmax)-1] <- pChl

  
  Fmin <- vectorise(Fmin[,1:ncol(Fmin)-1]) #### vectorise function outputs all non-zero elements as a vector (excluding chl column)
  Fmax <- vectorise(Fmax[,1:ncol(Fmax)-1]) #### vectorise function outputs all non-zero elements as a vector (excluding chl column)
  SE <- vectorise(Fd[,1:ncol(Fd)-1]) #### vectorise function outputs all non-zero elements as a vector (excluding chl column)
  Pro <- 0
  chlv <- c(chlv,Pro)

  res <- list(Fmin, Fmax, SE,chlv)
}
