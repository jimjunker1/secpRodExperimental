#'
#'
#'
#'

expMortSim <- function(t = NULL,n_init = NULL, z = NULL){
  if(is.null(n_init) | is.null(z)) stop("Error: 'n_init' and 'z' must be provided.")
  
  nVec = n_init*exp(-z*t)
  
  return(nVec)
}