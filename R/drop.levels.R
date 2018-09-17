######functions and libraries
################################################################################
########libraries and functions
drop.levels <- function(dat){
  # Drop unused factor levels from all factors in a data.frame
  # Author: Kevin Wright.  Idea by Brian Ripley.
  dat[] <- lapply(dat, function(x) x[,drop=TRUE])
  return(dat)
}
################################################################################