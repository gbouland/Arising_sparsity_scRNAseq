## Binarizes single-cell count matrix.
## Also works for log-normalzied matrix
binarize <- function(counts){
  binary <- (counts>0)*1
  return(binary)
}




