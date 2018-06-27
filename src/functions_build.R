check.replicates.number <- function(replicate.multiplexing, replicates) {
  number.of.replicates <- length(unique(replicates)) 
  if (replicate.multiplexing == FALSE &&
      number.of.replicates <= 1) {
    return (FALSE)
  } else {
    return (TRUE)
  }
}
