check.replicates.number <- function(replicate.multiplexing, replicates) {
  # Check if the replicates are more than 1, providing that we do not utilize replicates multiplexing
  #
  # Args:
  #   replicate.multiplexing: Boolean variable regarding the existence or not of replicates multiplexing
  #   replicates:             A vector of replicates 
  # Returns:
  #   TRUE if the number of replicates is correct, otherwise FALSE
  number.of.replicates <- length(unique(replicates)) 
  if (replicate.multiplexing == FALSE &&
      number.of.replicates <= 1) {
    return (FALSE)
  } else {
    return (TRUE)
  }
}

replicates.per.condition <- function(biological.replicates, technical.replicates, conditions) {
  # Finds out which biological and technical replicates belong to each condition
  #
  # Args:
  #   biological.replicates:A vector with the biological replicates
  #   technical.replicates: A vector with the technical replicates
  #   conditions:           A vector with the conditions
  # Returns:
  #   A list of lists. Each element contains the 2 biological and technical replicates 
  #   vectors, belonging to each condition
  #
  #   condition: The condition name 
  #     biological: A vector with the biological replicates for the parent condition
  #     technical:  A vector with the technical replicates for the parent condition
  
  # Make empty list for replicates per condition
  replicates.per.condition <- list()
  
  # Find the unique conditions
  unique.conditions <- unique(conditions)
  
  for (condition in unique.conditions) {
    
    # How many replicates do I have for this condition?
    raw.files.per.condition <- which(condition == conditions)
    
    # Store the biological and technical indexes of these replicates
    biological.replicates.per.condition <- biological.replicates[raw.files.per.condition]
    technical.replicates.per.condition <- technical.replicates[raw.files.per.condition]
    
    # Wrap them in a list element
    new.element <- list("biological" = biological.replicates.per.condition,
                        "technical"= technical.replicates.per.condition)
    
    # Add put the in the list
    if (length(replicates.per.condition) == 0 ) {
      replicates.per.condition <- new.element
    } else {
      replicates.per.condition <- list(replicates.per.condition, new.element)
    }
  }
  
  # Finally put the conditions names for each element
  names(replicates.per.condition) <-  unique.conditions
  
  return (replicates.per.condition)
}

check.replicates <- function(biological.replicates, technical.replicates) {
  # Checks if the replicates numbering is correct e.g. 1,2,3 is correct and 1,2,4 is wrong
  #
  # Args:
  #   biological.replicates:A vector of biological replicates
  #   technical.replicates: A vector of technical replicates
  #
  # Returns:
  #   A list with 2 booleans for each replicate type (biological/technical) regarding 
  #   if the numbering is correct or not
  #   biological: TRUE if the numbering is correct, otherwise FALSE
  #   technical: TRUE if the numbering is correct, otherwise FALSE
  
  # Find the unique numbers for each vector and sort them
  unique.biological.replicates <-  unique(biological.replicates)
  unique.technical.replicates <-  sort(unique(technical.replicates))
  
  expected.biological.replicates <- c(1:length(unique.biological.replicates))
  expected.technical.replicates <- c(1:length(unique.technical.replicates))
  
  # If the length of each vector equals the unique numbers,
  # then the replicate vector is correct, otherwise it is incorrect 
  biological.replicates.are.ok <- all(unique.biological.replicates == expected.biological.replicates)
  technical.replicates.are.ok <- all(unique.technical.replicates == expected.technical.replicates)
  
  # Now put them in a list
  replicates.status <- list( "biological" = biological.replicates.are.ok,
                             "technical" = technical.replicates.are.ok)
  
  return (replicates.status)
}

replicates.status.per.condition <- function(replicates.per.condition) {
  # Makes a list of list, in which each element list corresponds to a condition
  # and contains 2 booleans for the correctness of the biological and technical
  # replicates numbering for this condition
  #
  # Args:
  #   replicates.per.condition: A list of lists in which each element correspond to a condition
  #                             and the biological and technical replicates that are assossiated
  #                             with it.
  # Returns:
  #   A list of lists. Each element contains the 2 booleans regarding the correctness of the
  #   biological and technical replicates belonging to each condition
  #
  #   condition: The condition name 
  #     biological: A boolean for the correctness of the biological replicates
  #     technical:  A boolean for the correctness of the technical replicates
  
  # Make the list
  replicates.status.per.condition <- list()
  
  for (condition in replicates.per.condition) {
    biological.replicates <- condition[[1]]
    technical.replicates <- condition[[2]]
    
    # The the correctness of the biological/technical replicates of the 
    # current status
    new.status <- check.replicates(biological.replicates,
                                   technical.replicates)
    
    # Add them to a list
    if (length(replicates.status.per.condition) == 0) {
      replicates.status.per.condition <- new.status
    } else {
      replicates.status.per.condition <- list(replicates.status.per.condition, new.status)
    }
  }
  # Finally add the condition name
  names(replicates.status.per.condition) <- names(replicates.per.condition)
  return (replicates.status.per.condition)
}