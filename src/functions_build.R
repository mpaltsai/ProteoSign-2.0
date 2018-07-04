check.replicates.number <- function(replicate.multiplexing, replicates) {
  #
  # Check if the replicates are more than 1, providing that we do not utilize replicates multiplexing
  #
  # Args:
  #   replicate.multiplexing: Boolean variable regarding the existence or not of replicates multiplexing
  #   replicates:             A vector of replicates 
  #
  # Returns:
  #   TRUE:   The number of replicates is correct
  #   FALSE:  The number of replicates is incorrect
  #
  
  # Get number of unique replicates 
  number.of.replicates <- length(unique(replicates))
  
  # Assure that number of replicates is more than one
  # (We only allow this if there is replicate multiplexing)
  if (replicate.multiplexing == FALSE &&
      number.of.replicates <= 1) {
    return (FALSE)
  } else {
    return (TRUE)
  }
}

replicates.per.condition <- function(biological.replicates, technical.replicates, conditions) {
  #
  # Finds out which biological and technical replicates belong to each condition
  #
  # Args:
  #   biological.replicates:A vector with the biological replicates
  #   technical.replicates: A vector with the technical replicates
  #   conditions:           A vector with the conditions
  #
  # Returns:
  #   A list of lists. Each element contains the 2 biological and technical replicates 
  #   vectors, belonging to each condition
  #
  #   condition: The condition name 
  #     biological: A vector with the biological replicates for the parent condition
  #     technical:  A vector with the technical replicates for the parent condition
  #
  
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
  #
  # Checks if the replicates numbering is correct e.g. 1,2,3 is correct and 1,2,4 is wrong
  #
  # Args:
  #   biological.replicates:A vector of biological replicates
  #   technical.replicates: A vector of technical replicates
  #
  # Returns:
  #   A list with 2 booleans for each replicate type (biological/technical) regarding 
  #   if the numbering is correct or not
  #
  #   biological: TRUE if the numbering is correct, otherwise FALSE
  #   technical:  TRUE if the numbering is correct, otherwise FALSE
  #
  
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
  #
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
  #
  
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

fix.replicates <- function(biological.replicates, technical.replicates) {
  #
  # Fix the biological and technical replicates of a particular condition
  #
  # Args:
  #   biological.replicates:A vector with the biological replicates of a condition
  #   technical.replicates: A vector with the technical replicates of a condition
  #
  # Returns:
  #   A list with the corrected biological/technical replicates vectors
  #
  
  # Get the unique biological replicates
  unique.biological <- unique(biological.replicates)
  
  # Initialize the starting bioligical sample id
  biological.sample.id <- 1
  
  # Traverse over the unique biological replicates
  for (biorep in unique.biological) {
    
    # Find the same biological replicates
    same.samples.biological <- which(biological.replicates == biorep)
    
    # Correct the numbering of them
    biological.replicates[same.samples.biological] <- biological.sample.id
    
    # Increase the numbering by 1
    biological.sample.id <- biological.sample.id + 1
    
    # Now get the unique technical replicates for a 
    # particular biological replicate
    unique.technical <- unique(technical.replicates[same.samples.biological])
    
    # Initialize the starting technical sample id
    technical.sample.id <- 1
    
    # Traverse though the unique technical replicates for a particular
    # biological replicate
    for (techrep in unique.technical) {
      
      # Find the technical replicates sharing the same technical replicate numbering
      same.samples.technical <- which(technical.replicates == techrep &
                                        biological.replicates == biological.sample.id - 1)
      
      # Correct their numbering
      technical.replicates[same.samples.technical] <- technical.sample.id
      
      # Increase the numbering by 1
      technical.sample.id <- technical.sample.id +1
    }
  }
  
  # Finally make the list of the corrected vectors
  fixed.replicates <- list(biological.replicates, technical.replicates)
  names(fixed.replicates) <- c("biological", "technical")
  return (fixed.replicates)
}

fix.replicates.per.condition <- function (replicates.per.condition, replicates.status.per.condition) {
  #
  # Fixes the bad replicates of replicates.per.condition
  #
  # Args:
  #   replicates.per.condition: A list of lists where each element corresponds
  #   corresponds to a conditions and it contains 2 vectors regarding the
  #   biological/technical replicates IDs
  #
  #   replicates.status.per.condition: A list of lists where each element 
  #   corresponds to a conditions and it contains 2 booleans regarding 
  #   the correctness of the biological/technical replicates for this
  #   conditions.
  #
  # Returns:
  #   replicates.per.condition: The corrected replicates.per.condition
  #   vector with the correct numbering of biological/technical
  #   replicates
  #
  
  # Traverse over replicates per condition
  for (index in 1:length(replicates.status.per.condition)) {
    condition <- names(replicates.status.per.condition)[index]
    biological.status <- replicates.status.per.condition[[index]]$biological
    technical.status  <- replicates.status.per.condition[[index]]$technical
    
    # Make a code using the biological/technical replicates status
    # for faster comparison than if/else
    case <- (biological.status * 2  + technical.status * 1) + 1
    
    switch(case,
           cat("Condition", condition, "has bad biological and technical replicates.\n") ,
           cat("Condition", condition, "has bad biological replicates.\n") ,
           cat("Condition", condition, "has bad technical replicates.\n") ,
           cat("Condition", condition, "has good replicates.\n"))
    
    # If replicates are not good
    if (case != 4) {
      cat("Fixing condition ", condition," ...\n")
      
      # Get the biological/technical replicates vector for a particural
      # to be fixed
      biological.replicates <- replicates.per.condition[[index]]$biological
      technical.replicates <- replicates.per.condition[[index]]$technical
      
      # Now fix the bad replicates
      fixed.replicates <- fix.replicates(biological.replicates,
                                         technical.replicates)
      
      # Finally repair the replicates pern condition list with the 
      # corrected vectors
      replicates.per.condition[[index]] <- fixed.replicates
    }
  }
  return (replicates.per.condition)
}

restore.replicates <- function(replicates.per.condition) {
  #
  # Creates a list containing 2 vectors in roder to restore the initial
  # experimental dataset columns with the correct biological/technical
  # replicate vectors.
  #
  # Args:
  #   replicates.per.condition: A list of lists where each list element
  #   corresponds to a condition and its biological/technical replicates
  #   numbering.
  #
  # Returns:
  #   A list with 2 vectors containing the corrected biological/technical
  #   replicates vector
  #
  
  # Create the new vectors
  biological <- c()
  technical <- c()
  
  # Iterate over the conditions
  for (rep in replicates.per.condition) {
    
    # Add the biological/technical replicates of each condition to the vector
    biological <- c(biological, rep$biological)
    technical <- c(technical, rep$technical)
  }
  # Finally wrap them in a list
  restored.replicates <- list("biological" = biological, "technical" = technical)
  return (restored.replicates)
}

make.experimental.description <- function(experimental.setup.id, biological.replicates, technical.replicates, fractions) {
  #
  # Creates a vector with the experimental description of each sample
  # as the concatenation of BXTYFZ where XYZ are IDs for the biological/
  # technical/fraction
  # 
  # Args:
  #   experimental.setup.id:  An id indicating the experimental setup type
  #   biological.replicates:  A vector with the biological replicates 
  #   technical.replicates:   A vector with the technical replicates
  #   fractions:              A vector with the fractions
  #
  # Returns:
  #   A vector with the experimental description
  #
  
  # Make the empty description vector
  experimental.description <- c()
  
  # Paste B,T,F to each biological/technical/fractions vector
  experimental.description.biological <-  paste0("B", biological.replicates)
  experimental.description.technical <-   paste0("T", technical.replicates)
  experimental.description.fractions <-   paste0("F", fractions)
  
  # And pasted the appropriate columns based on the experimental id
  switch(experimental.setup.id,
         {
           cat("We have bioreps.\n")
           experimental.description <- experimental.description.biological
         },
         {
           cat("No bioreps, no fractions? Really? \n")
           cat("Invalid experimental description id in make.experimental.description function.\n")  
           return (FALSE)
         },
         {
           cat("We have bioreps and techreps.\n")
           experimental.description <- paste0(experimental.description.biological,
                                              experimental.description.technical)
         },
         {
           cat("We have bioreps and fractions.\n")
           experimental.description <- paste0(experimental.description.biological,
                                              experimental.description.fractions)
         },
         {
           cat("No bioreps? Really?\n")
           cat("Invalid experimental description id in make.experimental.description function.\n")  
           return (FALSE)
         },
         {
           cat("We have bioreps, techreps and fractions.\n")
           experimental.description <- paste0(experimental.description.biological,
                                              experimental.description.technical,
                                              experimental.description.fractions)
         })
  cat(experimental.description,"\n")
  return (experimental.description)
}

