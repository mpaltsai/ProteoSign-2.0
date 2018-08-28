assign.condition.to.raw.files <- function(condition, raw.files, conditions.to.raw.files, is.label.free) {
  #
  # Add a new condition along side with their raw files or terminates the analysis
  #
  # Args:
  #   condition:                The string of the new condition to add.
  #   raw.files:                A vector of strings with the raw file names
  #   conditions.to.raw.files:  The conditions.to.raw.files list
  #   is.label.free:            The boolean to say whenever we are on a labeled or not experiment
  #
  # Returns:
  #   The updated condition.to.raw.files list or in case of error terminates the analysis
  #
  
  # The new element to add
  condition.to.raw.files <- list()
  
  # Initialize experimental definition. We have conditions for label-free experiments and labels for labeled experiments
  experimental.definition <- ""
  
  # Set the experimental definition depending on the experiment
  if (is.label.free == TRUE) {
    experimental.definition <- "condition"
  }
  else {
    experimental.definition <- "label"
  }
  
  # Does the condition we are trying to add exists?
  conditions.exists <- condition %in% names(conditions.to.raw.files)
  
  # Lets assume that we add the raw files
  future.list <- c(unlist(conditions.to.raw.files, use.names = FALSE), raw.files)
  
  # Will we end up with same raw files assigned to multiple conditions?
  conditions.have.duplicates <- length(which(duplicated(future.list) == TRUE)) != 0
  
  # Now set a status code
  # 1:  Everything is ok
  # 2:  The condition already exists
  # 3:  There are raw files assigned to multiple conditions
  # 4:  The condition already exists and the raw files assigned to multiple conditions
  status.code <-    conditions.exists * 1 + 
                    conditions.have.duplicates * 2 + 1
  
  # Print the appropriate error message and terminate the execution
  # or add the new condition with its raw files
  switch ( status.code,
           {
             # 1:  Everything is ok
             
             # Set the new condition and connect it with the raw files 
             condition.to.raw.files[[condition]] <- raw.files
           },
           {
             # 2:  The condition already exists
             
             # Set the message
             message <- paste0("assign.condition.to.raw.files: Error adding ",
                               experimental.definition,
                               " '",
                               condition,
                               "': An existing ",
                               experimental.definition,
                               " with name '",
                               condition,
                               "' already exists. Please try a different name.\n")
             
             # Add stop the execution of the analysis
             stop(message)
           },
           {
             # 3:  There are raw files assigned to multiple conditions
             
             # Set the message
             message <- paste0("assign.condition.to.raw.files: Error adding ",
                               experimental.definition,
                               " '",
                               condition,
                               "': The ",
                               experimental.definition,
                               " with name '",
                               condition,
                               "' share the same raw files ",
                               paste(raw.files, collapse=" "),
                               " with another ",
                               experimental.definition,
                               ". Please recheck your data.\n")
             
             # Add stop the execution of the analysis
             stop(message)
           },
           {
             # 4:  The condition already exists and the raw files assigned to multiple conditions
             
             # Set the message
             message <- paste0("assign.condition.to.raw.files: Error adding ",
                               experimental.definition,
                               " '",
                               condition,
                               "': An existing ",
                               experimental.definition,
                               " with name '",
                               condition,
                               "' already exists and share the same raw files ",
                               paste(raw.files, collapse=" "),
                               " with another ",
                               experimental.definition,
                               ". Please recheck your data and try a different condition name.\n")
             
             # Add stop the execution of the analysis
             stop(message)
           })
  
  return (condition.to.raw.files)
}

build.condition.to.raw.files.from.matrix <- function(raw.files.to.coditions.matrix, is.label.free) {
  #
  # Builds a conditions.to.raw.files.list from the experimental structure matrix
  #
  # Args:
  #   raw.files.to.coditions.matrix:  The subset matrix of experimental structrure table, 
  #                                   containing only the raw files and condition columns 
  #   is.label.free:                  Boolean, is the analysis label-free or not?
  #
  # Returns:
  #   A new conditions.to.raw.files.list containing the conditions and the raw files that 
  #   belong to each condition
  #

  # Initialize the new conditions.to.raw.files.list
  conditions.to.raw.files.list <- list()

  # Get all the conditions
  conditions <- raw.files.to.coditions.matrix$condition

  # Get all the raw.files
  raw.files <-  raw.files.to.coditions.matrix$raw.file

  # Get the unique conditions
  unique.conditions <- unique(conditions)
  # Now for each condition get the raw files that belong to itself
  for (condition in unique.conditions) {
    raw.files.group.indexes <- which(raw.files.to.coditions.matrix$condition == condition)
    raw.files.group <- raw.files[raw.files.group.indexes]

    # Make the new element
    conditions.to.raw.files.element <- assign.condition.to.raw.files(condition, 
                                                                  raw.files.group,
                                                                  conditions.to.raw.files.list,
                                                                  is.label.free)

    conditions.to.raw.files.list <- c(conditions.to.raw.files.list, conditions.to.raw.files.element)
  }
  return (conditions.to.raw.files.list)
}

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
  #     biological.replicates: A vector with the biological replicates for the parent condition
  #     technical.replicates:  A vector with the technical replicates for the parent condition
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
    new.element <- list("biological.replicates" = biological.replicates.per.condition,
                        "technical.replicates"= technical.replicates.per.condition)
    
    # Add element to list
    replicates.per.condition[[condition]] <- new.element
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
  #   biological.replicates: TRUE if the numbering is correct, otherwise FALSE
  #   technical.replicates:  TRUE if the numbering is correct, otherwise FALSE
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
  replicates.status <- list( "biological.replicates" = biological.replicates.are.ok,
                             "technical.replicates" = technical.replicates.are.ok)
  
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
  #     biological.replicates: A boolean for the correctness of the biological replicates
  #     technical.replicates:  A boolean for the correctness of the technical replicates
  #
  
  # Make the list
  replicates.status.per.condition <- list()
  
  for (index in 1:length(replicates.per.condition)) {
    biological.replicates <- replicates.per.condition[[index]][[1]]
    technical.replicates <- replicates.per.condition[[index]][[2]]
    
    # The the correctness of the biological/technical replicates of the 
    # current status
    new.status <- check.replicates(biological.replicates,
                                   technical.replicates)
    
    # Add them to a list
    replicates.status.per.condition[[index]] <- new.status
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
  names(fixed.replicates) <- c("biological.replicates", "technical.replicates")
  return (fixed.replicates)
}

fix.replicates.per.condition <- function (replicates.per.condition, replicates.status.per.condition, is.label.free = TRUE) {
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
  #   is.label.free:  Default is TRUE. The experiment is label-free
  #
  # Returns:
  #   replicates.per.condition: The corrected replicates.per.condition
  #   vector with the correct numbering of biological/technical
  #   replicates
  #
  
  # Traverse over replicates per condition
  for (index in 1:length(replicates.status.per.condition)) {
    condition <- names(replicates.status.per.condition)[index]
    biological.status <- replicates.status.per.condition[[index]]$biological.replicates
    technical.status  <- replicates.status.per.condition[[index]]$technical.replicates
    
    # Make a code using the biological/technical replicates status
    # for faster comparison than if/else
    case <- (biological.status * 2  + technical.status * 1) + 1 
    
    # If the experiment is labeled change the print message for the cases
    if (is.label.free == FALSE) {
      case <- case + 4
    }
    
    switch(case,
         cat("Condition", condition, "for the label-free experiment has bad biological and technical replicates.\n") ,
         cat("Condition", condition, "for the label-free experiment has bad biological replicates.\n") ,
         cat("Condition", condition, "for the label-free experiment has bad technical replicates.\n") ,
         cat("Condition", condition, "for the label-free experiment has good replicates.\n"),
         cat("Experimental setup for the labeled experiment has bad biological and technical replicates.\n") ,
         cat("Experimental setup for the labeled experiment has bad biological replicates.\n") ,
         cat("Experimental setup for the labeled experiment has bad technical replicates.\n") ,
         cat("Experimental setup for the labeled experiment has good biological and technical replicates.\n"))
    
    # If the replicates are not good for both experimental setups
    if (case != 4 & case != 8) {
      
      # Print the appropriate message depending on the experiment type
      if (is.label.free == TRUE) {
        cat("Fixing condition", condition,"...\n")  
      } else {
        cat("Fixing", condition,"...\n")
      }  
      # Get the biological/technical replicates vector for a particural
      # to be fixed
      biological.replicates <- replicates.per.condition[[index]]$biological.replicates
      technical.replicates <- replicates.per.condition[[index]]$technical.replicates
      
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
  biological.replicates <- c()
  technical.replicates <- c()
  
  # Iterate over the conditions
  for (condition in replicates.per.condition) {
    
    # Add the biological/technical replicates of each condition to the vector
    biological.replicates <- c(biological.replicates, condition$biological.replicates)
    technical.replicates <- c(technical.replicates, condition$technical.replicates)
  }
  # Finally wrap them in a list
  restored.replicates <- list("biological.replicates" = biological.replicates, "technical.replicates" = technical.replicates)
  
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
           experimental.description <- paste0( experimental.description.biological,
                                               paste0("T", 1),
                                               paste0("F", 1))
         },
         {
           cat("No bioreps, no fractions? Really? \n")
           cat("Invalid experimental description id in make.experimental.description function.\n")  
           return (FALSE)
         },
         {
           cat("We have bioreps and techreps.\n")
           experimental.description <- paste0(experimental.description.biological,
                                              experimental.description.technical,
                                              paste0("F", 1))
         },
         {
           cat("We have bioreps and fractions.\n")
           experimental.description <- paste0(experimental.description.biological,
                                              paste0("T", 1),
                                              experimental.description.fractions)
         },
         {
           stop("No bioreps? Really? Invalid experimental description id in make.experimental.description function.\n")  
         },
         {
           cat("We have bioreps, techreps and fractions.\n")
           experimental.description <- paste0(experimental.description.biological,
                                              experimental.description.technical,
                                              experimental.description.fractions)
         })
  
  return (experimental.description)
}

get.experiment.metadata <- function(experimental.structure) {
  #
  # Makes a list of lists where each element corresponds to a conditions, and the number of biological/technical
  # replicates and fractions for this condition 
  #
  # Args:
  #   experimental.structure: The experimenta structure data.table
  #
  # Retutns:
  #   A list of lists where each element corresponds to a conditions, and the number of biological/technical
  #   replicates and fractions for this condition 
  #
  
  # Make an empty experimental metadata list
  experiment.metadata <- list()
  
  # Get the conditions 
  conditions <- unique(experimental.structure$condition)
  
  # Now traverse through the condition
  for (condition in conditions) {
    
    # For each condition
    # Get the biological replicates 
    unique.biological.replicates <- unique(experimental.structure[condition == condition, biological.replicate])
    
    # Get the technical replicates 
    unique.technical.replicates <-  unique(experimental.structure[condition == condition, technical.replicate])
    
    # Get the fractions
    unique.fractions <- unique(experimental.structure[condition == condition, fraction])
    
    # Now wrap the number of each a list
    condition.metadata <- list("number.of.biological.replicates" = length(unique.biological.replicates),
                               "number.of.technical.replicates" = length(unique.technical.replicates),
                               "number.of.fractions" = length(unique.fractions))
    
    # And add that element on the metadata list
    experiment.metadata[[condition]] <- condition.metadata
  }
  
  return (experiment.metadata)
}

find.proteome.discoverer.protein.column <- function(evidence.data) {
  #
  # Find the protein groups accession ids column from the Proteome Discoverer evidence file 
  #
  # Args:
  #   evidence.data: The evidence file
  #
  # Returns:
  #   The protein group accession id column name or the empty string
  #
  
  # Initialize the protein groups column
  protein.groups.column <- ""
  
  # And find the correct protein groups column depending on the version
  if ("Protein Group Accessions" %in% colnames(evidence.data) == TRUE) {
    protein.groups.column <- "Protein Group Accessions"
  } else if ("Protein Accessions" %in% colnames(evidence.data) == TRUE) {
    protein.groups.column <- "Protein Accessions"
  } else {
    # If the is no such column, abort the scripts execution
    stop("The evidence dataset does not contain the columns 'Protein Group Accessions' or 'Protein Accessions'. Aborting...\n")
  }
  return (protein.groups.column)
}

correct.maxquant.files <- function(protein.groups.data, evidence.data) {
  #
  # Constructs a Protein Names column, if it is abscent, using the FASTA headers,
  # cleans the protein groups of empty lines  and ad the columns id, Protein IDs, Protein Name
  # to the evidence file.
  #
  # Args:
  #   protein.groups.data:  The protein groups file
  #   evidence.data:        The evidence file
  #
  # Returns:
  #   A list with the corrected protein.groups file and the corrected evidence file
  #
  
  # Store the column names of the evidence and the protein groups files
  protein.groups.column.names <- colnames(protein.groups.data)
  evidence.column.names <- colnames(evidence.data)
  
  # Find the number of occurences of the Protein Names/Protein names (depending on the version)
  number.of.protein.names.columns <- length(which(protein.groups.column.names == "protein.names"))
  
  # If there is no Protein Names column
  if (number.of.protein.names.columns == 0) {
    
    # Construnt one using the FASTA headers e.g. the fasta header
    # >PROT123 SWISS-PROT: PROT123 my favourite protein;>PROT123 SWISS-PROT: PROT123 my second favourite protein
    # will only keep the first entry (everything after ';' will be ingored) and it will also remove
    # the >PROT123 part
    protein.groups.data[, protein.names := gsub("(^>[[:alnum:]]+[[:punct:][:blank:]])|(;>.*)|(;[[:alnum:]]+)",
                                                 "",
                                                 protein.groups.data[, fasta.headers],
                                                 perl = TRUE)]
  }
  
  # Reset the column names of the protein groups file
  protein.groups.column.names <- colnames(protein.groups.data)
  
  # Order the table by the Evidence IDs for fast conditional empty line remove
  setkey(protein.groups.data, evidence.ids)
  
  # Remove lines with empty Evidence IDs
  protein.groups.data <- protein.groups.data[!""]
  
  # In the protein groups file does not contain the column Protein IDs, but contains the column Peptide IDs,
  # rename the Peptide IDs to Protein IDs
  if (!"protein.ids" %in% protein.groups.column.names &
       "peptide.ids" %in% protein.groups.column.names)
  {
    setnames(protein.groups.data, "peptide.ids", "protein.ids")
  }
  
  # Subset the protein groups data.table keeping only the Protein IDs, Protein Names and Evidence IDs columns
  protein.groups.subset <- protein.groups.data[ , 
                                               .SD,
                                               .SDcols = c("protein.ids",
                                                            "protein.names",
                                                            "evidence.ids")]
  # Break the each Evidence IDs cell in multiple by the ';' and merge them by their corresponding Protein ID and
  # Protein Name
  protein.groups.subset.multiline.evidence <- cSplit(protein.groups.subset, "evidence.ids", ";", "long")
  
 
  # Rename Evidence IDs column to ID in the multiple evidence table
  setnames(protein.groups.subset.multiline.evidence,
           "evidence.ids",
           "id")
  
  # Change the class of the column from string to integer
  class(protein.groups.subset.multiline.evidence$id) <- 'integer'
  
  # Order the protein.groups.subset.multiline.evidence by the id column
  setkey(protein.groups.subset.multiline.evidence, id)
  
  # Now order the id column in the evidence file
  setkey(evidence.data, id)
  
  # Make a unique code id the if the Protein IDs or Protein Names column exist
  column.delete.case <- (any(grepl("protein.ids", evidence.column.names, perl = TRUE)) * 1) +
                        (any(grepl("protein.names", evidence.column.names, perl = TRUE)) * 2) + 1
  
  # Depending on the code, delete the appropriate columns
  # Unique Codes:
  #   1:  Eveythins is OK
  #   2:  Remove Protein IDs
  #   3:  Remove Protein Names
  #   4   Remove Protein IDs and Protein Names
  #
  switch(column.delete.case,
         cat("No need for column deletion...\n"),
         evidence.data[, "protein.ids" := NULL],
         evidence.data[, "protein.names" := NULL],
         evidence.data[, c("protein.ids", "protein.names") := NULL])
  
  # Reset the column names of the evidence file
  evidence.column.names <- colnames(evidence.data)
 
  # Join the multiline protein groups table with the evidence table in order to
  # make the data.table that we should have in the first place
  merged.evidence.table <- merge(protein.groups.subset.multiline.evidence,
                                 evidence.data,
                                 by = "id",
                                 all = TRUE)
  
  # Order the corrected evidence data.table by id
  setkey(merged.evidence.table, id)
  
  # Wrap the corrected files in a list
  corrected.files <- list("protein.groups.data"= protein.groups.data, "evidence.data" = merged.evidence.table)
  
  return (corrected.files)
}

set.evidence.metadata <- function(raw.file.column, condition.column) {
  #
  # Makes a list with the raw file columns and the condition column 
  #
  # Args:
  #   raw.file.column:  The raw file column string
  #   condition.column: The condition column string
  #
  # Returns:
  #   The list of the metadata
  #
  metadata <- list("raw.file" = raw.file.column, "condition" = condition.column)
  return (metadata)
}

get.evidence.metadata <- function(evidence.columns, data.origin, is.label.free, is.isobaric) {
  #
  # Return a list with metadata (raw file column/ condition column) of the evidence dataset
  #
  # Args:
  #   evidence.columns:   A string vector with the column names of the evidence file
  #   data.origin:        A string 'MaxQuant' or 'Proteome-Discoverer' for the origin of the data
  #   is.label.free:      TRUE or FALSE depending on the experiment type label-free or labeled
  #   is.isobaric:        TRUE or FALSE depending on the experiment type isobaric or not
  #
  # Returns:
  #   A list with the metadata
  #
  
  # Initialize the variables using the analysis.parameters
  origin.is.maxquant <- data.origin == "MaxQuant"
  
  # Set the code depending on the combination of the analysis
  status.code <-  (origin.is.maxquant * 1)  +
                  (is.label.free * 2)       + 1
  
  # 1: Proteome Discoverer and labeled experiment
  # 2: MaxQuant and labeled experiment
  # 3: Proteome Discoverer and label-free experiment OR Isobaric!     ### BEWARE  ###
  # 4: MaxQuant and label-free experiment OR Isobaric!                ### BEWARE  ###
  # TODO check real data!
  
  # Depending on the MaxQuant version there are might be slight 
  # differences in the column names
  if (status.code == 2 | status.code == 4) {
      maxquant.raw.file.column <- "raw.file"
  }
  
  switch( status.code,
          {
            # 1: Proteome Discoverer and labeled experiment
            metadata <- set.evidence.metadata("spectrum.file", "modifications")
          },
          {
            # 2: MaxQuant and labeled experiment
            metadata <- set.evidence.metadata(maxquant.raw.file.column, "labeling.state")
          },
          {
            # 3: Proteome Discoverer and label-free experiment OR Isobaric!     ### BEWARE  ### 
            if (is.isobaric == TRUE) {
             metadata <- set.evidence.metadata("spectrum.file", "modifications")
           } else {
             metadata <- set.evidence.metadata("spectrum.file", "spectrum.file")
           }
          },
          {
            # 4: MaxQuant and label-free experiment OR Isobaric!                ### BEWARE  ###
            # Correct for isobaric
            if (is.isobaric == TRUE) {
              metadata <- set.evidence.metadata(maxquant.raw.file.column, "labeling.state")
            } else {
              metadata <- set.evidence.metadata(maxquant.raw.file.column, maxquant.raw.file.column)
            }
          })
  
  return (metadata)

}

reform.evidence.isobaric.to.label.free <- function(evidence.data, data.origin) {
  #
  # Reforms the evidence.data table in an isobaric experiment turning it from wide format
  # into long format in order to treat it as a label free experiment
  #
  # Args:
  # evidence.data:  The evidence.data table
  # data.origin:    "MaxQuant" or "Proteome Discoverer" depending on the origin of the data
  #
  # Return:
  #   The reformated evidence.data table
  #
  
  # TODO  
  # Get the column names of the file
  evidence.columns <- colnames(evidence.data)
  
  if (data.origin == "MaxQuant") {
    
    # Get the columns names that we will remove
    columns.to.melt <- grep("^Reporter.intensity.[[:digit:]]",
                            evidence.columns,
                            value = TRUE,
                            perl = TRUE)
    
    columns.to.remove <- grep(paste(c("Intensity", columns.to.melt), collapse = "|"), evidence.columns)
    
    # Id.vars are the column names that will be left untouched by the melt function of data.table
    # Essentially the contain all the column names but Intensity and Reporter intensity X where X = 1, 2, 3 etc
    id.vars <- evidence.columns[-c(columns.to.remove)]
    
    # Now melt the data.table to turn it from wide to long format
    evidence.data <- melt(evidence.data,
                          id.vars = id.vars,
                          measure.vars = columns.to.melt,
                          variable.name = "Labeling.State",
                          value.name = "Intensity")
  
    # conditions.labels<<-sub("^X", "Reporter.intensity.", conditions.labels)
    
    # if (AllowLabelRename == T)
    # {
    #   Rename_Array$old_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1", Rename_Array$old_label)
    #   Rename_Array$new_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1", Rename_Array$new_label)
    # }
    # if (AllowLS == T)
    # {
    #   Ls_array$first_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1",  Ls_array$first_label)
    #   Ls_array$second_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1",  Ls_array$second_label)
    # }
    
    is.label.free <- TRUE;
    # filterL_lbl <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1",  filterL_lbl)
    
    # CMBACK    
    # if(RMisused){
    #   RMtagsdata$name <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1",  RMtagsdata$name)
    # }
  } else {
    # Get the columns names that we will remove
    columns.to.melt <- grep("^Abundance..[[:digit:]]*[[:alpha:]]?$",
                            evidence.columns,
                            value = TRUE,
                            perl = TRUE)
    
    columns.to.remove <- grep(paste(c("Intensity", columns.to.melt), collapse = "|"), evidence.columns)
    
    # Id.vars are the column names that will be left untouched by the melt function of data.table
    # Essentially the contain all the column names but Intensity and Reporter intensity X where X = 1, 2, 3 etc
    id.vars <- evidence.columns[-c(columns.to.remove)]
    
    # Now melt the data.table to turn it from wide to long format
    evidence.data <- melt(evidence.data,
                          id.vars = id.vars,
                          measure.vars = columns.to.melt,
                          variable.name = "Modifications",
                          value.name = "Intensity")
    
 
  #   LabelFree<-T;
  #   if (AllowLabelRename == T)
  #   {
  #     Rename_Array$old_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", Rename_Array$old_label)
  #     Rename_Array$new_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", Rename_Array$new_label)
  #   }
  #   if (AllowLS == T)
  #   {
  #     Ls_array$first_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", Ls_array$first_label)
  #     Ls_array$second_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", Ls_array$second_label)
  #   }
  #   filterL_lbl <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", filterL_lbl)
  #   if(RMisused){
  #     RMtagsdata$name <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", RMtagsdata$name)
  #   }
  }
  return (evidence.data)  
}

trim.evidence.data.protein.descriptions <- function(evidence.data, protein.description.column) {
  #
  # Pastes and trims the evidence protein ids  with the appropriate protein description column
  # e.g. 'ABC123' with 'ABC123 [DATABASEID:123 Tax_id=12345 Gene_Symbol=Abc123]...'
  # 
  # Args:
  #   evidence.data:              The evidence.data table
  #   protein.description.column: The appropriate description column depending on the data.origin 'Protein Descriptions'
  #                               for Proteome Discoverer or Proteins Names for MaxQuant. 
  #
  # Returns:
  #   The pasted and trimmed description column
  #
  
  # Subset the evidence data and keep only the Protein IDs column and the Description Column
  evidence.subset <- evidence.data[,
                                   .SD,
                                   .SDcols = c("protein.ids",
                                               protein.description.column)]
  
  # Add an index column for ordering
  evidence.subset[, index := c(1:nrow(evidence.data))]
  
  # Set the key to Protein IDs for later merging
  setkey(evidence.subset, protein.ids)
  
  # Now make a table with the number of the occurences of each protein
  protein.occurences <- evidence.subset[,
                                        .(n=.N),
                                        by=protein.ids]
  
  # Set the key to Protein IDs for later merging
  setkey(protein.occurences, protein.ids)
  
  # Now paste the columns Protein IDs/ Description and separate them with a ' ['
  evidence.subset[, 
                  trimmed.protein.description := do.call( paste,
                                                            c(  .SD,
                                                                sep = " [")), 
                  .SDcols = c("protein.ids",
                              protein.description.column)]
  
  # Keep only the first 60 characters of the description
  evidence.subset[,
                  tmp.trimmed.protein.description:= gsub("^(.{60}).*",
                                                     "\\1",
                                                     as.character(trimmed.protein.description),
                                                     perl = TRUE)]
  
  # And finally paste the ']...' at the end of the trimmed description
  evidence.subset[, 
                  trimmed.protein.description := do.call( paste,
                                                            c( .SD,
                                                               "]...",
                                                               sep = "")), 
                  .SDcols = "trimmed.protein.description"]

  # Order the table and get the trimmed protein descriptions
  new.protein.ids <- protein.occurences[evidence.subset][order(index), trimmed.protein.description]
  
  return (new.protein.ids)
}

add.user.condition.column.to.evidence <- function(evidence.data, conditions.to.compare, conditions.to.raw.files.list,
                                           evidence.metadata, data.origin, is.label.free, is.isobaric) {
  #
  # Adds a condition column to the evidence.data table with the condition name given by the user
  #
  # Args:
  #   evidence.data:                The evidence data table
  #   conditions.to.compare:        A vector of the conditions to compare e.g. c("Wild", "Mutant")
  #   conditions.to.raw.files.list: The list of conditions to raw.files
  #   evidence.metadata:            The evidence metadata regarding the column names for raw.file definition
  #                                 and condition column definition
  #   data.origin:                  The origin of the data 'MaxQuant' or 'Proteome-Discoverer'
  #   is.label.free:                TRUE or FALSE depending on the experiment
  #   is.isobaric:                  TRUE or FALSE depending on the experiment
  # 
  # Returns:
  #   The altered data.table
  #
  
  # Do the data origin from MaxQuant?
  origin.is.maxquant <- data.origin == "MaxQuant"
  
  # Copy the evidence data in order to dodge the data.table bug
  evidence.data.reformed <- copy(evidence.data)
  
  # TODO maybe this is not needed anymore
  # # Make the pattern of the conditions that we want to compare
  # conditions.to.compare.pattern <- paste(conditions.to.compare, collapse = "|")
  # 
  # # Find the ones to keep e.g. if I have the conditions Mutant1, Mutant2, Wild and 
  # # I want to focus on Mutant1 vs Wild
  # conditions.to.keep <- grepl(conditions.to.compare.pattern,
  #                             names(conditions.to.raw.files.list),
  #                             perl = TRUE)
  # 
  # # Keep only the ones I need
  # trimmed.conditions.to.raw.files.list <- conditions.to.raw.files.list[conditions.to.keep]
  # 
  # # Get the names of the conditions
  # trimmed.conditions <- names(trimmed.conditions.to.raw.files.list)
  
  # Make a status code for fast comparisons
  # 1:  # Proteome-Discoverer Isotopic
  # 2:  # MaxQuant Isotopic
  # 3:  # Proteome-Discoverer Label-Free
  # 4:  # MaxQuant Label-Free
  # 5:  --- Illegal State ---
  # 6:  # Proteome-Discoverer Isobaric
  # 7:  # MaxQuant isobaric
  status.code <-  (origin.is.maxquant * 1) +
                  (is.label.free * 2) +
                  (is.isobaric * 3) + 1
  
  
  
  # If we are on an isotopic experiment, there is no need to add the condition column
  if (status.code == 1  | status.code == 2) {
    return (evidence.data.reformed)
  }
  
  # Add the new column and initialize it to ''
  evidence.data.reformed[, "condition":= ""]
  
  # TODO maybe for labeled PD the column should be Quan channel
  # Get the evidence metadata depending on the data origin and the experiment type
  condition.column <- evidence.metadata$condition
  
  # Now add the conditions to the condition column 
  # TODO changed trimmed.conditions to inside the for
  for (condition.to.compare in conditions.to.compare) {
    
    switch( status.code,
            { 
              # Proteome-Discoverer Labeled
              
              stop("add.user.condition illegan case in switch!\n")
            },
            {
              # MaxQuant Labeled
              
              stop("add.user.condition illegan case in switch!\n")
            },
            {
              # Proteome-Discoverer label-free
              
              # Get the raw files
              condition.raw.files <- conditions.to.raw.files.list[[condition.to.compare]]
              
             
              condition.raw.files.pattern <- paste(condition.raw.files, collapse = "|")
              
              # Find which rows should should have the corresponding condition 
              indexes.of.condition <- grepl(condition.raw.files.pattern,
                                            evidence.data.reformed[, get(condition.column)],
                                            perl = TRUE)
              
              # Add the condition on these rows
              evidence.data.reformed[indexes.of.condition, condition := condition.to.compare] 
            },
            {
              # MaxQuant label-free
              
              # Get the raw files
              condition.raw.files <- conditions.to.raw.files.list[[condition.to.compare]]
              
              cat("Case MaxQuant Label free add compare column\n")
              
              condition.raw.files.pattern <- paste(condition.raw.files, collapse = "|")
              
              # Find which rows should should have the corresponding condition 
              indexes.of.condition <- grepl(condition.raw.files.pattern,
                                            evidence.data.reformed[, get(condition.column)],
                                            perl = TRUE)
              
              # Add the condition on these rows
              evidence.data.reformed[indexes.of.condition, condition := condition.to.compare] 
            },
            {
              cat("Incorrect case in add.compare.column.to.evidence...\n")
            },
            {
              # Proteome-Discoverer Isobaric
              
              # Add the condition on these rows
              evidence.data.reformed[,condition := get(condition.column)]
              
            },
            {
             
              # MaxQuant isobaric
              
              # Add the condition on these rows
              evidence.data.reformed[, condition:=get(condition.column)] 
            })
  }
  
  return (evidence.data.reformed)
}

discard.useless.conditions.per.experiment <- function(evidence.data, conditions.to.compare, is.label.free, is.isobaric) {
  #
  # Removes the useless rows for the conditions in which we are not interested, dependingon the experimental setup
  #
  # Args:
  #   evidence.data:          The evidence data.table
  #   conditions.to.compare:  The vector with the condition I want to compare
  #   is.label.free:          TRUE or FALSE depending on if it is an label-free experiment or not.
  #   is.isobaric:            TRUE or FALSE depending on if it is an isobaric experiment or not.
  #
  # Returns:
  #   The cleaned evidence data.table
  #
  
  # Copy the data that i want to clean
  cleaned.evidence.data <- copy(evidence.data)
  
  # Make an id to determine the experimental setup
  case.id <-  (is.label.free * 1) + 
              (is.isobaric * 2) + 1
  
  # By default clean all the blanks and the NAs
  
  # Now depending on the experimental setup remove the appropropriate rows
  switch(case.id,
         {
            
           cleaned.evidence.data <- clear.user.condition.rows(cleaned.evidence.data, "BLANK", "raw.file")
           cleaned.evidence.data <- clear.user.condition.rows(cleaned.evidence.data, "NA", "raw.file")
            
         },
         {
            # Label-free Experiment
           cleaned.evidence.data <- clear.user.condition.rows(cleaned.evidence.data, "BLANK")
           cleaned.evidence.data <- clear.user.condition.rows(cleaned.evidence.data, "NA")
         },
         {
           # Isobaric Labeled Experiment
           cleaned.evidence.data <- clear.user.condition.rows(cleaned.evidence.data, "BLANK")
           cleaned.evidence.data <- clear.user.condition.rows(cleaned.evidence.data, "NA")
          })
  
  return (cleaned.evidence.data)
}

clear.user.condition.rows <- function(evidence.data, condition.to.remove, by.column = "condition") {
  #
  # Clears the evidence data from rows with unassigned condition
  #
  # Args:
  #   evidence.data:        The evidence data table
  #   condition.to.remove:  The condition to remove
  #   by.column:            Default value is "condition" column.  The column that the filtering will be based on
  # Returns:
  #   The cleaned data.table
  #
  
  # Order the table by Condition or Raw.file depending on the experiment
  if (by.column == "condition") {
    setkey(evidence.data, condition)
  } else {
    setkey(evidence.data, raw.file)
  }
  
  switch(condition.to.remove,
         "BLANK" = {
           # Clear the empty condition rows. We put BLANK instead of "" because the empty string is not valid!
           clear.evidence.data <- evidence.data[!""]
         },
         "-1" = {
           # Clear the empty condition rows
           clear.evidence.data <- evidence.data[!"-1"]
         },
         "0" = {
           # Clear the empty condition rows
           clear.evidence.data <- evidence.data[!"0"]
         },
         "1" = {
           # Clear the empty condition rows
           clear.evidence.data <- evidence.data[!"1"]
         },
         "2" = {
           # Clear the empty condition rows
           clear.evidence.data <- evidence.data[!"2"]
         },
         "NA" = {
           if (by.column == "condition") {
             # Clear the NA condition rows
             clear.evidence.data <- na.omit(evidence.data, cols = "condition")
           } else {
             # Clear the NA raw.file rows
             clear.evidence.data <- na.omit(evidence.data, cols = "raw.file")
           }
         })
  
  # And reset the order to the Protein IDs
  setkey(clear.evidence.data, "protein.ids")
  
  return (clear.evidence.data)
}

bring.data.to.common.format <- function(evidence.data, data.origin, is.label.free, is.isobaric) {
  #
  # Bring both evidence files from MaxQuant or Proteome-Discoverer to a common format for analysis
  #
  # Args:
  #   evidence.data:  The evidence.data data.table
  #   data.origin:    Where the data come from, 'MaxQuant' or 'Proteome-Dsicoverer'
  #   is.label.free:  TRUE or FALSE depending on the experiment
  #   is.isobaric:    TRUE or FALSE depending on the experiment
  #
  # Returns:
  #   The subset of the evidence.data table
  #
  
  # Copy the evidence.data
  evidence.data <- copy(evidence.data)
  
  # Get the evidence columns
  evidence.columns <- colnames(evidence.data)
  
  # If data come from the MaxQuant  rename the 'Peptide ID' to 'Unique Sequence ID'
  # If data come form Proteome-Discoverer the column maybe already exists
  if (data.origin == "MaxQuant") {
    
    if ("peptide.id" %in% evidence.columns == TRUE) {
      setnames(evidence.data,"peptide.id", "unique.sequence.id")
    }
    
    # Renew the evidence columns vector
    evidence.columns <- colnames(evidence.data)
    
    # Clean 'Intensity' part from 'Intensity XYZ' columns
    # Not sure...
    if (is.isobaric == FALSE) {
      # colnames(evidence.data) <- sub('intensity\\.(.+)',
      #                                "\\1",
      #                                evidence.columns,
      #                                perl = TRUE)
    } else {
      # yorgodilo had something commented here...
    }
  }
  
  # Renew the evidence columns vector
  evidence.columns <- colnames(evidence.data)
  
  # If data come form Proteome-Discoverer and the the column 'Unique Sequence ID'
  # still does not exist, rename the 'Annotated Sequence'
  # Not sure...
  if ("unique.sequence.id" %in% evidence.columns == FALSE &
      "annotated.sequence" %in% evidence.columns == TRUE) {
    
    setnames(evidence.data, "annotated.sequence", "unique.sequence.id")
    
    evidence[, unique.sequence.id := sub(".*? (.*?) .*",
                                           "\\1",
                                           evidence$unique.sequence.id)]
  }
  
  # TODO check with real data!
  evidence.columns <- colnames(evidence.data)
  
  if ("quan.usage" %in% evidence.columns == FALSE &
      "peptide.quan.usage" %in% evidence.columns == TRUE) {
    setnames(evidence.data, "peptide.quan.usage", "quan.usage")
  }
  
  # Renew the evidence columns vector
  evidence.columns <- colnames(evidence.data)
  
  # Do the data come from MaxQuant?
  origin.is.maxquant <- data.origin == "MaxQuant"
  
  
  # 1:  Proteome-Discoverer Labeled
  # 2:  MaxQuant Labeled
  # 3:  Proteome-Discoverer Label-free
  # 4:  MaxQuant Label-free
  # Make a status code for fast case evaluation
  status.code <-  (origin.is.maxquant * 1) +
                  (is.label.free * 2) + 1
  
  switch( status.code,
          {
            # 1:  Proteome-Discoverer Labeled
            # TODO
            evidence.data.subset <- evidence.data[, .SD,
                                                    .SDcols = c("quan.usage",
                                                                "protein.ids",
                                                                "unique.sequence.id",
                                                                intensity.columns,
                                                                "condition",
                                                                "description")]
          
            setkey(evidence.data.subset, 
                   "description",
                   "protein.ids",
                   "unique.sequence.id")  
          },
          {
            # 2:  MaxQuant Labeled
            intensity.columns <- grep("^intensity", colnames(evidence.data), perl = TRUE, value = TRUE)
            
            
            evidence.data.subset <- evidence.data[, .SD,
                                                  .SDcols = c("protein.ids",
                                                              "unique.sequence.id",
                                                              intensity.columns,
                                                              "condition",
                                                              "description")]
            
            setkey(evidence.data.subset, 
                   "description",
                   "protein.ids",
                   "unique.sequence.id")
            
            # evidence.data.subset <- evidence.data.subset[, intensity.columns[2]:=NULL]
            
            intensity.column.1 <- intensity.columns[1]
            
            # Get maximum PSM intensity per peptide/protein/[(rep_desc/label) = raw_file]
            evidence.data.subset.1 <- evidence.data.subset[, 
                                                         .("intensity.1" = list(get(intensity.column.1))),
                                                         by=.(description,
                                                              protein.ids,
                                                              unique.sequence.id,
                                                              condition)]
            
            setnames(evidence.data.subset.1, "intensity.1", intensity.column.1)
            
            intensity.column.2 <- intensity.columns[2]
            
            # Get maximum PSM intensity per peptide/protein/[(rep_desc/label) = raw_file]
            evidence.data.subset.2 <- evidence.data.subset[, 
                                                           .("intensity.2" = list(get(intensity.column.2))),
                                                           by=.(description,
                                                                protein.ids,
                                                                unique.sequence.id,
                                                                condition)]
            
            setnames(evidence.data.subset.2, "intensity.2", intensity.column.2)
            
            evidence.data.subset <- merge(evidence.data.subset.1,
                                          evidence.data.subset.2)
          },
          {
            # 3:  Proteome-Discoverer Label-free
            
            # TODO Here we should use Precursor Area is unfortunately buggy (sometimes 0/NA), so we are
            # left with Intensity to work with. 
            # intensity.column <- "Precursor Area"            
            intensity.column <- "intensity"
            
            evidence.data.subset <- evidence.data[, .SD,
                                                  .SDcols = c(  "protein.ids",
                                                                "unique.sequence.id",
                                                                intensity.column,
                                                                "condition",
                                                                "description")]
            
            setkey( evidence.data.subset,
                    description,
                    protein.ids,
                    unique.sequence.id,
                    condition)
            
            # Get maximum PSM intensity per peptide/protein/[(rep_desc/label) = raw_file]
            evidence.data.subset <- evidence.data.subset[, 
                                                         .("intensities" = paste(get(intensity.column),
                                                                                 collapse = ";")),
                                                         by=.(description,
                                                              protein.ids,
                                                              unique.sequence.id,
                                                              condition)]
          },
          {
            # 4:  MaxQuant Label-free
            
            intensity.column <- "intensity"
            
            evidence.data.subset <- evidence.data[, .SD,
                                                    .SDcols = c(  "protein.ids",
                                                                  "unique.sequence.id",
                                                                  intensity.column,
                                                                  "condition",
                                                                  "description")]
            
            # evidence.data.subset <- na.omit(evidence.data.subset, "Intensity")
            
            setkey(evidence.data.subset,
                   description,
                   protein.ids,
                   unique.sequence.id,
                   condition)
            
            # Get maximum PSM intensity per peptide/protein/[(rep_desc/label) = raw_file]
            evidence.data.subset <- evidence.data.subset[, 
                                                         # .("Test Intensities" = list(get(intensity.column))),
                                                         .("intensities" = list(get(intensity.column))),
                                                          by=.(description,
                                                               protein.ids,
                                                               unique.sequence.id,
                                                               condition)]

            })
  
  return (evidence.data.subset)
}

zeros.to.nas <- function(evidence.data) {
  #
  #
  #
  #
  #
  #
  #
  #
  
  data <- copy(evidence.data)
  
  # Find the Intensity column names
  intensity.columns <- grep("^intensity", colnames(data), perl = TRUE, value = TRUE)
  
  zeros.A <- which(data[, get(intensity.columns[1])] == 0)
  zeros.B <- which(data[, get(intensity.columns[2])] == 0)
  
  data <- data[zeros.A, eval(intensity.columns[1]) := NA]
  data <- data[zeros.B, eval(intensity.columns[2]) := NA]
  
  return (data)
}

build.analysis.data <- function(protein.groups.data, evidence.data, data.origin, is.isobaric, is.label.free) {
  #
  # Combines the information from "protein.groups.data" and the "evidence.data" tables in one matrix,
  # depending on the experimenta set up and transforms in in a common format
  #
  # Args:
  #   protein.groups.data:  The proteinGroups.txt file from the MaxQuant or NULL if data comes from Proteome-Discoverer
  #   evidence.data:        The evidence.txt data.table for data coming from the MaxQuant
  #                         or the PSM file data.table for data coming from the Proteome Discoverer
  #   data.origin:          String, "MaxQuant" or "Proteome-Discoverer" depending on the data source
  #   is.isobaric:          TRUE or FALSE depending on the experiment
  #   is.label.free:        TRUE or FALSE depending on the experiment
  #
  # Returns:
  #   A tranformed data.table with only the needed column for the analysis
  #
  
  # # TODO Remove test
  # if (project.variables$development.stage == TRUE) {
  #   protein.groups.data <- global.variables$protein.groups.data
  #   evidence.data <- global.variables$evidence.data
  #   data.origin <- global.variables$dataset.origin
  #   is.label.free <- global.variables$is.label.free
  #   is.isobaric <- global.variables$is.isobaric
  # }
  
  # Initialize the protein groups column
  protein.groups.column <- ""
  
  # Initialize evidence column names
  evidence.column.names <- colnames(evidence.data)
  
  # Pick the right protein groups column depending on the data origin
  if (data.origin == "Proteome-Discoverer") {
    
    # Explicit handling for the Proteome Discoverer as there are differences between versions
    protein.groups.column <- find.proteome.discoverer.protein.column(evidence.data)
  } else {
    protein.groups.column <- "proteins"
  }
  
  # Get the index of the protein groups accessions column
  protein.groups.column.position <- which(colnames(evidence.data) == protein.groups.column)
  
  # And rename it as Protein IDs
  evidence.column.names[protein.groups.column.position] <- "protein.ids"
  
  # Reset evidence column names
  evidence.column.names <- colnames(evidence.data)
  
  # If the data come from the MaxQuant correct the evidence and proteinGroups files
  if (data.origin == "MaxQuant") {
    corrected.files <- correct.maxquant.files(protein.groups.data, evidence.data)
    protein.groups.data <- corrected.files$protein.groups.data
    evidence.data <- corrected.files$evidence.data
  }
  
  # Order the evidence.data by the Protein IDs
  setkey(evidence.data, protein.ids)
  
  # And remove any row with empty Protein IDs
  evidence.data <- evidence.data[!""]
  
  # In the case of an isobaric labeled experiment we can treat is as if it was label-free,
  # after some reformating 
  if (is.isobaric == TRUE) {
    
    cat("We have an isobaric label file!\n")
    # Reform the evidence data table
    evidence.data <- reform.evidence.isobaric.to.label.free(evidence.data, data.origin)
    
    # Set the label.free flag to TRUE
    is.label.free <- TRUE
  }
  
  # Reset evidence column names
  evidence.column.names <- colnames(evidence.data)
  
  if (data.origin == "MaxQuant") {
    protein.description.column <- "protein.names"
  } else {
    protein.description.column <- "protein.descriptions"
    if (! protein.description.column %in% evidence.column.names) {
      evidence.data[, protein.descriptions := ""]
    }
  }

  # Paste and trim the evidence protein ids  with the appropriate protein description column
  # e.g. 'ABC123' with 'ABC123 [DATABASEID:123 Tax_id=12345 Gene_Symbol=Abc123]...'
  evidence.data$protein.ids <- trim.evidence.data.protein.descriptions( evidence.data,
                                                                        protein.description.column)
  
  # Store the raw.file column and the condition/label column depending on the data origin
  evidence.metadata <- get.evidence.metadata(colnames(evidence.data),
                                             data.origin,
                                             is.label.free,
                                             is.isobaric)
  
  # Add a column with the user defined conditions to compare
  evidence.data <- add.user.condition.column.to.evidence(evidence.data,
                                                    conditions.to.compare,
                                                    conditions.to.raw.files.list,
                                                    evidence.metadata,
                                                    data.origin,
                                                    is.label.free,
                                                    is.isobaric)
  
  # Clear the data from the rows with specific condition depending on the experiment
  evidence.data <- discard.useless.conditions.per.experiment( evidence.data,
                                                              conditions.to.compare,
                                                              is.label.free,
                                                              is.isobaric)
  
  
  evidence.data <- merge(evidence.data,
                         global.variables$experimental.structure,
                         by.x = c(evidence.metadata$raw.file),
                         by.y = "raw.file")
  
  
  # Finally prepare the evidence to have a common format across experimental setups
  evidence.data <- bring.data.to.common.format( evidence.data,
                                                data.origin,
                                                is.label.free,
                                                is.isobaric)
  return (evidence.data) 
}

