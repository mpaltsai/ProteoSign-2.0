make.limma.folder <- function() {
  #
  # Makes the limma output folder 
  #
  # Args:
  #
  # Returns:
  #   Prints a message if the folder was created successfully, otherwise terminates the analysis.
  #
  # Set limma output path
  data.output.path <- here("data-output")
  
  # Set limma output folder name
  limma.folder <- "limma-output"
  
  limma.output <- paste(data.output.path, limma.folder, sep = "/")
  # If folder already exist, remove it
  if( file.exists(limma.output) == TRUE) {
    unlink(limma.output, recursive=T, force=T)
  }
  
  # And then create the new folder
  tryCatch(
    {
      dir.create(limma.output)
      success.message <- paste0("Folder ", limma.folder, " was created under ", data.output.path, "/!\n")
      cat(success.message)
    },
    error = function(cond) {
      stop("Error in make.limma.folder:", cond$message)
    },
    warning = function(cond) {
      stop("Warning in make.limma.folder:", cond$message)
    })
}

make.Venn.diagram <- function(analysis.data, conditions.to.compare, analysis.title) {
  #
  # Makes a Venn diagram between 2 conditions
  #
  # Args:
  #   analysis.data:          The data.table with the evidence data
  #   conditions.to.compare:  A vector with the 2 conditions to compare
  #   analysis.title:         The analysis title provided by the user
  #
  # Returns:
  #   A venn diagram in png format and the data needed for its construction
  #
  
  # Subset the evidence data keeping only the "Protein IDs", "Condition","biological replicate",
  # "technical replicate", "fraction" columns
  venn.table <- evidence.data[,
                              .SD,
                              .SDcols =  c("Protein IDs",
                                           "Condition",
                                          "biological replicate",
                                          "technical replicate",
                                          "fraction")]
  
  # Order the subseted data.table
  setkey(venn.table,  "Protein IDs",
         "Condition",
         "biological replicate",
         "technical replicate",
         "fraction")
  
  # Set limma output path
  data.output.path <- here("data-output")
  
  # Set limma output folder name
  limma.folder <- "limma-output"
  
  # Prepare the limma output path
  limma.output <- paste(data.output.path, limma.folder, sep = "/")
  
  # Change directory to limma output
  setwd(limma.output)
  
  # Count the occurences of each protein in each condition
  venn.table <- venn.table[, .("Occurences"=.N), 
                           by=c( "Protein IDs",
                                  "Condition")]
  
  # Order the data.table
  setkey(venn.table, Condition)
  
  # Write the data in a TSV file
  write.table(venn.table,
              file = paste0(analysis.title, 
                          "-venn-data.txt"),
              sep="\t",
              row.names = FALSE,
              quote = FALSE)
  
  # Now store the 2 conditions
  conditions <- unique(venn.table[, Condition])
  
  # Store in the variable the condition A
  condition.1.table <- venn.table[conditions[1]]
  
  # Store in the variable the condition B
  condition.2.table <- venn.table[conditions[2]]
  
  # Take the Protein IDs and make the result a vector
  condition.1.proteins <- unlist(condition.1.table[, "Protein IDs"])
  condition.2.proteins <- unlist(condition.2.table[, "Protein IDs"])
  
  # Make the set
  sets <- list(condition.1.proteins, condition.2.proteins)
  
  # Give the names of the conditions to the list's set
  names(sets) <- conditions
  
  # By default the VennDiagram package makes a log file, that we want to supperss
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  
  # Now draw the venn diagram
  venn.diagram <- venn.diagram(sets,
                               "test.venn.png",
                               imagetype = "png",
                               main = paste0("Venn diagram between the conditions ",
                                            conditions[1],
                                            " and ",
                                            conditions[2],
                                            "."),
                               category = conditions,
                               scaled = FALSE,
                               alpha = c(0.7, 0.7),
                               fill = c("blue", "red"),
                               cat.default.pos = "text",
                               cat.pos = c(1,1),
                               cex =  2,
                               cat.cex = 2,
                               cat.dist = c(0.05, 0.05))
  
  # And finally reset the working directory to the scripts folder
  setwd(here("src"))  
}

filter.out.reverse.and.contaminants <- function(analysis.data) {
  #
  # Removes Contamintants and Reverse flagged peptides from the evidence data
  # 
  # Args:
  #   analysis.data: The analysis.data 
  #
  # Returns:
  #   The cleaned analysis.data data.table
  #
  
  # Copy the analysis data
  data <- copy(analysis.data)
  
  # Clean them from the Protein IDs starting with CON__ regarding the contaminants
  data.no.contaminants <- subset( data,
                                  grepl("^CON__",
                                        `Protein IDs`,
                                        perl = TRUE) == FALSE)
  
  # Clean them from the Protein IDs starting with REV__ regarding the Reverse sequences
  data.no.contaminants.no.reverse <- subset(data.no.contaminants,
                                            grepl("^REV__",
                                            `Protein IDs`,
                                            perl = TRUE) == FALSE)
  
  return (data.no.contaminants.no.reverse)
}

use.peptides.median <- function (intensities, minimum.detections = 2) {
  intensities <- unlist(intensities)
  intensities <- intensities[!is.na(intensities)]
  if (length(intensities) < minimum.detections) {
    return (list(NA))
  } else {
    return (list(median(intensities)))  
  }
  
}

do.vsn.normalization <- function(filtered.data, condition.names, minimum.detections = 2) {
  #
  # Does variance stabilizing normalization (VSN) on the peptides intensities
  # 
  # Args:
  #   filtered.data:      The filtered data data.table
  #   condition.names:    The conditions to compare
  #   minimum.detections: Default is 2. A strictness parameter regarding how many times is a peptide measured
  #                       If it is measures less than N times, the peptide is considered as detected but not measures
  #                       hence its value is set to NA
  #
  # Returns:
  #   The normalized analysis data.table
  #
  #
  
  # Initialize the variable to return 
  vsn.normalized.data <- NULL
  
  # Make a copy of the analysis data
  data <- copy(filtered.data)
  
  # Keep only the information regarding the biological replicate for each peptide
  data[, description:= gsub("^B|T.*$","", description, perl = TRUE)]
  
  # Rename the "description" column to "Biological Replicate"
  setnames(data, "description", "Biological Replicate")
  
  # Reorder the data.table
  setkey(data, `Biological Replicate`,`Protein IDs`,`Unique Sequence ID` )
  condition.names <- global.variables$conditions.to.compare
  # Iterate over conditions
  for (condition in condition.names)  {
    
    # condition <- global.variables$conditions.to.compare[1]
    # Replace the NaN values with NA
    set(data, which(is.nan(data[, eval(condition)])), condition, NA)
    
    # Transfrorm the data from long format to wide where we have the peptides on the y axis and the biological
    # replicates on the x axis
    wide.data <- dcast(data,
                   `Protein IDs` + `Unique Sequence ID` ~ `Biological Replicate`,
                   value.var = condition,
                   fun.aggregate = use.peptides.median)
    
    # Get only the part of the data.table with the intensities of the peptides
    vsn.matrix <- as.matrix(wide.data[, .SD, .SDcols = -c(1, 2)])
    
    # Do the normalization and suppress the warnings regardind the NA removal the the vsn package throws by default
    suppressWarnings(vsn.matrix.normalized <- justvsn(vsn.matrix))
    
    # Replace the old values of the matrix with the new
    for (index in 1:ncol(vsn.matrix)) {
      wide.data[, index + 2] <- vsn.matrix.normalized[, index]  
    }
    
    # Return the data.table to long format
    long.data <- melt(wide.data,
                 id.vars = c( "Protein IDs", "Unique Sequence ID"),
                 variable.name = "Biological Replicate",
                 value.name = condition)
    
    # Stupid renaming workaround to set the condition name as the next command had trouble to get the name
    # from the condition variable
    setnames(long.data, condition, "V1")
    
    # Reorder the columns to match the format of the initial data.table
    setcolorder(long.data,
             c("Biological Replicate",
             "Protein IDs",
             "Unique Sequence ID",
             "V1"))
    
    # End of stupid renaming workaround to set the condition name as the next command had trouble to get 
    # the name from the condition variable
    setnames(long.data, "V1", condition)
    
    # If the return table is empty just copy the first data.table, otherwise add the second condition column
    if (is.null(vsn.normalized.data) == TRUE) {
      vsn.normalized.data <- copy(long.data)
    } else {
      vsn.normalized.data[, eval(condition) := long.data[, get(condition)]]
    }
  }
  
  return (vsn.normalized.data)
}
