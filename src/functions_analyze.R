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

make.Venn.diagram <- function(evidence.data, conditions.to.compare, analysis.title) {
  #
  # Makes a Venn diagram between 2 conditions
  #
  # Args:
  #   evidence.data:          The data.table with the evidence data
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
  venn.table <- venn.table[, .("Occurences"=.N), by=c( "Protein IDs",
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
