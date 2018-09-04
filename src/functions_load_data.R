add.analysis.parameters.to.global.variables <- function(analysis.metadata) {
  #
  # Makes the global.variables list with all the needed variables that I want to
  # be anytime reachable
  #
  # Args:
  #   analysis.metadata: The analysis parameters, alongside with the values of 
  #   each parameter
  #
  # Returns:
  #   The global.variables list
  #
  
  # Turn the analysis parameters value vector to list
  global.variables <- as.list(analysis.metadata[,2])
  
  # Add the names to the list
  names(global.variables) <- analysis.metadata[,1]
  
  # Make a vector with the boolean variables
  boolean.variables <- c("replicate.multiplexing.is.used",
                         "is.label.free",
                         "is.isobaric")
  
  # Make a vector with the numeric variables
  numeric.variables <- c("plots.format",
                         "minimum.peptide.detections",
                         "minimum.peptides.per.protein")
  
  # Make a vector with the double variables
  double.varables <- c("min.valid.values.percentance",
                       "fold.change.cut.off",
                       "FDR")
  
  # Now turn them to their correct variable type
  global.variables[boolean.variables] <- as.logical(global.variables[boolean.variables])
  global.variables[numeric.variables] <- as.numeric(global.variables[numeric.variables])
  global.variables[double.varables] <- as.double(global.variables[double.varables])
  
  # And finally split and unlist the conditions to compare string
  conditions.to.compare <- unlist(strsplit(global.variables[["conditions.to.compare"]],
                                           split = ","))
  
  # Trim them in case of whitespaces
  global.variables[["conditions.to.compare"]] <- trimws(conditions.to.compare)
  
  return (global.variables)
  
}

trim.and.lowercase.column.names <- function(old.column.names) {
  #
  # Trims the column names from dots and whitespaces and lowercases
  # for compatibility across softares and versions
  #
  # Args:
  #   old.column.names: The column names in the original dataset
  #
  # Returns:
  #   The trimmed and lowercased column names
  #
  
  # Trim the column names from dots and whitespaces
  trimmed.column.names <- unlist(lapply(old.column.names,
                                        gsub,
                                        pattern = "[[:space:].]+",
                                        replacement = "."))
  
  # Now lowercase the column names
  trimmed.and.lowercased.column.names <- tolower(trimmed.column.names)
  
  return (trimmed.and.lowercased.column.names)
}

keep.only.specific.timestamps.or.cultures <- function(evidence.data, 
                                                      timestamp.to.keep = "", subset.to.keep = "") {
  #
  # Does filtering of the raw file column based on the wanted timestamp or subset
  #
  # Args:
  #   evidence.data:      The evidence data.table
  #   timestamp.to.keep:  Default is NA. The timestamp i want to investigate
  #   subset.to.keep:     Default is NA. The subproteome on which i want to focus
  #
  # Returns:
  #   The cleaned evidence data.table with only the rows corresponding to specific timestamp and/or subset
  #
  
  # Make a copy of the original data
  data.to.clean <- copy(evidence.data)
  
  # Keep only wanted timestamp
  if (timestamp.to.keep != "") {
    # Make the pattern of the timestamp to keep
    timestamp.to.keep.pattern <- paste0(".*", timestamp.to.keep, ".*")
    
    # Find the pattern
    timestamp.to.keep.positions <- grep(timestamp.to.keep.pattern,
                                data.to.clean$raw.file,
                                perl = TRUE)
    
    # And keep only these rows for the specific pattern
    data.to.clean <- data.to.clean[timestamp.to.keep.positions,]
  }
  
  # Keep only the wanted subset of the proteome
  if (subset.to.keep != "") {
    # Make the pattern of the timestamp to keep
    subset.to.keep.pattern <- paste0(".*", subset.to.keep, ".*")
    
    # Find the pattern
    subset.to.keep.positions <- grep(subset.to.keep.pattern,
                                data.to.clean$raw.file,
                                perl = TRUE)
    
    # And keep only these rows for the specific pattern
    data.to.clean <- data.to.clean[subset.to.keep.positions,]
  }
  
  # Now our data are cleaned
  data.cleaned <- data.to.clean
  
  return (data.cleaned)
}

remove.and.rename.raw.files <- function(evidence.data, raw.files.to.remove, raw.files.to.rename) {
  #
  # Removes duplicate raw files in case of problematic run, and keeps only
  # the latter run depending on the input e.g. for some reason rawfile1 was
  # run 2 times due to electric sortage so the evidence file has a
  # rawfile1.1 run and a rawfile1.2 run from which we want to hold only the
  # rawfile1.2.
  #
  # Args:
  #   evidence.data:        The evidence data.table
  #   raw.files.to.remove:  The vector of the raw files we want to remove
  #   raw.files.to.rename:  The vector of the second-run raw files   
  #
  # Returns:
  #   The clean evidence file with only the needed raw files
  #
  
  # Get the evidence data
  data.to.clean <- copy(evidence.data)
  
  # If the length of the input vectors is different, stop the execution
  # of the analysis
  if (length(raw.files.to.remove) != length(raw.files.to.rename) |
      raw.files.to.remove == "") {
    stop("Invalid raw.files.to.remove, raw.files.to.rename lengths. The 2 vectors should have the same length.\n")
  }
  
  # Now for each raw file we want to remove
  for (index in length(raw.files.to.remove)) {
    
    # Make the pattern e.g. "^raw.file$"
    raw.file.to.remove.pattern <- paste0("^", raw.files.to.remove[index], "$")
    
    # Find the rows in which the pattern we want to remove, does not occures
    raw.files.to.keep.positions <- grep(raw.file.to.remove.pattern,
                                        data.to.clean$raw.file,
                                        perl = TRUE,
                                        invert = TRUE)
    
    # And keep only these rows for the specific pattern
    data.cleaned <- data.to.clean[raw.files.to.keep.positions,]
    
    # Now make the pattern we want to rename to the previous raw.file
    raw.file.to.rename.pattern <- paste0("^", raw.files.to.rename[index], "$")
    
    # Finally rename the each occurence of "raw.file1.2" to "raw.file1.1"
    data.cleaned[, raw.file := gsub(raw.file.to.rename.pattern, raw.files.to.remove[index], raw.file)]
  }
  
  return (data.cleaned)
}