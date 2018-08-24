keep.only.specific.timestamps.or.cultures <- function(evidence.data, 
                                                      timestamp.to.keep = NA, subset.to.keep = NA) {
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
  if (is.na(timestamp.to.keep) == FALSE) {
    # Make the pattern of the timestamp to keep
    timestamp.to.keep.pattern <- paste0(timestamp.to.keep, ".*")
    
    # Find the pattern
    timestamps.to.clean <- grep(timestamp.to.keep.pattern,
                                data.to.clean$`Raw file`,
                                perl = TRUE)
    
    # And keep only these rows for the specific pattern
    data.cleaned <- data.to.clean[timestamps.to.clean,]
  }
  
  
  # Keep only the wanted subset of the proteome
  if (is.na(subset.to.keep) == FALSE) {
    # Make the pattern of the timestamp to keep
    subset.to.keep.pattern <- paste0(subset.to.keep, ".*")
    
    # Find the pattern
    subset.to.clean <- grep(subset.to.keep.pattern,
                                data.to.clean$`Raw file`,
                                perl = TRUE)
    
    # And keep only these rows for the specific pattern
    data.cleaned <- data.to.clean[subset.to.clean,]
  }
  
  return (data.cleaned)
}

remove.and.rename.raw.files <- function(evidence.data, raw.files.to.remove, raw.files.to.rename) {
  
  data.to.clean <- copy(evidence.data)
  
  if (length(raw.files.to.remove) != length(raw.files.to.rename)) {
    stop("Invalid raw.files.to.remove, raw.files.to.rename lengths. The 2 vectors should have the same length.\n")
  }
  
  for (index in length(raw.files.to.remove)) {
    
    raw.file.to.remove.pattern <- paste0("^", raw.files.to.remove[index], "$")
    
    # Find the pattern
    raw.files.to.keep.positions <- grep(raw.file.to.remove.pattern,
                                        data.to.clean$`Raw file`,
                                        perl = TRUE,
                                        invert = TRUE)
    
    # And keep only these rows for the specific pattern
    data.cleaned <- data.to.clean[raw.files.to.keep.positions,]
    
    
    raw.file.to.rename.pattern <- paste0("^", raw.files.to.rename[index], "$")
    
    data.cleaned[, `Raw file`:=gsub(raw.file.to.rename.pattern, raw.files.to.remove[index], `Raw file`)]
  }
  
  return (data.cleaned)
}