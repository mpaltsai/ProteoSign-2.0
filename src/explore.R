
addLabel<-function(lblname, lbl.Modifications, LabelFree = T){
  
  #levellog("", change=1)
  #If label name is a number some routines won't work, it has to be converted to some acceptable variable name
  lblname<-make.names(lblname)
  labeltxt <- "label";
  if(!LabelFree){
    lbl.Modifications<-gsub("\\(","\\\\(",lbl.Modifications)
    lbl.Modifications<-gsub("\\)","\\\\)",lbl.Modifications)
    unmod_idx<-which(lbl.Modifications == "")
    if(length(unmod_idx) > 0){
      rest_idx<-which(lbl.Modifications[-unmod_idx] != "")
      if(length(rest_idx) > 0){
        lbl.Modifications<-c(lbl.Modifications[unmod_idx],paste(lbl.Modifications[-unmod_idx],"\\)",sep=""))
      }
    }else{
      lbl.Modifications<-paste(lbl.Modifications,"\\)",sep="")
    }
  }else{
    labeltxt <- "condition";
  }
  cat("condition labels", conditions.labels, "<<<\n")
  lblname_i<-which(grepl(paste("^",lblname,"$",sep=""),conditions.labels))
  print(lblname_i)
  if(length(lblname_i) != 0){
    cat(paste("addLabel: Error adding ",
              labeltxt," '",
              lblname,
              "': An existing ",
              labeltxt,
              " with name '",
              lblname,
              "' (specification: ",
              paste(unlist(conditions.labels.Modifications[lblname_i]),
                    collapse=", "),
              ") already exists. Please try a different name.",
              sep=""))
    return(FALSE)
  }
  
  conditions.labels<<-c(conditions.labels, lblname)
  j<-length(conditions.labels.Modifications)+1
  conditions.labels.Modifications[[j]]<<-lbl.Modifications
  nConditions<<-length(conditions.labels)
  #cat("", change=-1)
}

clearLabels<-function(){
  #levellog("", change=1)
  conditions.labels<<-c()
  conditions.labels.Modifications<<-list()
  nConditions<<-length(conditions.labels)
  #levellog("", change=-1)
}

conditions.to.raw.files <- list()

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
  
  #TODO labeled experiments
  
  # Initialize experimental definition. We have conditions for label-free experiments and labels for labeled experiments
  experimental.definition <- ""
  
  # Get the length of the conditions.to.raw.files list
  conditions.number <- length(conditions.to.raw.files)
  
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
             conditions.to.raw.files[[condition]] <- raw.files
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
  
  return (conditions.to.raw.files)
}
