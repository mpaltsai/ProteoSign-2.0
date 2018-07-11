make.limma.folder <- function() {
  #
  # TODO error handling?
  # Makes the limma output folder 
  #
  # Args:
  #
  # Returns:
  #   TRUE if the folder was created successfully, otherwise FALSE.
  #
  # Set limma output path
  data.output.path <- here("data-output2")
  
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
      return (TRUE)
    },
    error = function(cond) {
      cat("Error in make.limma.folder:", cond$message)
      return (FALSE)
    },
    warning = function(cond) {
      cat("Warning in make.limma.folder:", cond$message)
      return (FALSE)
    })
}
