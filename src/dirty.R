##########################  CLEAN  #########################################
logger <- NULL

initialize.logger <- function(){
  logger <<- create.logger(level = "DEBUG")
}

do.logging <- function(message, level){
  memory <- humanReadable(mem_used(), standard = "SI")
  message <- paste(message," ", memory, ".", sep = "")
  levellog(logger, level, message)
}

