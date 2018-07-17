# Sets the current working directory and calls the following scripts:
# initialize.R, load_data.R, pull_data_from_DB.R, build.R, analyze.R

# Clear enviroment
rm(list = ls())

# Set project variables
project.variables <- list("development.stage" = TRUE)

check.packages <- function(packages, origin.is.github = FALSE) {
  #
  # Checks to see if desired packages are installed and if they are not, it installs them.
  #
  # Args:
  #   package:          A vector of the desired packages
  #   origin.is.github: Default is FALSE. If TRUE the origin of the package is github
  #
  # Returns:
  #   TRUE by default
  #
  
  if (origin.is.github == TRUE) {
    
    # For github packages each package is in the form repositoryName/packageName,
    # so in order to see if the package exists, I split the input on the slash
    # and get every second element
    splited.packages <- unlist(strsplit(packages, "/"))[c(FALSE, TRUE)]
    
    new.packages <- packages[ !(splited.packages %in% installed.packages()[, "Package"])]
  } else {
    new.packages <- packages[ !(packages %in% installed.packages()[, "Package"])]
  }
  
  if (length(new.packages) > 0) {
    if (origin.is.github == TRUE) {
      install_github(new.packages)
    } else {
      install.packages(new.packages, dependencies = TRUE)
    }
  }  
  return (TRUE)
}

# Project Packages to be installed
cran.packages <- c("here", "devtools")
github.packages <- c("dokato/todor")

if (dir.exists("packrat") == FALSE & 
    project.variables[["development.stage"]] == TRUE) {
  
  # Set mirrors and repositories
  chooseCRANmirror(graphics = FALSE, ind = 1)
  
  chooseBioCmirror(graphics = FALSE, ind = 1)
  
  setRepositories(graphics = FALSE, ind = 1:4)
  
  # Install packrat as the base package manager
  check.packages("packrat")
  
  # Load Packrat
  library(packrat)
  
  # # Remove old packrat
  # system("rm -rf packrat/ .Rprofile")
  
  
  # Initialize packarat
  init(getwd())
  
  # Set packrat mode ON
  packrat_mode(on = TRUE)
  
  # Install the Project CRAN Packages needed for development
  check.packages(cran.packages)
  
  # Load the project CRAN packages
  invisible(lapply(cran.packages, require, character.only = TRUE))
  
  # Load github packages
  check.packages(github.packages, origin.is.github = TRUE)
  
  # Take a packrat snapshot
  snapshot()
  
} else {
  
  library("packrat")
  
  # Set packrat mode ON
  packrat_mode(on = TRUE)
  
  # Load all the packages (packrat, here, devtools)
  invisible(lapply(cran.packages, library, character.only = TRUE))
  
}

# Scripts to call
files.to.load <- c( "initialize.R",
                    "load_data.R",
                    "build.R")#,
                    #"analyze.R")

# Set the currenct working directory
setwd(here("src"))

# Try to run the scripts, if something is wrong, stop the analysis with an error message
tryCatch(lapply(files.to.load, source),
         error = function(error) {
           cat(error[[1]],"\n")
           }
         )
