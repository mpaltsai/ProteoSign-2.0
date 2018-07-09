# Sets the current working directory and calls the following scripts:
# initialize.R, load_data.R, pull_data_from_DB.R, build.R, analyze.R

# Clear enviroment
rm(list = ls())

# Set project variables
project.variables <- list("development.stage" = TRUE)

check.packages <- function(packages) {
  # Checks to see if desired packages are installed and if they are not, it installs them.
  # Args:
  #   package: A vector of the desired packages
  # Returns:
  #   TRUE by default
  #
  new.packages <- packages[ !(packages %in% installed.packages()[, "Package"])]
  if ( length(new.packages) )
    install.packages(new.packages, dependencies = TRUE)
  return (TRUE)
}

# Set mirrors and repositories
if ( project.variables[["development.stage"]] == TRUE ) {
  chooseCRANmirror(graphics = FALSE, ind = 1)
  
  chooseBioCmirror(graphics = FALSE, ind = 1)
  
  setRepositories(graphics = FALSE, ind = 1:4)
}

# Project Packages to be installed
packages <- c("here", "devtools")


if (dir.exists("packrat") == FALSE & 
    project.variables[["development.stage"]] == TRUE) {
  
  # Install packrat as the base package manager
  check.packages("packrat")
  
  # # Remove old packrat
  # system("rm -rf packrat/ .Rprofile")
  
  library("packrat")
  
  # Initialize packarat
  init(getwd())
  
  # Set packrat mode ON
  packrat_mode(on = TRUE)
  
  # Install the Project Packages needed for development
  check.packages(packages)
  
  # Load the project packages
  invisible(lapply(packages, require, character.only = TRUE))
  
} else {
  # Load all packages
  invisible(lapply(c("packrat", packages), library, character.only = TRUE))
  
}



if (project.variables[["development.stage"]] == TRUE) {
  # Install TODOR addin for TODOs
  install_github("dokato/todor")
  
  # Take a packrat snapshot
  snapshot()
}

# Scripts to call
files.to.load <- c( "initialize.R",
                    "load_data.R")#,
                    #"build.R",
                    #"analyze.R")


# Set the currenct working directory
setwd(here("src"))


# Run the scripts
successful.run <- unlist(lapply(files.to.load, source))
