# This file loads all the packages,libraries and data needed regarding the
# workspace, loads the functions.R script and sets the global variables.

# Clear enviroment
rm(list = grep("^global.variables", ls(), value = TRUE, invert = TRUE))

global.variables[["development.stage"]] <- TRUE

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

# Project Packages to be installed
packages <- c("data.table")

# Add packages used during development only
if(development.enviroment == TRUE) {
  packages <- c(packages, "rbenchmark")
}

# Install missing packages
check.packages(packages)

# Load all packages
invisible(lapply(packages, require, character.only = TRUE))

source("functions.R")