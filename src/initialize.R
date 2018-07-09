# This file loads all the packages,libraries and data needed regarding the
# workspace, loads the functions.R script and sets the global variables.

# Clear enviroment
rm(list = grep("^global.variables", ls(), value = TRUE, invert = TRUE))

global.variables <- list("development.stage" = TRUE,
                         "quantitation.type" = "Proteins",
                         "replicate.multiplexing.is.used" = FALSE,
                         "is.proteome.discoverer.data" = FALSE)

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
if(global.variables[["development.stage"]] == TRUE) {
  packages <- c(packages, "rbenchmark")
}

# Temporary reset current working directory
# in order to work packrat package installation
setwd(here())

# Install missing packages
check.packages(packages)

# Reset current working
setwd(here("src"))

# Load all packages
invisible(lapply(packages, require, character.only = TRUE))

# Save loaded packages in packrat
snapshot()

source("functions.R")