# Sets the current working directory and calls the following scripts:
# initialize.R, load_data.R, pull_data_from_DB.R, build.R, analyze.R

# Clear enviroment
rm(list = ls())

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
packages <- c("here")

# Install missing packages
check.packages(packages)

# Load all packages
invisible(lapply(packages, require, character.only = TRUE))

#Do I need to use a database for data retrieve
no.DB <- TRUE

# Set it as global so other scipts can also see this
global.variables <- list("no.DB" = TRUE)

# Scripts to call
files.to.load <- c( "initialize.R",
                    "load_data.R",
                    "pull_data_from_DB.R",
                    "build.R",
                    "analyze.R")

# Remove database script if not used
if (no.DB == TRUE) {
  files.to.load <- files.to.load[-c(3)]
}

# Set the currenct working directory
setwd(here("src"))

# Run the scripts
successful.run <- unlist(lapply(files.to.load, source))
