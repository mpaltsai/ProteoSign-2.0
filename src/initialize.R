# This file loads all the packages,libraries and data needed regarding the
# workspace, loads the functions.R script and sets the global variables.

# Clear enviroment
rm(list = grep("^project.variables|^check.packages", ls(), value = TRUE, invert = TRUE))

global.variables <- list("quantitation.type" =              "Proteins",
                         "replicate.multiplexing.is.used" =  FALSE,
                         "dataset.origin" =                  "MaxQuant",
                         "is.label.free" =                   TRUE)


# Project Packages to be installed
cran.packages <- c("data.table")

# Add packages used during development only
if(project.variables[["development.stage"]] == TRUE) {
  cran.packages <- c(cran.packages, "rbenchmark")
  
  # Temporary reset current working directory
  # in order to work packrat package installation
  setwd(here())
  
  # Install missing packages
  check.packages(cran.packages)
  
  # Reset current working
  setwd(here("src"))
}

# Load all packages
invisible(lapply(cran.packages, library, character.only = TRUE))

# Save loaded packages in packrat
if (project.variables[["development.stage"]] == TRUE) {
  
  # Temporary reset current working directory
  # in order to work packrat package installation
  setwd(here())
  
  # Check the packrat packages status
  packages.status <- status()
  
  # Are the packages in the last snapshot, the same as in the local downloaded packages?
  packages.are.up.to.date <- all(packages.status$packrat.version == packages.status$library.version)
  
  # If not, snapshot the loaded packages
  if (packages.are.up.to.date == FALSE) {
    snapshot()  
  }
  
  # Reset current working
  setwd(here("src"))
}

source("functions.R")

cat("========== End of initialize.R ==========\n")
