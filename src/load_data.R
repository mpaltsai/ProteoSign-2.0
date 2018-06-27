# This file loads all the csv/txt/xlsx/RDS/etc files needed and displays which
# files have been created/loaded in the workspace.

# Clear enviroment
rm(list = grep("^global.variables", ls(), value = TRUE, invert = TRUE))

# Which software do the data come from
if( global.variables[["is.proteome.discoverer.data"]] == FALSE) {
  cat("Data origin: MaxQuant.\n")
} else {
  cat("Data origin: Proteome Discoverer.\n")
}

# Path to the experimental structure file
experimental.structrure.file <- paste(here(), "data-input", "experimental-structure.csv", sep = "/")

# Path to the experimental structure file
evidence.file <- paste(here(), "data-input", "evidence.txt", sep = "/")

cat("Reading experimental structure file...\n")

# Read the experimental structure file
experimental.structure.table <- read.csv(experimental.structrure.file)

cat("Experimental structure file loaded!\n")

cat("Reading evidence file: ", evidence.file,"...\n", sep = "")

# Read the evidence data
evidence.data <- fread(evidence.file, integer64 = "numeric")

cat("Evidence file loaded!\n")

# Add them in the global variables list
global.variables[["evidence.data"]] <- evidence.data
global.variables[["protein.groups.data"]] <- protein.groups.data

# If data come from MaxQuant Software, find the proteinGroups file and read it
if( global.variables[["is.proteome.discoverer.data"]] == FALSE) {
  protein.groups.file <- paste(here(), "data-input", "proteinGroups.txt", sep = "/")
  cat("Reading proteinGroups file: ", protein.groups.file,"...\n", sep = "")
  protein.groups.data <- fread(protein.groups.file, integer64 = "numeric")
  cat("ProteinGroups file loaded!\n")
}


