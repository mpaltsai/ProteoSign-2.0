# This file loads all the csv/txt/xlsx/RDS/etc files needed and displays which
# files have been created/loaded in the workspace.

# Clear enviroment and only keep functions and global/project variables
rm(list = grep(paste(c("^global.variables",
                       "^project.variables",
                       lsf.str()),
                     collapse = "|"),
               ls(),
               value = TRUE,
               invert = TRUE))

# Which software do the data come from
if( global.variables[["dataset.origin"]] == "MaxQuant") {
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
experimental.structure.table <- read.csv(experimental.structrure.file,
                                         stringsAsFactors = FALSE,
                                         check.names = FALSE)

# Add the experimental structrue to the global.variables list
global.variables[["experimental.structure"]] <- experimental.structure.table

cat("Experimental structure file loaded!\n")

cat("Reading evidence file: ", evidence.file,"...\n", sep = "")

# Fast cleaning of the evidence file using the tr unix command
cleaning.command <-  paste("tr -d \'\"\\\\\"\' <", evidence.file)

# Read the evidence data
evidence.data <- fread(cleaning.command, integer64 = "numeric")

cat("Evidence file loaded!\n")

# Add the evidence data to the global variables list
global.variables[["evidence.data"]] <- evidence.data

# If data come from MaxQuant Software, find the proteinGroups file and read it
if( global.variables[["dataset.origin"]] == "MaxQuant") {
  
  # The proteinGroups  file path
  protein.groups.file <- paste(here(), "data-input", "proteinGroups.txt", sep = "/")
  
  cat("Reading proteinGroups file: ", protein.groups.file,"...\n", sep = "")
  
  # Fast cleaning of the proteinGroups file using the tr unix command
  cleaning.command <-  paste("tr -d \'\"\\\\\"\' <", protein.groups.file)
  
  # Read the proteinGroups data
  protein.groups.data <- fread(cleaning.command, integer64 = "numeric")
  
  cat("ProteinGroups file loaded!\n")
  
  # Store it in a global variable
  global.variables[["protein.groups.data"]] <- protein.groups.data
  
}


