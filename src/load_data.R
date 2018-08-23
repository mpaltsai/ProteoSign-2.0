# This file loads all the csv/txt/xlsx/RDS/etc files needed and displays which
# files have been created/loaded in the workspace.

# Clear enviroment and only keep functions and global/project variables
rm(list = grep(paste(c("^project.variables",
                       lsf.str()),
                     collapse = "|"),
               ls(),
               value = TRUE,
               invert = TRUE))

# Return the memory to the OS
gc(verbose = FALSE,
   reset = TRUE)

# Path to the analysis metadata file
analysis.metadata.file <- paste(here(), "data-input/test-case", "analysis-metadata.csv", sep = "/")

# Read the analysis metadata file
analysis.metadata <- read.csv(analysis.metadata.file,
                              stringsAsFactors = FALSE,
                              check.names = FALSE)

# Make an empty global.variables list to store any global variable
global.variables <- list()

# Initialize the global variables with the analysis parameters
global.variables[["replicate.multiplexing.is.used"]] <- analysis.metadata$replicate.multiplexing.is.used
global.variables[["dataset.origin"]] <- analysis.metadata$dataset.origin
global.variables[["is.label.free"]]  <- analysis.metadata$is.label.free
global.variables[["is.isobaric"]]    <- analysis.metadata$is.isobaric


# Which software do the data come from
if( global.variables[["dataset.origin"]] == "MaxQuant") {
  cat("Data origin: MaxQuant.\n")
} else {
  cat("Data origin: Proteome Discoverer.\n")
}

is.label.free <- global.variables[["is.label.free"]]

# Path to the experimental structure file
experimental.structrure.file <- paste(here(), "data-input/test-case", "experimental-structure.csv", sep = "/")

# Path to the experimental structure file
evidence.file <- paste(here(), "data-input/test-case", "evidence.txt", sep = "/")

# Path to the proteinGroups  file 
protein.groups.file <- paste(here(), "data-input/test-case", "proteinGroups.txt", sep = "/")

cat("Reading experimental structure file...\n")

# Read the experimental structure file
experimental.structure.table <- read.csv(experimental.structrure.file,
                                         stringsAsFactors = FALSE,
                                         check.names = FALSE)

# Lowercase all the columns
colnames(experimental.structure.table) <- tolower(colnames(experimental.structure.table))

# Load the raw-file-to-condition matrix for label free experiments 
if (is.label.free == TRUE) {
  
  # Path to the raw.files.to.condition.matrix file for label free experiments 
  raw.files.condition.matrix <- paste(here(), "data-input/test-case", "raw-files-to-conditions.csv", sep = "/")
  
  # Read the raw.files.to.condition.matrix file
  label.free.raw.files.condition.matrix <- read.csv(raw.files.condition.matrix,
                                                        stringsAsFactors = FALSE,
                                                        check.names = FALSE)
  
  cat("Label-free experiment: raw.files.to.condition.matrix file loaded!\n")
  
  # Add the file to the global data
  
  global.variables[["raw.files.condition"]] <- label.free.raw.files.condition.matrix
}

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
  
  cat("Reading proteinGroups file: ", protein.groups.file,"...\n", sep = "")
  
  # Fast cleaning of the proteinGroups file using the tr unix command
  cleaning.command <-  paste("tr -d \'\"\\\\\"\' <", protein.groups.file)
  
  # Read the proteinGroups data
  protein.groups.data <- fread(cleaning.command, integer64 = "numeric")
  
  cat("ProteinGroups file loaded!\n")
  
  # Store it in a global variable
  global.variables[["protein.groups.data"]] <- protein.groups.data

  cat("========== End of load_data.R ==========\n")
}


