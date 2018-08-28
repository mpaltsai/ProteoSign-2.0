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
global.variables[["analysis.name"]] <- analysis.metadata$analysis.name
global.variables[["replicate.multiplexing.is.used"]] <- analysis.metadata$replicate.multiplexing.is.used
global.variables[["dataset.origin"]] <- analysis.metadata$dataset.origin
global.variables[["is.label.free"]]  <- analysis.metadata$is.label.free
global.variables[["is.isobaric"]]    <- analysis.metadata$is.isobaric
global.variables[["timestamp.to.keep"]] <- analysis.metadata$timestamp.to.keep
global.variables[["culture.to.keep"]] <- analysis.metadata$culture.to.keep
global.variables[["raw.files.to.remove"]] <- analysis.metadata$raw.files.to.remove
global.variables[["raw.files.to.rename"]] <- analysis.metadata$raw.files.to.rename
global.variables[["conditions.to.compare"]] <- unlist(strsplit(analysis.metadata$conditions.to.compare, split = ","))

# Get the dataset origin
dataset.origin <- global.variables[["dataset.origin"]]

# Which software do the data come from
if (dataset.origin == "MaxQuant") {
  cat("Data origin: MaxQuant.\n")
} else {
  cat("Data origin: Proteome Discoverer.\n")
}

is.label.free <- global.variables[["is.label.free"]]

# Path to the experimental structure file
experimental.structrure.file <- paste(here(), "data-input/test-case", "experimental-structure.csv", sep = "/")

# Path to the experimental structure file
evidence.file <- paste(here(), "data-input/test-case", "evidence.txt", sep = "/")

cat("Reading experimental structure file...\n")

# Read the experimental structure file
experimental.structure.table <- read.csv(experimental.structrure.file,
                                         stringsAsFactors = FALSE,
                                         check.names = FALSE)

# Lowercase all the columns of the experimental structure file
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

# Read the evidence data but only for the needed columns
evidence.data <- fread(cleaning.command,
                       integer64 = "numeric")

# Get the evidence data column names
evidence.column.names <- colnames(evidence.data)

# Trim and lowercase the column names for universal handling across
# MaxQuant/Proteome Discoverer versions 
trimmed.and.lowercased.column.names <- trim.and.lowercase.column.names(evidence.column.names)

# Set the trimmed and lowercase the evidence columns
colnames(evidence.data) <- trimmed.and.lowercased.column.names

evidence.data.column.names <- colnames(evidence.data)
# The needed columns for the analysis depending on the software
if (dataset.origin == "MaxQuant") {
  
  # Get the "intensity.X" columns where X stands for l, m, h respectivelly
  intensity.columns <- grep("^intensity\\.", evidence.data.column.names,perl = TRUE, value = TRUE)
  
  # Get the conditions to compare
  conditions.to.compare <- global.variables$conditions.to.compare
  
  # Make a pattern e.g "(h|m)$"
  conditions.to.compare.pattern <- paste0("(", tolower(conditions.to.compare), ")$", collapse = "|")
  
  # And keep only the intensity columns I want to compare
  intensity.columns <- grep(conditions.to.compare.pattern, intensity.columns, perl = TRUE, value = TRUE)
  
  evindence.columns.to.keep <- c("proteins",
                                 "raw.file",
                                 "protein.ids",
                                 "protein.names",
                                 "id",
                                 "protein.descriptions",
                                 "peptide.id",
                                 "unique.sequence.id",
                                 intensity.columns)
  
  # In any case, we take the intersection
  evindence.columns.subset <- intersect(evidence.data.column.names,
                                             evindence.columns.to.keep)
  
  # Now subset the columns to keep only the needed, in order to make
  # the data.table as light-weight as possible
  evidence.data <- evidence.data[, .SD, .SDcols = evindence.columns.subset]
  
  
} else {
  # Case proteome discoverer
  evindence.columns.to.keep <- c("unique.sequence.id",
                                 "annotated.sequence")
}

cat("Evidence file loaded!\n")

# Get the timestamp and the culture to analyze
timestamp.to.keep <- global.variables[["timestamp.to.keep"]]
culture.to.keep <- global.variables[["culture.to.keep"]]

# Clean evidence file on specific timestamp or culture
evidence.data <- keep.only.specific.timestamps.or.cultures(evidence.data,
                                                           timestamp.to.keep,
                                                           culture.to.keep)

# If there ar raw file to be renamed, rename them
if (is.na(global.variables[["raw.files.to.remove"]]) == FALSE & 
    is.na(global.variables[["raw.files.to.rename"]]) == FALSE) {
  
  raw.files.to.remove <- unlist(strsplit(global.variables[["raw.files.to.remove"]], split = ","))
  raw.files.to.rename <- unlist(strsplit(global.variables[["raw.files.to.rename"]], split = ","))
  
  evidence.data <- remove.and.rename.raw.files(evidence.data, raw.files.to.remove, raw.files.to.rename)
  
}

# Add the evidence data to the global variables list
global.variables[["evidence.data"]] <- evidence.data

# If data come from MaxQuant Software, find the proteinGroups file and read it
if (dataset.origin == "MaxQuant") {
  
  # Path to the proteinGroups  file 
  protein.groups.file <- paste(here(), "data-input/test-case", "proteinGroups.txt", sep = "/")
  
  cat("Reading proteinGroups file: ", protein.groups.file,"...\n", sep = "")
  
  # Fast cleaning of the proteinGroups file using the tr unix command
  cleaning.command <-  paste("tr -d \'\"\\\\\"\' <", protein.groups.file)
  
  # Read the proteinGroups data
  protein.groups.data <- fread(cleaning.command,
                               integer64 = "numeric")
  
  # Get the protein.groups.data column names
  protein.column.names <- colnames(protein.groups.data)
  
  # Trim and lowercase the column names for universal handling across
  # MaxQuant/Proteome Discoverer versions 
  trimmed.and.lowercased.column.names <- trim.and.lowercase.column.names(protein.column.names)
  
  # Set the trimmed and lowercase the evidence columns
  colnames(protein.groups.data) <- trimmed.and.lowercased.column.names
  
  # Ideally the protein groups file should contain these files
  protein.groups.columns.to.keep <- c("protein.ids",
                                      "fasta.headers",
                                      "protein.names",
                                      "only.identified.by.site",
                                      "reverse",
                                      "contaminant",
                                      "id",
                                      "peptide.ids",
                                      "evidence.ids")
  
  # In any case, we take the intersection
  protein.groups.columns.subset <- intersect(colnames(protein.groups.data),
                                             protein.groups.columns.to.keep)
  
  # Now subset the columns to keep only the needed, in order to make
  # the data.table as light-weight as possible
  protein.groups.data <- protein.groups.data[, .SD, .SDcols = protein.groups.columns.subset]
  
  cat("ProteinGroups file loaded!\n")
  
  # Store it in a global variable
  global.variables[["protein.groups.data"]] <- protein.groups.data

}

cat("========== End of load_data.R ==========\n")
# 
# # Remove the load_data.R and functions_load_data.R functions from the enviroment
# functions.in.load_data.R <- list.functions.in.file("load_data.R")
# functions.in.load_data.R <- functions.in.load_data.R$.GlobalEnv
# 
# 
# functions.in.functions_load_data.R <- list.functions.in.file("functions_load_data.R")
# functions.in.functions_load_data.R <- functions.in.functions_load_data.R$.GlobalEnv
# 
# rm(list = c(functions.in.load_data.R, functions.in.functions_load_data.R))
# 
