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

# Initialize the global variables with the analysis parameters
global.variables <- add.analysis.parameters.to.global.variables(analysis.metadata)

# Get the dataset origin
dataset.origin <- global.variables[["dataset.origin"]]

# Which software do the data come from
if (dataset.origin == "MaxQuant") {
  cat("Data origin: MaxQuant.\n")
} else {
  cat("Data origin: Proteome Discoverer.\n")
}

# Get the experiment type
is.label.free <- global.variables[["is.label.free"]]
is.isobaric <- global.variables[["is.isobaric"]]

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

# Is the experiment isobaric
if (is.isobaric == TRUE) {
  
  # Path to the tags.to.conditions file for isotopic experiments 
  tags.to.conditions.path <- paste(here(),
                                   "data-input/test-case",
                                   "tags-to-conditions.csv",
                                   sep = "/")
  
  # Read the tags.to.conditions file
  tags.to.conditions.matrix <- read.csv(tags.to.conditions.path,
                                        stringsAsFactors = FALSE,
                                        check.names = FALSE)
  
  cat("Isotopic experiment: tags-to-conditions file loaded!\n")
  
  # Add the "reporter.intensity." prefix to the tag column
  tags.to.conditions.matrix$tag <- paste0("reporter.intensity.",
                                          tags.to.conditions.matrix$tag)
  
  # Match columns to the conditions
  tags.to.conditions <- make.tags.to.conditions.list(tags.to.conditions.matrix)
  
  # Add the list to the global data
  global.variables[["tags.to.conditions"]] <- tags.to.conditions
  
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

# Get the evidence data column names
evidence.data.column.names <- colnames(evidence.data)

# 1:  For Isotopic experiment e.g. SILAC
# 2:  For Label-free
# 3:  For Isobaric e.g. TMT
# Make an experimental setup code for fast case switch
experimental.type <-  (is.label.free * 1) +
                      (is.isobaric * 2)   + 1

# The needed columns for the analysis depending on the software
if (dataset.origin == "MaxQuant") {
  switch( experimental.type,
         {
           # Case Isotopic
           intensity.columns <- grep("^intensity",
                                     evidence.data.column.names,
                                     perl = TRUE,
                                     value = TRUE)
           
            # Get the conditions to compare parameter value
            conditions.to.compare <- global.variables[["conditions.to.compare"]]
            
            # Lowercase them
            conditions.to.compare <- tolower(conditions.to.compare)
            
            # Make a pattern e.g "\\.(h|m)$"
            conditions.to.compare.pattern <- paste0("\\.(",
                                                    conditions.to.compare[1],
                                                    "|",
                                                    conditions.to.compare[2],
                                                    ")$")
            
            # Get the "intensity.X" columns where X stands for l, m, h respectivelly from the
            # intensity columns
            intensity.columns <- grep(conditions.to.compare.pattern,
                                      intensity.columns,
                                      perl = TRUE,
                                      value = TRUE)
          },
         {
           # Case Label-Free
           intensity.columns <- grep("^intensity",
                                     evidence.data.column.names,
                                     perl = TRUE,
                                     value = TRUE) 
         },
         {
           # Case Isobaric
           intensity.columns <- grep("^reporter\\.intensity\\.[[:digit:]]+$",
                                     evidence.data.column.names,
                                     perl = TRUE,
                                     value = TRUE)
         })
  
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
subset.to.keep <- global.variables[["subset.to.keep"]]

# If we have to select only a subset based on a timestamp or a culture
if (timestamp.to.keep != "" | 
    subset.to.keep != "") {

  # Clean evidence file on specific timestamp or culture
  evidence.data <- keep.only.specific.timestamps.or.cultures(evidence.data,
                                                             timestamp.to.keep,
                                                             subset.to.keep)
}

# If there ar raw file to be renamed, rename them
if (global.variables[["raw.files.to.remove"]] != "" & 
    global.variables[["raw.files.to.rename"]] != "") {

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
  protein.groups.data <- protein.groups.data[, .SD,
                                             .SDcols = protein.groups.columns.subset]
  
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
