# Data wrangling with dplyr/tidyr/etc. All the magic happens here. Newly created
# files in workspace should be displayed. A build boolean variable can be used for
# data reload from RDS for faster data reload.

# Clear enviroment and only keep functions and global/project variables
rm(list = grep(paste(c("^global.variables",
                       "^project.variables",
                       lsf.str()),
                     collapse = "|"),
               ls(),
               value = TRUE,
               invert = TRUE))

# Return the memory to the OS
gc(verbose = FALSE,
   reset = TRUE)

### Load all the needed global variables ###

# Read the experimental structure from the global variables list
experimental.structure <- global.variables[["experimental.structure"]]

# Read the analysis arguments
conditions.to.compare <- global.variables[["conditions.to.compare"]]
is.label.free <- global.variables[["is.label.free"]]
is.isobaric <- global.variables[["is.isobaric"]]
dataset.origin <- global.variables[["dataset.origin"]]
replicates.multiplexing <- global.variables[["replicate.multiplexing.is.used"]]
evidence.data <- global.variables[["evidence.data"]]
protein.groups.data <- global.variables[["protein.groups.data"]]


# Order the experimental structure by raw.file name
experimental.structure <- experimental.structure[order(experimental.structure$raw.file),]

# Be sure that we have 2 conditions to compare
if (length(conditions.to.compare) != 2) {
  stop("Invalid conditions.to.compare arguments. Conditions has to be exactly 2.\n")
}


if (is.label.free == TRUE) {
  # Initialize conditions.to.raw.files list
  conditions.to.raw.files.list <- list()
  
  # Get the raw.files.condition.matrix in case of label-free experiments
  raw.files.condition.matrix <- global.variables[["raw.files.condition"]]
  
  # Order the raw.files.condition.matrix structure by raw.file name
  raw.files.condition.matrix <- raw.files.condition.matrix[order(raw.files.condition.matrix$raw.file),]
  
  # Build the conditions.to.raw.files list from the experimental structure matrix
  conditions.to.raw.files.list <- build.condition.to.raw.files.from.matrix(   raw.files.condition.matrix,
                                                                              is.label.free = TRUE)
  
  # # The conditions that we want to focus on
  # conditions.to.compare <- names(conditions.to.raw.files.list[conditions.to.compare])
  
  # Add conditions.to.raw.files.list to the global variables list
  global.variables[["conditions.to.raw.files.list"]] <- conditions.to.raw.files.list
  
  print("Experiment type: Label-free Experiment")
} 

# In case of isobaric experiment, fix the evidence columns
if (is.isobaric == TRUE) {
  
  # The the tags to conditions
  tags.to.conditions <- global.variables[["tags.to.conditions"]]
  
  # And combine the "reporter.intensity.X" columns into one for each condition and convert 0s to NAs
  evidence.data <- merge.reporter.intensity.columns(evidence.data, tags.to.conditions)
  
  print("Experiment type: Isobaric Experiment (TMT or iTRAQ)")
}

# If we are on a label-free experiment, add the conditions, otherwise just add a generic description
if (is.label.free == TRUE) {
  experimental.structure$condition <- raw.files.condition.matrix$condition 
} else {
  experimental.structure$condition <- "Labeled Experiment"
}

# Reorder the experimental structure based on conditions/biological replicates
# /technical replicates/ fractions
experimental.structure <- experimental.structure[
  order(experimental.structure$condition,
        experimental.structure$biological.replicate,
        experimental.structure$technical.replicate,
        experimental.structure$fraction), ]

# Correct the rownames
rownames(experimental.structure) <- c(1:length(experimental.structure$raw.file))

# Store each column on a separate variable
biological.replicates.list <- experimental.structure$biological.replicate
technical.replicates.list <- experimental.structure$technical.replicate
fraction.list <- experimental.structure$fraction
experimental.conditions.list <- experimental.structure$condition

# Do I have more than 1 replicate 
# but without replicate multiplexing? COMBAK 
biological.replicates.number.status <- check.replicates.number(replicates.multiplexing,
                                                                biological.replicates.list)

# If not inform user and abort analysis
if (biological.replicates.number.status  == FALSE) {
  stop("Cannot accept dataset with just one biological replicate. Aborting...")
}

# Make a list of list where each element is a condition paired
# with its biological and technical replicates
replicates.per.condition.list <- replicates.per.condition(biological.replicates.list,
                                                     technical.replicates.list,
                                                     experimental.conditions.list)

# Make a list of list where each element is a condition paired
# with booleans regarding the correctness of biological and technical replicates
replicates.status.per.condition.list <- replicates.status.per.condition(replicates.per.condition.list)

# Find the problematic replicates and fix them
fixed.replicates.per.condition <- fix.replicates.per.condition(replicates.per.condition.list,
                                                               replicates.status.per.condition.list,
                                                               is.label.free)

# Reset the with the corrected replicates
replicates.per.condition.list <- fixed.replicates.per.condition

# Make a list with the corrected replicates concatenated
restored.replicates <- restore.replicates(replicates.per.condition.list)

# Are the biological and the technical replicates different from the initial data, provided by the user?
biological.replicates.are.the.same <- Reduce("&", biological.replicates.list == restored.replicates$biological.replicates)
technical.replicates.are.the.same <- Reduce("&", technical.replicates.list == restored.replicates$technical.replicates)

# If they are different, inform the user!
if (biological.replicates.are.the.same == FALSE |
    technical.replicates.are.the.same == FALSE) {
    message("Attetion: The experimental structure has been corrected due to mistyped replicates (gaps in the numbering e.g. 1, 2, 4 instead of 1, 2, 3).")  
}

# Restore the biological replicates list
biological.replicates.list <- restored.replicates$biological.replicates

# Restore the technical replicates list
technical.replicates.list <- restored.replicates$technical.replicates

# Restore the biological column on the experimenta structure file
experimental.structure$biological.replicate <- biological.replicates.list

# Restore the technical column on the experimenta structure file
experimental.structure$technical.replicate <- technical.replicates.list

# Now exploit the check.number.of.replicates to find if there are
# replicates of each type in our experiment and construct the setup id
# ID = 1: We have only biological replicates 
# ID = 2: Illegal state! We cannot have only technical replicates...
# ID = 3: We have biological and technical replicates
# ID = 4: We have biological replicates and fractions
# ID = 5: Illegal state! There are no biological replicates...
# ID = 6: We have biological and technical replicates as well as fractions
# TODO ensure that this works with replicate multiplexing too,
# and if yes  then maybe the ID 5 is then correct?
experimental.setup.id <-  1 * check.replicates.number(replicates.multiplexing,
                                                         biological.replicates.list) +
                          
                          2 * check.replicates.number(replicates.multiplexing,
                                                         technical.replicates.list) +
                          
                          3 * check.replicates.number(replicates.multiplexing,
                                                         fraction.list)

experimental.description <- make.experimental.description(experimental.setup.id,
                                                          biological.replicates.list,
                                                          technical.replicates.list,
                                                          fraction.list)

# Add a description column to the experimental structure matrix
experimental.structure$description <- experimental.description

# Turn experimental structure into a data.table
experimental.structure <- data.table(experimental.structure)

# Add the number of biological/ technical replicates
# and fraction per condition as metadata in a variable
experimental.metadata <- get.experiment.metadata(experimental.structure) 

# And put them as well in a global variable
global.variables[["experimental.metadata"]] <- experimental.metadata

# Update the experimental.structure global variable 
global.variables[["experimental.structure"]] <- experimental.structure

# If we are not in a isobaric experiment, we replace 0s with NAs here, otherwise, we do it earlier as it is easier
if (is.isobaric == FALSE) {
  # Replace the 0s with NAs
  evidence.data <- zeros.to.nas(evidence.data, prefix = "^intensity")
}

# Now build the analysis data 
analysis.data <- build.analysis.data(protein.groups.data = protein.groups.data,
                                     evidence.data       = evidence.data,
                                     dataset.origin      = dataset.origin,
                                     is.label.free       = is.label.free,
                                     is.isobaric         = is.isobaric)

# Store the data in a global variable
global.variables[["analysis.data"]] <- analysis.data

cat("after analysis data\n")

# Now depending on the condition column type we will do different operations on the evidence.data
experiment.type <- unique(analysis.data$condition)

# If we are on an Isotopic Experiment
if (length(experiment.type) == 1) {
  
  # Keep the evidence for the venn diagram
  global.variables[["evidence.data"]] <- copy(analysis.data)
  
  # Find the Intensity column names
  intensity.columns <- grep("^intensity", colnames(analysis.data), perl = TRUE, value = TRUE)
  
  # Find the rows on those columns that contain NA
  na.rows.conditon.A <- which(is.na(analysis.data[[intensity.columns[1]]]) == TRUE)
  na.rows.conditon.B <- which(is.na(analysis.data[[intensity.columns[2]]]) == TRUE)
  
  # Take their union
  na.rows <- union(na.rows.conditon.A, na.rows.conditon.B)
  
  # Remove NA rows only if they exist!
  if (length(na.rows) != 0) {
    # And remove those rows in order to have rows where we have measurements in both conditions
    analysis.data <- analysis.data[-c(na.rows), ]
  }
  
  
  # Discard the condition columns as now there is no need for it
  analysis.data <- analysis.data[, condition := NULL]
  
  # Find the position of the intensity columns
  intensity.columns.position <- grep("^intensity", colnames(analysis.data), perl = TRUE)
  
  # Discard the "intensity." part of each "intensity.x"
  intensity.columns <- gsub(".*\\.", "", intensity.columns, perl = TRUE)
  
  # Uppercase the condition column names
  intensity.columns <- toupper(intensity.columns)
  
  # And restore the condition column names
  colnames(analysis.data)[intensity.columns.position] <- intensity.columns
  
} else {
  
  # Cast the multiline conditions to 2 columns Condition1 Condition2 with their peptides' intensities
  analysis.data.mixed <- dcast.data.table(analysis.data,
                                          description + protein.ids + unique.sequence.id ~ condition,
                                          value.var = "intensities")
  
  # Score with 1 or 0 (True or False) if we have at least one intensity
  # for each peptide
  peptides.present.in.analysis.data <- dcast.data.table(analysis.data,
                                                        description + protein.ids + unique.sequence.id ~ condition,
                                                        value.var = "intensities",
                                                        fun.aggregate = length)
  
  # Find the peptides detected in condition A
  peptides.present.in.condition.A <- which(peptides.present.in.analysis.data[[experiment.type[1]]] == 1)
  
  # Find the peptides detected in condition B
  peptides.present.in.condition.B <- which(peptides.present.in.analysis.data[[experiment.type[2]]] == 1)
  
  # Now get the intersection
  peptides.present.in.both.conditions <- intersect(peptides.present.in.condition.A, 
                                                   peptides.present.in.condition.B)
  
  # And subset the analysis.data to get the peptides that exist in both conditions
  analysis.data <- analysis.data.mixed[peptides.present.in.both.conditions, ]
  
  # Remove NA proteins
  analysis.data <- na.omit(analysis.data, cols = "protein.ids")
  
  global.variables[["evidence.data"]] <- copy(analysis.data)
  
}

# Split multi preotein ids
analysis.data <- cSplit(analysis.data, "protein.ids", sep=";", direction = "long")

# Remove descriptions from protein ids columns
analysis.data$protein.ids <- gsub(" \\[.*", "",
                                  analysis.data$protein.ids,
                                  perl = TRUE,
                                  fixed = FALSE)


# Store the data in a global variable
global.variables[["analysis.data"]] <- analysis.data

# # Remove the protein groups data as it is not needed any more
# global.variables[["protein.groups.data"]] <- NULL

cat("========== End of build.R ==========\n")

# 
# # Remove the build.R and functions_build.R from the enviroment
# functions.in.build.R <- list.functions.in.file("build.R")
# functions.in.build.R <- functions.in.build.R$.GlobalEnv
# 
# 
# functions.in.functions_build.R <- list.functions.in.file("functions_build.R")
# functions.in.functions_build.R <- functions.in.functions_build.R$.GlobalEnv
# rm(list = c(functions.in.build.R, functions.in.functions_build.R))
# 
# 
