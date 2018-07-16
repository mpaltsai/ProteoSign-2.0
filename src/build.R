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

# Read the experimental structure from the global variables list
experimental.structure <- global.variables[["experimental.structure"]]

# Order the experimental structure by raw.file name
experimental.structure <- experimental.structure[order(experimental.structure$raw.file),]

# Initialize conditions.to.raw.files list
conditions.to.raw.files.list <- list()

# Get the raw.files.condition.matrix
raw.files.condition.matrix <- global.variables[["raw.files.condition"]]

# Order the raw.files.condition.matrix structure by raw.file name
raw.files.condition.matrix <- raw.files.condition.matrix[order(raw.files.condition.matrix$raw.file),]


# Build the conditions.to.raw.files list from the experimental structure matrix
conditions.to.raw.files.list <- build.condition.to.raw.files.from.matrix(   raw.files.condition.matrix,
                                                                            conditions.to.raw.files.list,
                                                                            is.label.free = TRUE)

# Add conditions.to.raw.files.list to the global variables list
global.variables[["conditions.to.raw.files.list"]] <- conditions.to.raw.files.list

# Read parameters from the global variables list
replicates.multiplexing <- global.variables[["replicate.multiplexing.is.used"]]

experimental.structure$condition <- raw.files.condition.matrix$condition 
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
# but without replicate multiplexing?
biological.replicates.number.status <- check.replicates.number(replicates.multiplexing,
                                                                biological.replicates.list)

# If not inform user and abort analysis
if (biological.replicates.number.status  == FALSE) {
  stop("Cannot accept dataset with just one biological replicate. Aborting...")
}

# Make a list of list where each element is a condition paired
# with its biological and technical replicates
replicates.per.condition <- replicates.per.condition(biological.replicates.list,
                                                     technical.replicates.list,
                                                     experimental.conditions.list)

# Make a list of list where each element is a condition paired
# with booleans regarding the correctness of biological and technical replicates
replicates.status.per.condition <- replicates.status.per.condition(replicates.per.condition)

# Find the problematic replicates and fix them
fixed.replicates.per.condition <- fix.replicates.per.condition(replicates.per.condition,
                                                               replicates.status.per.condition)

# Reset the with the corrected replicates
replicates.per.condition <- fixed.replicates.per.condition

# Make a list with the corrected replicates concatenated
restored.replicates <- restore.replicates(replicates.per.condition)

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

# Update the experimental.structure global variable 
global.variables[["experimental.structure"]] <- experimental.structure

# Store max biological replicates for duplicates handling from limma
global.variables[["max.biological.replicates"]] <- experimental.structure[,
                                                                          which.max(biological.replicate)]

# Store the minimum number of technical replicates
global.variables[["min.technical.replicates"]] <- min(experimental.structure[,
                                                           .SD[which.max(technical.replicate)],
                                                           by = biological.replicate]$technical.replicate)

cat("end of build.R\n")


