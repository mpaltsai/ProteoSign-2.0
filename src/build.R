# Data wrangling with dplyr/tidyr/etc. All the magic happens here. Newly created
# files in workspace should be displayed. A build boolean variable can be used for
# data reload from RDS for faster data reload.

# Read the experimental structure from the global variables list
experimental.structure <- global.variables[["experimental.structure"]]

# Read parameters from the global variables list
replicates.multiplexing <- global.variables[["replicate.multiplexing.is.used"]]

# Reorder the experimental structure based on conditions/biological replicates
# /technical replicates/ fractions
experimental.structure <- experimental.structure[
  order(experimental.structure$conditions,
        experimental.structure$biorep,
        experimental.structure$techrep,
        experimental.structure$fractions),]

# Correct the rownames
rownames(experimental.structure) <- c(1:length(experimental.structure$raw.file))

# Store each column on a separate variable
biological.replicates.list <- experimental.structure$biorep
technical.replicates.list <- experimental.structure$techrep
experimental.fraction.list <- experimental.structure$fraction
experimental.conditions.list <- experimental.structure$condition

# Do I have more than 1 replicate 
# but without replicate multiplexing?
biological.replicates.number.status <- check.replicates.number(replicates.multiplexing,
                                                                biological.replicates.list)

# If not inform user and abort analysis
if (biological.replicates.number.status  == FALSE) {
  cat("Error User: Cannot accept dataset with just one biological replicate. Aborting ...\n")
  # TODO Handle the sourcing in order to stop
  return (FALSE)
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

# Restore the biological replicates list
biological.replicates.list <- restored.replicates$biological

# Restore the technical replicates list
technical.replicates.list <- restored.replicates$technical

# Restore the biological column on the experimenta structure file
experimental.structure$biorep <- biological.replicates.list

# Restore the technical column on the experimenta structure file
experimental.structure$techrep <- technical.replicates.list

# TODO If bioreps and techreps are changed inform the user!

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
experimental.setup.id <-  1 * check.number.of.replicates(replicates.multiplexing,
                                                         biological.replicates.list) +
                          
                          2 * check.number.of.replicates(replicates.multiplexing,
                                                         technical.replicates.list) +
                          
                          3 * check.number.of.replicates(replicates.multiplexing,
                                                         experimental.fraction.list)

experimental.description <- make.experimental.description(experimental.setup.id,
                                                          biological.replicates.list,
                                                          technical.replicates.list,
                                                          experimental.fraction.list)

# Add a description column to the experimental structure matrix
experimental.structure$description <- experimental.description

# Update the global variable and turn in to data.table
experimental.structure <- data.table(experimental.structure)
global.variables[["experimental.structure"]] <- experimental.structure

# Store max biological replicates for duplicates handling from limma
global.variables[["max.biological.replicates"]] <- experimental.structure[,
                                                                          which.max(biological)]

# Turn experimental structure to data.table
global.variables[["min.technical.replicates"]] <- min(data[,
                                                           .SD[which.max(technical)],
                                                           by=biological]$technical)