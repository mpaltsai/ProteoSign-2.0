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


fixed.replicates.per.condition <- fix.replicates.per.condition(replicates.per.condition,
                                                               replicates.status.per.condition)

replicates.per.condition <- fixed.replicates.per.condition

restored.replicates <- restore.replicates(replicates.per.condition)

biological.replicates.list <- restored.replicates$biological

technical.replicates.list <- restored.replicates$technical

experimental.structure$biorep <- biological.replicates.list

experimental.structure$techrep <- technical.replicates.list
