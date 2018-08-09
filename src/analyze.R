# Here does all the analysis steps. Results are reported in xlsx/csv media format
# in the designated folders.

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

# Make the limma output folder
# make.limma.folder()

# 
# make.Venn.diagram(evidence.data)

### DATA IMPORT STEP ###

# Get the analysis data from the global.variables list
analysis.data <- global.variables[["analysis.data"]]

# Get the conditions to compare vector from the global.variables list
conditions.to.compare <- global.variables[["conditions.to.compare"]]

### FILTERING STEP ###
cat("peos1\n")
# Filter out contaminants and reverse sequences
filtered.data <- filter.out.reverse.and.contaminants(analysis.data)
cat("peos2\n")
### NORMALIZATION STEP ###

# Select the median of the peptide intensities
vsn.normalized.data <- do.vsn.normalization(filtered.data, conditions.to.compare)
cat("peos3\n")
### IMPUTATION ###

imputed.data <- do.LCMD.imputation(vsn.normalized.data)
cat("peos4\n")
### AGGREGATION ### 




cat("========== End of analyze.R ==========\n")

cat(date(),"end \n")

