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
make.data.output.folders()

# 
# make.Venn.diagram(evidence.data)

### DATA IMPORT STEP ###

# Get the analysis data from the global.variables list
analysis.data <- global.variables[["analysis.data"]]

# Get the conditions to compare vector from the global.variables list
conditions.to.compare <- global.variables[["conditions.to.compare"]]

experimental.metadata <- global.variables[["experimental.metadata"]]

# Save the data.table to the intermediate-data
save.intermediate.data.tables(analysis.data, deparse(substitute(analysis.data)))

### FILTERING STEP ###

# Filter out contaminants and reverse sequences
filtered.data <- filter.out.reverse.and.contaminants(analysis.data)

# Save the data.table to the intermediate-data
save.intermediate.data.tables(filtered.data, deparse(substitute(filtered.data)))

### NORMALIZATION STEP ###

# Select the median of the peptide intensities
vsn.normalized.data <- do.vsn.normalization(filtered.data, conditions.to.compare)
not.normalized.data <- do.vsn.normalization(filtered.data,conditions.to.compare, do.norm = FALSE)

# Plot the intensities before and after the normalizations
do.peptide.intensities.plots(not.normalized.data, vsn.normalized.data)

# Save the data.table to the intermediate-data
save.intermediate.data.tables(vsn.normalized.data, deparse(substitute(vsn.normalized.data)))

### IMPUTATION ###

imputed.data <- do.LCMD.imputation(vsn.normalized.data)

# Save the data.table to the intermediate-data
save.intermediate.data.tables(imputed.data, deparse(substitute(imputed.data)))

### AGGREGATION ### 

aggregated.data <- do.peptides.aggregation(imputed.data)

# Save the data.table to the intermediate-data
save.intermediate.data.tables(aggregated.data, deparse(substitute(aggregated.data)))

# Do the QQ plots
do.QQ.plots(aggregated.data, conditions.to.compare)

### DIFFERENTIAL EXPRESSION ###

limma.results <- do.limma.analysis(aggregated.data, conditions.to.compare, experimental.metadata, error.correction.method = "BH")

# Save the data.table to the intermediate-data
save.intermediate.data.tables(limma.results, deparse(substitute(limma.results)), output.folder = "limma-output")

### PLOTS ###

# Do the volcano plots
do.volcano.plots(limma.results, conditions.to.compare, plots.format = 5)

# Do the MA plots
do.MA.plots(limma.results, conditions.to.compare)

# Do value order plots
do.value.ordered.ratio.plot(limma.results, conditions.to.compare)

cat("========== End of analyze.R ==========\n")

cat(date(),"end \n")

