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

# Get the analysis name
analysis.name <- global.variables$analysis.name

# If the analysis name is bad fix itis bad, fix it 
analysis.name <- gsub("[[:space:]|[:punct:]]", "-", analysis.name)

# Make the limma output folder
make.data.output.folders(analysis.name)

# Get the evidence raw data
evidence.data <- global.variables[["evidence.data"]]

# Make tje Venn diagram
make.Venn.diagram(evidence.data, conditions.to.compare, analysis.name)

# And remove them from the global.variables
global.variables[["evidence.data"]] <- NULL

### DATA IMPORT STEP ###

# Get the analysis data from the global.variables list
analysis.data <- global.variables[["analysis.data"]]

# Get the conditions to compare vector from the global.variables list
conditions.to.compare <- global.variables[["conditions.to.compare"]]

experimental.metadata <- global.variables[["experimental.metadata"]]

# Save the data.table to the intermediate-data
save.intermediate.data.tables(analysis.data, deparse(substitute(analysis.data)), analysis.name)

### FILTERING STEP ###

# Filter out contaminants, reverse sequences and only identified by site
filtered.data <- filter.out.reverse.and.contaminants(analysis.data)

# Save the data.table to the intermediate-data
save.intermediate.data.tables(filtered.data, deparse(substitute(filtered.data)), analysis.name)

### NORMALIZATION STEP ###

# Select the median of the peptide intensities
vsn.normalized.data <- do.vsn.normalization(filtered.data, conditions.to.compare)
not.normalized.data <- do.vsn.normalization(filtered.data,conditions.to.compare, do.norm = FALSE)

# Plot the intensities before and after the normalizations
do.peptide.intensities.plots(not.normalized.data, vsn.normalized.data, analysis.name)

# Save the data.table to the intermediate-data
save.intermediate.data.tables(vsn.normalized.data, deparse(substitute(vsn.normalized.data)), analysis.name)

### IMPUTATION ###

imputed.data <- do.LCMD.imputation(vsn.normalized.data)

# Save the data.table to the intermediate-data
save.intermediate.data.tables(imputed.data, deparse(substitute(imputed.data)))

### AGGREGATION ### 

aggregated.data <- do.peptides.aggregation(imputed.data)

# Save the data.table to the intermediate-data
save.intermediate.data.tables(aggregated.data, deparse(substitute(aggregated.data)), analysis.name)

# Do the QQ plots
do.QQ.plots(aggregated.data, conditions.to.compare, analysis.name)

### DIFFERENTIAL EXPRESSION ###

limma.results <- do.limma.analysis(aggregated.data, conditions.to.compare, experimental.metadata, error.correction.method = "B")

# Save the data.table to the intermediate-data
save.intermediate.data.tables(limma.results, deparse(substitute(limma.results)), analysis.name, output.folder = "limma-output")

### PLOTS ###

# Do the volcano plots
do.volcano.plots(limma.results, conditions.to.compare, analysis.name, plots.format = 5, error.correction.method = "B")

# Do the fold change histograms
do.fold.change.histogram(limma.results, conditions.to.compare, analysis.name)

# Do the MA plots
do.MA.plots(limma.results, conditions.to.compare, analysis.name)

# Do value order plots
do.value.ordered.ratio.plot(limma.results, conditions.to.compare, analysis.name)


# Remove the build.R and functions_build.R from the enviroment
functions.in.analyze.R <- list.functions.in.file("analyze.R")
functions.in.analyze.R <- functions.in.analyze.R$.GlobalEnv

functions.in.functions_analyze.R <- list.functions.in.file("functions_analyze.R")
functions.in.functions_analyze.R <- functions.in.functions_analyze.R$.GlobalEnv

rm(list = c(functions.in.analyze.R, functions.in.functions_analyze.R))

cat("========== End of analyze.R ==========\n")

cat(date(),"end \n")

