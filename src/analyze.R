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

### Get the analysis data from the global.variables list
analysis.data <- global.variables[["analysis.data"]]
minimum.peptide.detections <- global.variables[["minimum.peptide.detections"]]
minimum.peptides.per.protein <- global.variables[["minimum.peptides.per.protein"]]
min.valid.values.percentance <- global.variables[["min.valid.values.percentance"]]
knn.neighbors <- global.variables[["knn.neighbors"]]

# Get the experiment metadata
experimental.metadata <- global.variables[["experimental.metadata"]]
is.label.free <- global.variables[["is.label.free"]]

# Get the plots parameters
analysis.name <- global.variables[["analysis.name"]]
plots.format <- global.variables[["plots.format"]]

# Get the limma analysis parameters
error.correction.method <- global.variables[["error.correction.method"]]
fold.change.cut.off <- global.variables[["fold.change.cut.off"]]
FDR <- global.variables[["FDR"]]

# Get the conditions to compare vector from the global.variables list
conditions.to.compare <- global.variables[["conditions.to.compare"]]

# If the analysis name is bad fix itis bad, fix it 
analysis.name <- gsub("[[:space:]|[:punct:]]", "-", analysis.name)

# Make the limma output folder
make.data.output.folders(analysis.name)

# Get the evidence raw data
evidence.data <- global.variables[["evidence.data"]]

# Make the Venn diagram
make.Venn.diagram(evidence.data,
                  conditions.to.compare,
                  analysis.name,
                  is.label.free,
                  minimum.peptide.detections = minimum.peptide.detections,
                  min.valid.values.percentance = min.valid.values.percentance,
                  minimum.peptides.per.protein = minimum.peptides.per.protein,
                  plots.format = plots.format)

### DATA IMPORT STEP ###

# Save the data.table to the intermediate-data
save.intermediate.data.tables(analysis.data, deparse(substitute(analysis.data)), analysis.name)

### FILTERING STEP ###

# Filter out contaminants, reverse sequences and only identified by site
filtered.data <- filter.out.reverse.and.contaminants(analysis.data)

# Save the data.table to the intermediate-data
save.intermediate.data.tables(filtered.data,
                              deparse(substitute(filtered.data)),
                              analysis.name)

### NORMALIZATION STEP ###

# Select the median of the peptide intensities
vsn.normalized.data <- do.vsn.normalization(filtered.data,
                                            conditions.to.compare,
                                            minimum.peptide.detections = minimum.peptide.detections)

not.normalized.data <- do.vsn.normalization(filtered.data,
                                            conditions.to.compare,
                                            minimum.peptide.detections = minimum.peptide.detections,
                                            do.norm = FALSE)

# Plot the intensities before and after the normalizations
do.peptide.intensities.plots(not.normalized.data,
                             vsn.normalized.data,
                             analysis.name,
                             plots.format = plots.format)

# Save the data.table to the intermediate-data
save.intermediate.data.tables(vsn.normalized.data,
                              deparse(substitute(vsn.normalized.data)),
                              analysis.name)

### IMPUTATION ###

imputed.data <- do.knn.imputation(vsn.normalized.data,
                                  conditions.to.compare,
                                  knn.neighbors = knn.neighbors,
                                  min.valid.values.percentance = min.valid.values.percentance)

# Save the data.table to the intermediate-data
save.intermediate.data.tables(imputed.data,
                              deparse(substitute(imputed.data)),
                              analysis.name)

### AGGREGATION ### 

aggregated.data <- do.peptides.aggregation(imputed.data)

# Save the data.table to the intermediate-data
save.intermediate.data.tables(aggregated.data,
                              deparse(substitute(aggregated.data)),
                              analysis.name)

# Do the QQ plots
do.QQ.plots(aggregated.data,
            conditions.to.compare,
            analysis.name,
            experimental.metadata,
            plots.format = plots.format)

### DIFFERENTIAL EXPRESSION ###

limma.results <- do.limma.analysis(aggregated.data,
                                   conditions.to.compare,
                                   experimental.metadata,
                                   error.correction.method = error.correction.method,
                                   fold.change.cut.off = fold.change.cut.off,
                                   FDR = FDR)

# Save the data.table to the intermediate-data
save.intermediate.data.tables(limma.results,
                              deparse(substitute(limma.results)),
                              analysis.name,
                              output.folder = "limma-output")

### FINAL PLOTS ###

# Do the volcano plots
do.volcano.plots(limma.results,
                 conditions.to.compare,
                 analysis.name,
                 plots.format = plots.format,
                 error.correction.method = error.correction.method,
                 fold.change.cut.off = fold.change.cut.off,
                 FDR = FDR)

# Do the fold change histograms
do.fold.change.histogram(limma.results,
                         conditions.to.compare,
                         analysis.name,
                         plots.format = plots.format)

# Do the MA plots
do.MA.plots(limma.results,
            conditions.to.compare,
            analysis.name,
            plots.format = plots.format)

# Do value order plots
do.value.ordered.ratio.plot(limma.results,
                            conditions.to.compare,
                            analysis.name,
                            plots.format = plots.format)

# 
# # Remove the build.R and functions_build.R from the enviroment
# functions.in.analyze.R <- list.functions.in.file("analyze.R")
# functions.in.analyze.R <- functions.in.analyze.R$.GlobalEnv
# 
# functions.in.functions_analyze.R <- list.functions.in.file("functions_analyze.R")
# functions.in.functions_analyze.R <- functions.in.functions_analyze.R$.GlobalEnv
# 
# rm(list = c(functions.in.analyze.R, functions.in.functions_analyze.R))

cat("========== End of analyze.R ==========\n")

cat(date(),"end \n")

