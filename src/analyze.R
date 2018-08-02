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

# Get the analysis data from the global.variables list
analysis.data <- global.variables[["analysis.data"]]

# Get the conditions to compare vector from the global.variables list
conditions.to.compare <- global.variables[["conditions.to.compare"]]
  
# Filter out contaminants and reverse sequences
analysis.data <- filter.out.reverse.and.contaminants(analysis.data)

# Select the median of the peptide intensities
median.intensities <- use.peptides.median(analysis.data, conditions.to.compare)

# Replace the multiple peptide intensities with the median from the 2 conditions
analysis.data[, eval(conditions.to.compare[1]) := median.intensities[[1]]$`Median Intensity`]

analysis.data[, eval(conditions.to.compare[2]) := median.intensities[[2]]$`Median Intensity`]

# And rename the columns as "Condition X Median Intensity"
setnames(analysis.data,
         conditions.to.compare[1],
         paste(conditions.to.compare[1],
               "Median Intensity"))

setnames(analysis.data,
         conditions.to.compare[2],
         paste(conditions.to.compare[2],
               "Median Intensity"))

cat("========== End of analyze.R ==========\n")

cat(date(),"end \n")

