make.folder <- function(folder.name, analysis.name) {
  #
  # Makes a folder tree for the analysis, the root the analysis names and leafs
  # the "limma-output", "plots", "intermediate-data" folders
  #
  # Args:
  #   folder.name:    The folder name to create
  #   analysis.name:  The analysis folder name
  #
  # Returns:
  #   The newly created folders
  #
  
  # Get the data-output path
  data.output.path <- here("data-output")
  
  # Make the full path
  folder.output <- paste(data.output.path, analysis.name, folder.name, sep = "/")
  
  # If folder already exist, remove it
  if( file.exists(folder.output) == TRUE) {
    unlink(folder.output, recursive=T, force=T)
  }
  
  # And then create the new folder
  tryCatch(
    {
      dir.create(folder.output, recursive = TRUE)
      success.message <- paste0("Folder ", folder.output, " was created under ", data.output.path, "/!\n")
      cat(success.message)
    },
    error = function(cond) {
      stop(paste0("Error in make.",folder.output," folder: ", cond$message))
    },
    warning = function(cond) {
      stop(paste0("Warning in make.",folder.output," folder: ", cond$message))
    })
}

make.data.output.folders <- function(analysis.name) {
  #
  # Makes the data output folders limma-output, plots 
  #
  # Args:
  #   analysis.name: The name of the experiment e.g. "SILAC HUMAN"
  #
  # Returns:
  #   Prints a message if the folder was created successfully, otherwise terminates the analysis.
  #
  
  # Set limma output path
  data.output.path <- here("data-output")
  
  
  # Set limma output folder name
  folders.to.create <- c("limma-output", "plots", "intermediate-data")
  
  invisible(lapply(folders.to.create, make.folder, analysis.name))
  
}

make.Venn.diagram <- function(evidence.data, conditions.to.compare, analysis.name, is.label.free,
                              minimum.peptide.detections = 2, plots.format = 5) {
  #
  # Makes a Venn diagram between 2 conditions
  #
  # Args:p
  #   evidence.data:          The data.table with the evidence data
  #   conditions.to.compare:  A vector with the 2 conditions to compare
  #   analysis.name:          The analysis title provided by the user
  #   is.label.free:          The experiment type, is it label-free or not        
  #
  #   minimum.peptide.detections: Default is 2. A strictness parameter regarding how many times is a peptide measured
  #                               If it is measures less than N times, the peptide is considered as detected but not measures
  #                               hence its value is set to NA
  #
  #   plots.format:           Default is 5 (jpeg format). A numeric value indicating the format of the plots. 
  #                           The numbers correspont to:  1     2     3     4      5      6     7     8     9
  #                                                     "eps" "ps"  "tex" "pdf" "jpeg" "tiff" "png" "bmp" "svg"
  #
  # Returns:
  #   A venn diagram in png format and the data needed for its construction
  #
  
  # Get the evidence data
  data <- copy(evidence.data)
  
  # Remove the fraction information
  data <- data[, description := gsub("F.*$","", description, perl = TRUE)]
  
  if (is.label.free == TRUE) {
    
    # Make the condition pattern e.g. ".*.h$", ".*.wild$", ".*.mutant$" etc for each condition 
    condition.A.column <- conditions.to.compare[1]
    condition.B.column <- conditions.to.compare[2]
    
    # Rename the column
    setnames(data, "protein.ids", "proteins")
    
  } else {
    # Make the condition pattern e.g. ".*.h$" etc for each condition 
    condition.A.column.pattern <- paste0(".*\\.", tolower(conditions.to.compare[1]),"$")
    condition.B.column.pattern <- paste0(".*\\.", tolower(conditions.to.compare[2]),"$")
    
    # Find those columns
    condition.A.column <- grep(condition.A.column.pattern, colnames(data), perl = TRUE, value = TRUE) 
    condition.B.column <- grep(condition.B.column.pattern, colnames(data), perl = TRUE, value = TRUE) 
    
    # Keep all the columns except the opposite condition
    condition.A <- data[, .SD, .SDcols = !c(condition.B.column, "condition")]
    condition.B <- data[, .SD, .SDcols = !c(condition.A.column, "condition")]
  }
  
  # Filter out contaminants
  condition.A <- filter.out.reverse.and.contaminants(condition.A)
  condition.B <- filter.out.reverse.and.contaminants(condition.B)
  
  # Transfrorm the data from long format to wide where we have the peptides on the y axis and the biological
  # replicates on the x axis
  condition.A <- dcast(condition.A,
                     protein.ids + unique.sequence.id ~ description,
                     value.var = condition.A.column,
                     fun.aggregate = use.peptides.median)
  
  condition.B <- dcast(condition.B,
                       protein.ids + unique.sequence.id ~ description,
                       value.var = condition.B.column,
                       fun.aggregate = use.peptides.median)
  
  # Get the number of columns on each condition
  samples.column.number.A <- ncol(condition.A)
  samples.column.number.B <- ncol(condition.B)
  
  # Now for each row compute the NA percentance
  condition.A[, na.percentance := ((samples.column.number.A - 2 - Reduce("+", lapply(.SD, is.na))) * 100 )/samples.column.number.A,
               .SDcols = c(3:samples.column.number.A)]
  
  condition.B[, na.percentance := ((samples.column.number.B - 2 - Reduce("+", lapply(.SD, is.na))) * 100 )/samples.column.number.B,
              .SDcols = c(3:samples.column.number.B)]
  
  # Then find the rows with correct observations
  rows.with.good.observations.A <- which(condition.A$na.percentance > max.na.percentance.per.row)
  rows.with.good.observations.B <- which(condition.B$na.percentance > max.na.percentance.per.row)
  
  # Keep only them
  condition.A <- condition.A[rows.with.good.observations.A, ]
  condition.B <- condition.B[rows.with.good.observations.B, ]
  
  # Remove the useless na.percentance column
  condition.A[, na.percentance := NULL]
  condition.B[, na.percentance := NULL]
  
  # Then remove the intensities columns
  condition.A[, number.of.peptides := .N, by = "protein.ids"]
  condition.B[, number.of.peptides := .N, by = "protein.ids"]
  
  # Find the peptides with at least 1 intensity
  peptides.above.threshold.A <- which(condition.A$number.of.peptides >= minimum.peptide.detections)
  peptides.above.threshold.B <- which(condition.B$number.of.peptides >= minimum.peptide.detections)
  
  # Keep the ones with at least 1 intensity
  condition.A <- condition.A[peptides.above.threshold.A,]
  condition.B <- condition.B[peptides.above.threshold.B,]
  
  # Get the protein names on this sample
  condition.A.proteins <- unique(condition.A$protein.ids)
  condition.B.proteins <- unique(condition.B$protein.ids)
  
  # TODO
  
  # Find which vector has the maximum length
  max.proteins <- max(length(condition.A.proteins),
                      length(condition.B.proteins))
  
  # So if the 2 vectors have different legnths, fill the shortest of the 2 with NAs
  length(condition.A.proteins) <- max.proteins
  length(condition.B.proteins) <- max.proteins
  
  # Put them in a matrix
  venn.table <- cbind(condition.A.proteins, condition.B.proteins)
  
  # Correct the column names
  colnames(venn.table) <- conditions.to.compare
  
  # Wrap them to a data.table
  venn.table <- as.data.table(venn.table)
  
  # Save them in a csv
  save.intermediate.data.tables(venn.table, "venn-data", analysis.name)
  
  # Remove the previously added NAs
  condition.A.proteins <- condition.A.proteins[!is.na(condition.A.proteins)]
  condition.B.proteins <- condition.B.proteins[!is.na(condition.B.proteins)]
  
  # Make the set
  sets <- list(condition.A.proteins, condition.B.proteins)
  
  # Give the names of the conditions to the list's set
  names(sets) <- conditions.to.compare
  
  # By default the VennDiagram package makes a log file, that we want to supperss
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  
  # Initialize the plot format
  plot.format <- ""
  
  # In user has specified the he preffed different image format, change to the appropriate format
  switch(plots.format,
         plot.format <- "eps",
         plot.format <- "ps",
         plot.format <- "tex",
         plot.format <- "pdf",
         plot.format <- "jpeg",
         plot.format <- "tiff",
         plot.format <- "png",
         plot.format <- "bmp",
         plot.format <- "svg")
  
  # Prepare the plots path
  plots.path <- paste("data-output", analysis.name, "plots",sep = "/")
  
  # Change to the plots path
  setwd(here(plots.path))
  
  # Now draw the venn diagram
  venn.diagram <- venn.diagram(sets,
                               filename = NULL,
                               main = paste0("Venn diagram between the conditions ",
                                             conditions.to.compare[1],
                                             " and ",
                                             conditions.to.compare[2]),
                               category = conditions.to.compare,
                               scaled = FALSE,
                               alpha = c(0.7, 0.7),
                               fontfamily = "Helvetica",
                               main.fontface = "bold",
                               fill = c("#999999", "#E69F00"),
                               cat.default.pos = "text",
                               cat.pos = c(1,1),
                               cex =  2,
                               cat.cex = 2,
                               cat.dist = c(0.05, 0.05))
  
  # Save the boxplot after normalization with no description in high resolution 1920x1080 pixels, in the appropriate format
  ggsave(filename = paste0("venn-diagram",".",plot.format),
         plot = venn.diagram,
         device = plot.format,
         width = 6.4,
         height = 3.6,
         units = "in",
         dpi = "print",
         limitsize = FALSE)
  
  # And finally reset the working directory to the scripts folder
  setwd(here("src")) 
}

save.intermediate.data.tables <- function(data, file.name, analysis.name,
                                          output.folder = "intermediate-data") {
  #
  # Saves the intermediate data.table into csv format
  #
  # Args:
  #   data:           The data.table
  #   file.name:      The csv name
  #   analysis.name:  The name of the experiment e.g. "SILAC HUMAN"
  #   output.folder:  Default is "intermediate-data". The folder to store the csvs
  #
  # Returns:
  #   A csv with the data of the appropriate data.table
  #
  
  # Construct the output folder path
  output.folder.path <- paste("data-output", analysis.name, output.folder, sep = "/")
  
  # Change to the appropriate directory
  setwd(here(output.folder.path))
  
  # Substitute the dots to dashes
  dashed.file.name <- gsub("\\.", "-", file.name)
  
  # Add the csv extension to the name
  dashed.file.name <- paste0(dashed.file.name, ".csv")
  
  # Write the data.table to csv
  fwrite(data, file = dashed.file.name)
  
  # Inform the user
  cat(dashed.file.name, "saved under", here(output.folder.path), "\n")
  
  # Change back to the src directory
  setwd(here("src"))
}

filter.out.reverse.and.contaminants <- function(analysis.data) {
  #
  # Removes Contamintants, Reverse flagged peptides and only identified by site from the evidence data
  # 
  # Args:
  #   analysis.data: The analysis.data 
  #
  # Returns:
  #   The cleaned analysis.data data.table
  #
  
  # Copy the analysis data
  data <- copy(analysis.data)
  
  # Clean them from the Protein IDs starting with CON__ regarding the contaminants
  data.no.contaminants <- subset( data,
                                  grepl("^CON__",
                                        protein.ids,
                                        perl = TRUE) == FALSE)
  
  # Clean them from the Protein IDs starting with REV__ regarding the Reverse sequences
  data.no.contaminants.no.reverse <- subset(data.no.contaminants,
                                            grepl("^REV__",
                                                  protein.ids,
                                                  perl = TRUE) == FALSE)
  
  return (data.no.contaminants.no.reverse)
}

use.peptides.median <- function (intensities,
                                 do.norm = TRUE) {
  #
  # Get the median intensity of each peptide
  #
  # Args:
  #   intensities:        The list of the peptide's intensities between the technical replicates and fractions
  #   do.norm:            Default is TRUE. Should TODO
  #
  # Returns:
  #   The median intensity of each peptide 
  #
  
  # Unlist the intensities as it is a list of lists
  intensities <- unlist(intensities)
  
  # Remove the NAs
  intensities <- intensities[!is.na(intensities)]
  
  # If the number oif intensities is less than the threshold we assume that the peptide was not detected
  # Else  re return the median of the peptides intensity
  if (length(intensities) == 0) {
    return (list(NA))
  } else {
      return (list(median(intensities)))
  }
}

do.vsn.normalization <- function(filtered.data, conditions.to.compare,
                                 minimum.peptide.detections = 2, do.norm = TRUE) {
  #
  # Does variance stabilizing normalization (VSN) on the peptides intensities
  # 
  # Args:
  #   filtered.data:          The filtered data data.table
  #   conditions.to.compare:  The conditions to compare
  #   minimum.peptide.detections:     Default is 2. A strictness parameter regarding how many times is a peptide measured
  #                           If it is measures less than N times, the peptide is considered as detected but not measures
  #                           hence its value is set to NA
  #   do.norm:                Default is TRUE. Do the normalization or not? Only to use it for the plots before and 
  #                           after normalization
  #
  # Returns:
  #   The normalized analysis data.table
  #
  
  # Initialize the variable to return 
  vsn.normalized.data <- NULL
  
  # Make a copy of the analysis data
  data <- copy(filtered.data)
  
  # Keep only the information regarding the biological replicate for each peptide
  data[, description:= gsub("F.*$","", description, perl = TRUE)]
  
  # Rename the "description" column to "Biological Replicate"
  setnames(data, "description", "biological.replicate")
  
  # Reorder the data.table
  setkey(data, biological.replicate, protein.ids, unique.sequence.id)
  
  # Iterate over conditions
  for (condition in conditions.to.compare)  {
    
    # Replace the NaN values with NA
    set(data, which(is.nan(data[, eval(condition)])), condition, NA)
    
    # Transfrorm the data from long format to wide where we have the peptides on the y axis and the biological
    # replicates on the x axis
    wide.data <- dcast(data,
                       protein.ids + unique.sequence.id ~ biological.replicate,
                       value.var = condition,
                       fun.aggregate = use.peptides.median,
                       do.norm = do.norm)
    
    # Get only the part of the data.table with the intensities of the peptides
    vsn.matrix <- as.matrix(wide.data[, .SD, .SDcols = -c(1, 2)])
    
    # Do the normalization and suppress the warnings regardind the NA removal the the vsn package throws by default
    if (do.norm == TRUE) {
      suppressWarnings(vsn.matrix.normalized <- justvsn(vsn.matrix,verbose = FALSE))
    } else {
      vsn.matrix.normalized <- vsn.matrix
    }
    
    
    # Replace the old values of the matrix with the new
    for (index in 1:ncol(vsn.matrix)) {
      wide.data[, index + 2] <- vsn.matrix.normalized[, index]  
    }
    
    # Get the names of the columns
    no.condition.column.names <- colnames(wide.data)
    
    # Add the condition prefix
    condition.prefix.column.names <- paste(condition, no.condition.column.names[-c(1, 2)])
    
    # And reset the column names
    colnames(wide.data) <- c(no.condition.column.names[c(1,2)], condition.prefix.column.names)
    
    # If the return table is empty just copy the first data.table, otherwise add the second condition column
    if (is.null(vsn.normalized.data) == TRUE) {
      vsn.normalized.data <- copy(wide.data)
    } else {
      vsn.normalized.data <- merge(vsn.normalized.data,
                                   wide.data,
                                   by=c("protein.ids", "unique.sequence.id"))
    }
  }
  
  return (vsn.normalized.data)
}

do.peptide.intensities.plots <- function(not.normalized.data, vsn.normalized.data, analysis.name,
                                         plots.format = 5) {
  #
  # Makes 4 boxplot, before and after the normalization, with and without titles
  #
  # Args:
  #   intensities.before.normalization: The data.table of the non-normalized peptide intensities
  #   vsn.normalized.data:  The data.table of the normalized peptide intensities
  #   analysis.name:                    The name of the experiment e.g. "SILAC HUMAN"
  #   plots.format:           Default is 5 (jpeg format). A numeric value indicating the format of the plots. 
  #                           The numbers correspont to:  1     2     3     4      5      6     7     8     9
  #                                                     "eps" "ps"  "tex" "pdf" "jpeg" "tiff" "png" "bmp" "svg"
  #
  # Returns:
  #   The boxplots
  #
  
  # Copy the data.tables
  raw.data <- copy(not.normalized.data)
  normalized.data <- copy(vsn.normalized.data)
  
  # Remove the 'Protein IDs' columns
  raw.data[, protein.ids:= NULL]
  normalized.data[, protein.ids:= NULL]
  
  # Melt the data.tables
  raw.data.melted <- melt(raw.data,
                          id.vars = 1,
                          measure.vars = c(2:ncol(raw.data)) ,
                          variable.factor = FALSE)
  normalized.data.melted <- melt( normalized.data,
                                  id.vars = 1,
                                  measure.vars = c(2:ncol(normalized.data)) ,
                                  variable.factor = FALSE)
  
  # Substitute the space with newlines so the x-axis labels do not overlap
  raw.data.melted[, variable := gsub(" ","\n", variable)]
  normalized.data.melted[, variable := gsub(" ","\n", variable)]
  
  
  # Unlist the intensities column and for the non-normalized data log2 tranform them
  raw.data.melted[, value := log2(unlist(value))]
  normalized.data.melted[, value := unlist(value)]
  
  # For the raw data we may have Inf values for the 0 so in that case, we set the infs to NA
  infs <- which(is.finite(raw.data.melted$value) == FALSE)
  raw.data.melted$value[infs] <- NA 
  
  # Find the NAs
  raw.data.nas <- which(is.na(raw.data.melted$value) == TRUE)
  normalized.data.nas <- which(is.na(normalized.data.melted$value) == TRUE)
  
  # Remove the NAs
  raw.data.melted <- raw.data.melted[-c(raw.data.nas),]
  normalized.data.melted <- normalized.data.melted[-c(normalized.data.nas),]
  
  # Colorblind-friendly palette with grey:
  colorblind.palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # Add a column with the condition of each replicate
  raw.data.melted[, condition:= gsub("\n.*","", variable)]
  normalized.data.melted[, condition:= gsub("\n.*","", variable)]
  
  # Make the plots with no title
  intensities.before.normalization.plot.no.title <-  ggplot(data  = raw.data.melted,
                                                            aes(x = variable,
                                                                y = value)) +
    theme(legend.position = "none",
          plot.title      = element_text( hjust   = 0.5,
                                          family  = "Helvetica",
                                          size    = 12,
                                          face    = "bold"),
          axis.title.x    = element_text( family  = "Helvetica",
                                          size    = 8),
          axis.title.y    = element_text( family  = "Helvetica",
                                          size    = 8))+
    geom_boxplot(aes(fill  = as.factor(condition))) + 
    labs(x = "Replicates",
         y = "Peptide Intensities") +
    scale_fill_manual(values = colorblind.palette)
  
  intensities.after.normalization.plot.no.title <-  ggplot( data  = normalized.data.melted,
                                                            aes(x = variable,
                                                                y = value)) +
    theme(legend.position = "none",
          plot.title      = element_text( hjust   = 0.5,
                                          family  = "Helvetica",
                                          size    = 12,
                                          face    = "bold"),
          axis.title.x    = element_text( family  = "Helvetica",
                                          size    = 8),
          axis.title.y    = element_text( family  = "Helvetica",
                                          size    = 8))+
    geom_boxplot(aes(fill  = as.factor(condition))) + 
    labs(x = "Replicates",
         y = "Peptide Intensities") +
    scale_fill_manual(values = colorblind.palette)
  
  # Make the plots with titles
  intensities.before.normalization.plot.with.title <-  ggplot(data  = raw.data.melted,
                                                              aes(x = variable,
                                                                  y = value)) +
    theme(legend.position = "none",
          plot.title      = element_text( hjust   = 0.5,
                                          family  = "Helvetica",
                                          size    = 12,
                                          face    = "bold"),
          axis.title.x    = element_text( family  = "Helvetica",
                                          size    = 8),
          axis.title.y    = element_text( family  = "Helvetica",
                                          size    = 8))+
    geom_boxplot(aes(fill  = as.factor(condition))) + 
    labs(x = "Replicates",
         y = "Peptide Intensities",
         title = "Peptide Intensities Before Variance Stabilizing Normalization (VSN)") +
    scale_fill_manual(values = colorblind.palette)
  
  intensities.after.normalization.plot.with.title <-  ggplot( data  = normalized.data.melted,
                                                              aes(x = variable,
                                                                  y = value)) +
    theme(legend.position = "none",
          plot.title      = element_text( hjust   = 0.5,
                                          family  = "Helvetica",
                                          size    = 12,
                                          face    = "bold"),
          axis.title.x    = element_text( family  = "Helvetica",
                                          size    = 8),
          axis.title.y    = element_text( family  = "Helvetica",
                                          size    = 8))+
    geom_boxplot(aes(fill  = as.factor(condition))) + 
    labs(x = "Replicates",
         y = "Peptides' Intensities",
         title = "Peptide Intensities After Variance Stabilizing Normalization (VSN)") +
    scale_fill_manual(values = colorblind.palette)
  
  # Create the plots path for each analysis
  plots.path <- paste("data-output", analysis.name, "plots", sep = "/")
  
  # Change to the limma directory
  setwd(here(plots.path))
  
  # Initialize the plot format
  plot.format <- ""
  
  # In user has specified the he preffed different image format, change to the appropriate format
  switch(plots.format,
         plot.format <- "eps",
         plot.format <- "ps",
         plot.format <- "tex",
         plot.format <- "pdf",
         plot.format <- "jpeg",
         plot.format <- "tiff",
         plot.format <- "png",
         plot.format <- "bmp",
         plot.format <- "svg")
  
  # Save the boxplot before normalization with no description in high resolution 1920x1080 pixels, in the appropriate format
  ggsave(filename = paste0("plain-before-normalization-plot",".",plot.format),
         plot = intensities.before.normalization.plot.no.title,
         device = plot.format,
         width = 6.4,
         height = 3.6,
         units = "in",
         dpi = "print",
         limitsize = FALSE)
  
  # Save the boxplot after normalization with no description in high resolution 1920x1080 pixels, in the appropriate format
  ggsave(filename = paste0("plain-after-normalization-plot",".",plot.format),
         plot = intensities.after.normalization.plot.no.title,
         device = plot.format,
         width = 6.4,
         height = 3.6,
         units = "in",
         dpi = "print",
         limitsize = FALSE)
  
  # Save the boxplot before normalization with  description in high resolution 1920x1080 pixels, in the appropriate format
  ggsave(filename = paste0("descriptive-before-normalization-plot",".",plot.format),
         plot = intensities.before.normalization.plot.with.title,
         device = plot.format,
         width = 6.4,
         height = 3.6,
         units = "in",
         dpi = "print",
         limitsize = FALSE)
  
  # Save the boxplot after normalization with description in high resolution 1920x1080 pixels, in the appropriate format
  ggsave(filename = paste0("descriptive-after-normalization-plot",".",plot.format),
         plot = intensities.after.normalization.plot.with.title,
         device = plot.format,
         width = 6.4,
         height = 3.6,
         units = "in",
         dpi = "print",
         limitsize = FALSE)
  
  
  # Change back to the src directory
  setwd(here("src"))
}

do.knn.imputation <- function(vsn.normalized.data,
                              k.neighbours = 10, max.na.percentance.per.row = 50.0) {
  #
  # Does left-censored missing data imputation using quantile regression imputation of left-censored data method
  #
  # Args:
  #   vsn.normalized.data: The normalized data.table
  #
  # Returns:
  #   The imputated data.table
  #
  
  # Copy the initial normalized data
  imputed.data <- copy(vsn.normalized.data)
  
  # Get the column number 
  samples.column.number <- ncol(imputed.data)
  
  # Now for each row compute the NA percentance
  imputed.data[, na.percentance := ((samples.column.number - 2 - Reduce("+", lapply(.SD, is.na))) * 100 )/samples.column.number,
               .SDcols = c(3:samples.column.number)]
  
  # Then find the rows with correct observations
  rows.with.good.observations <- which(imputed.data$na.percentance > max.na.percentance.per.row)
  
  # Keep only them
  imputed.data <- imputed.data[rows.with.good.observations, ]
  
  # Remove the useless na.percentance column
  imputed.data <- imputed.data[, na.percentance := NULL]
  
  # Make a matrix with only the intensities
  data.to.impute <- as.matrix(imputed.data[, -c(1,2)])
  
  # Do the kNN imputation
  knn.imputation.results <- impute.knn(data.to.impute,
                                       k = k.neighbours, rowmax = max.na.percentance.per.row)
  
  # Get the new matrix
  imputed.matrix <- knn.imputation.results$data
  
  # And update the initial data.table
  for (index in 1:ncol(imputed.matrix)) {
    imputed.data[, index + 2] <- imputed.matrix[, index]  
  }
  
  return (imputed.data)  
}

do.peptides.aggregation <- function(imputed.datam
                                    minimum.peptide.detections = 2) {
  #
  # Aggregates the peptides' intensities in order to get the protein abundance in each condition/biological replicate
  #
  # Args:
  #   imputed.data; The data.table with the imputed data
  #
  # Returns:
  #   The data.table with the protein abundances for each protein in each condition/biological replicate
  #
  
  # Make a copy of the imputed data
  aggregated.data <- copy(imputed.data)
  
  # Count how many peptides where detected for the protein
  aggregated.data <- aggregated.data[, peptides.detected := .N, by = protein.ids]
  
  # Find which rows fulfill the minimum peptides rule
  proteins.the.threshold <- which(aggregated.data$peptides.detected > minimum.peptide.detections)
  
  # And keep only the proteins for which the minimum peptides threshold is fulfilled
  aggregated.data <- aggregated.data[proteins.the.threshold,]
  
  # Delete the column
  aggregated.data[, peptides.detected := NULL]
  
  # Do the aggregation of the peptides' intensity that belong to the same protein
  aggregated.data <- aggregated.data[, lapply(.SD,
                                              sum),
                                     .SDcols=!"unique.sequence.id",
                                     by="protein.ids"]
  
  # And rename the columns
  old.column.names <- colnames(aggregated.data)[-c(1)]
  
  # In order to  know that
  new.column.names <- paste("protein.abundance", old.column.names)
  
  # The columns now hold the protein abundances
  colnames(aggregated.data)[-c(1)] <- new.column.names
  
  return (aggregated.data)
}

do.QQ.plots <- function(aggregated.data, conditions.to.compare, analysis.name, experimental.metadata,
                        plots.format = 5) {
  #
  # Produces 2 QQ plots, one with tile and one with no title
  #
  # Args:
  #   aggregated.data:        The aggregated.data data.table
  #   conditions.to.compare:  The condition for the descriptive plot
  #   analysis.name:          The name of the experiment e.g. "SILAC HUMAN"
  #   experimental.metadata:  The metadata for each condition, number of biological/technical replicates
  #   plots.format:           Default is 5 (jpeg format). A numeric value indicating the format of the plots. 
  #                           The numbers correspont to:  1     2     3     4      5      6     7     8     9
  #                                                     "eps" "ps"  "tex" "pdf" "jpeg" "tiff" "png" "bmp" "svg"
  # Returns:
  #   The saved plots inside the data-out/plots folder
  #
  
  # Copy the data.table
  data <- copy(aggregated.data)
  
  # Remove the protein ids column
  data[, "protein.ids":= NULL]
  
  # Get the conditions' names
  condition.A <- conditions.to.compare[[1]]
  condition.B <- conditions.to.compare[[2]]
  
  # In case of an Isotopic Experiment we have only one element in the list
  if (length(names(experimental.metadata)) == 1) {
    condition.A <- "Labeled Experiment"
    condition.B <- "Labeled Experiment"
  }
  
  # Get the number of biological replicates for each condition
  condition.A.number.of.biological <- experimental.metadata[[condition.A]]$number.of.biological.replicates
  
  # Get the number of technical replicates for each condition
  condition.A.number.of.technical <- experimental.metadata[[condition.A]]$number.of.technical.replicates
  
  # Colorblind-friendly palette with grey:
  colorblind.palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # Get the data for condition A
  condition.A.data <- as.double(as.matrix(data[, 1:(condition.A.number.of.biological * condition.A.number.of.technical)]))
  
  # Tag them
  condition.A.data <- cbind(condition.A.data, 1)
  
  # Get the data for condition B
  condition.B.data <- as.double(as.matrix(data[, ((condition.A.number.of.biological * condition.A.number.of.technical) + 1): ncol(data)]))
  
  # Tag them
  condition.B.data <- cbind(condition.B.data, 2)
  
  # Wrap them up in a data.frame
  qq.plot.data <-  as.data.frame(rbind(condition.A.data, condition.B.data))  
  
  # Fix the column names
  colnames(qq.plot.data) <- c("abundancies", "condition")
  
  # Make a QQ plot tha does not contain any title or description
  qq.plot.no.title <- ggplot(data       = qq.plot.data,
                             aes(sample  = abundancies,
                                 group   = condition,
                                 color   = as.factor(condition))) +
    stat_qq(alpha  = 0.4,
            size   = 1.75) +
    stat_qq_line(alpha  = 0.4,
                 size   = 1.75) +
    theme(legend.position = "bottom",
          axis.title.x = element_text( family = "Helvetica",
                                       size   = 8),
          axis.title.y = element_text( family = "Helvetica",
                                       size   = 8)) +
    labs(x = "Theoretical Quantiles",
         y = "Sample Quantiles") +
    scale_color_manual(values = colorblind.palette,
                       name = "conditions",
                       labels = conditions.to.compare)
  
  # Now make a title for the descriptive QQ plot
  qq.plot.title <- paste0( "Q-Q plot of conditions ",
                           conditions.to.compare[[1]],
                           " and ",
                           conditions.to.compare[[2]],
                           " protein abundancies")
  
  # In case that the description is longer than 70 characters, break the title into multiple lines
  qq.plot.title <- paste0(strwrap(qq.plot.title, 70), collapse = "\n")      
  
  # Make a QQ plot that contains title or description
  qq.plot.with.title <- ggplot(data       = qq.plot.data,
                               aes(sample  = abundancies,
                                   group   = condition,
                                   color   = as.factor(condition))) +
    stat_qq(alpha  = 0.4,
            size   = 1.75) +
    stat_qq_line(alpha  = 0.4,
                 size   = 1.75) +
    theme(legend.position = "bottom",
          plot.title      = element_text( hjust   = 0.5,
                                          family  = "Helvetica",
                                          size    = 12,
                                          face    = "bold"),
          axis.title.x = element_text( family = "Helvetica",
                                       size   = 8),
          axis.title.y = element_text( family = "Helvetica",
                                       size   = 8)) +
    labs(x = "Theoretical Quantiles",
         y = "Sample Quantiles",
         title = qq.plot.title) +
    scale_color_manual(values = colorblind.palette,
                       name = "Conditions",
                       labels = conditions.to.compare)
  
  # Create the plots path for each analysis
  plots.path <- paste("data-output", analysis.name, "plots", sep = "/")
  
  # Change to the plots directory
  setwd(here(plots.path))
  
  # Initialize the plot format
  plot.format <- ""
  
  # In user has specified the he preffed different image format, change to the appropriate format
  switch(plots.format,
         plot.format <- "eps",
         plot.format <- "ps",
         plot.format <- "tex",
         plot.format <- "pdf",
         plot.format <- "jpeg",
         plot.format <- "tiff",
         plot.format <- "png",
         plot.format <- "bmp",
         plot.format <- "svg")
  
  # Save the volcano plot with no description in high resolution 1920x1080 pixels, in the appropriate format
  ggsave(filename = paste0("plain-qq-plot",".",plot.format),
         plot = qq.plot.no.title,
         device = plot.format,
         width = 6.4,
         height = 3.6,
         units = "in",
         dpi = "print",
         limitsize = FALSE)
  
  
  # Save the descriptive volcano plot in high resolution 1920x1080 pixels, in the appropriate format
  ggsave(filename = paste0("descriptive-qq-plot",".",plot.format),
         plot = qq.plot.with.title,
         device = plot.format,
         width = 6.4,
         height = 3.6,
         units = "in",
         dpi = "print",
         limitsize = FALSE)
  
  # Change back to the src directory
  setwd(here("src"))
}

do.fold.change.histogram <- function(limma.results, conditions.to.compare, analysis.name,
                                     plots.format = 5) {
  #
  # Produces 2 log fold change histograms, one with tile and one with no title
  #
  # Args:
  #   limma.results:          The data.frame with the limma results
  #   conditions.to.compare:  The condition for the descriptive plot
  #   analysis.name:          The name of the experiment e.g. "SILAC HUMAN"
  #   plots.format:           Default is 5 (jpeg format). A numeric value indicating the format of the plots. 
  #                           The numbers correspont to:  1     2     3     4      5      6     7     8     9
  #                                                     "eps" "ps"  "tex" "pdf" "jpeg" "tiff" "png" "bmp" "svg"
  # Returns:
  #   The saved plots inside the data-out/plots folder
  #
  
  # The histogram without title
  fold.change.histogram.no.title <-  ggplot(data = limma.results,
                                            aes(x = logFC)) +
                                     geom_histogram(color = "#999999",
                                                    fill  = "#E69F00",
                                                    bins  = 40) +
                                     labs(x = "Log Fold-Change",
                                          y = "Frequency")
  
  # The histogram with title
  fold.change.histogram.with.title <-   ggplot( data  = limma.results,
                                                aes(x = logFC)) +
                                        geom_histogram(color = "#999999",
                                                       fill  = "#E69F00",
                                                       bins  = 40) +
                                        theme(plot.title = element_text(hjust   = 0.5,
                                                                        family  = "Helvetica",
                                                                        size    = 12,
                                                                        face    = "bold"),
                                            axis.title.x = element_text( family = "Helvetica",
                                                                         size   = 8),
                                            axis.title.y = element_text( family = "Helvetica",
                                                                         size   = 8)) +
                                         labs(x = "Log2 Fold-Change",
                                              y = "Frequency",
                                              title = paste0("Log2 Fold-Change ",
                                                             conditions.to.compare[[1]],
                                                             "/",
                                                             conditions.to.compare[[2]],
                                                             " histogram"))
  
  # Create the plots path for each analysis
  plots.path <- paste("data-output", analysis.name, "plots", sep = "/")
  
  # Change to the plots directory
  setwd(here(plots.path))
  
  # Initialize the plot format
  plot.format <- ""
  
  # In user has specified the he preffed different image format, change to the appropriate format
  switch(plots.format,
         plot.format <- "eps",
         plot.format <- "ps",
         plot.format <- "tex",
         plot.format <- "pdf",
         plot.format <- "jpeg",
         plot.format <- "tiff",
         plot.format <- "png",
         plot.format <- "bmp",
         plot.format <- "svg")
  
  # Save the volcano plot with no description in high resolution 1920x1080 pixels, in the appropriate format
  ggsave(filename = paste0("plain-fold-change-histogram",".",plot.format),
         plot = fold.change.histogram.no.title,
         device = plot.format,
         width = 6.4,
         height = 3.6,
         units = "in",
         dpi = "print",
         limitsize = FALSE)
  
  
  # Save the descriptive volcano plot in high resolution 1920x1080 pixels, in the appropriate format
  ggsave(filename = paste0("descriptive-fold-change-histogram",".",plot.format),
         plot = fold.change.histogram.with.title,
         device = plot.format,
         width = 6.4,
         height = 3.6,
         units = "in",
         dpi = "print",
         limitsize = FALSE)
  
  # Change back to the src directory
  setwd(here("src"))
  
}

do.value.ordered.ratio.plot <- function(limma.results, conditions.to.compare, analysis.name,
                                        plots.format = 5) {
  #
  # Does 2 value ordered ratio plots, one with title and one with not.
  #
  # Args:
  #   limma.results:          The data.frame with the limma results
  #   conditions.to.compare:  The condition for the descriptive plot
  #   analysis.name:          The name of the experiment e.g. "SILAC HUMAN"
  #   plots.format:           Default is 5 (jpeg format). A numeric value indicating the format of the plots. 
  #                           The numbers correspont to:  1     2     3     4      5      6     7     8     9
  #                                                     "eps" "ps"  "tex" "pdf" "jpeg" "tiff" "png" "bmp" "svg"
  # Returns:
  #   The saved plots inside the data-out/plots folder
  
  
  # Get the limma results
  ratio.plot.data <- limma.results 
  
  # Order them by logFC
  ratio.plot.data <- ratio.plot.data[order(ratio.plot.data$logFC),]
  
  # Calculate the upper and the lower value
  ratio.plot.data$lower <- ratio.plot.data$logFC - log2(ratio.plot.data$AveExpr)
  ratio.plot.data$upper <- ratio.plot.data$logFC + log2(ratio.plot.data$AveExpr)
  
  # Colorblind-friendly palette with grey:
  colorblind.palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # Make the plot with no title
  ratio.plot.no.title <- ggplot( data      = ratio.plot.data,
                                 aes(x     = 1:length(Protein),
                                     y     = logFC, 
                                     color = is.diffenetially.expressed))+
    geom_point(alpha = 0.7,
               size  = 1.5) +
    geom_errorbar(aes(ymin=lower, ymax=upper),
                  width = 2.5,
                  size  = 1) +
    theme(legend.position="bottom",
          axis.title.x = element_text( family = "Helvetica",
                                       size   = 8),
          axis.title.y = element_text( family = "Helvetica",
                                       size   = 8)) +
    labs(x ="Protein IDs",
         y = paste0("Average logFC(",
                    conditions.to.compare[1],
                    "/",
                    conditions.to.compare[2],
                    ")")) +
    scale_color_manual(values = colorblind.palette,
                       name   = "Protein is Differentially Expressed")
  
  # Make the plot with the title
  ratio.plot.with.title <- ggplot( data      = ratio.plot.data,
                                   aes(x     = 1:length(Protein),
                                       y     = logFC, 
                                       color = is.diffenetially.expressed))+
    geom_point(alpha = 0.7,
               size  = 1.5) +
    geom_errorbar(aes(ymin = lower,
                      ymax = upper),
                  width = 2.5,
                  size  = 1) +
    theme(legend.position="bottom",
          plot.title      = element_text( hjust   = 0.5,
                                          family  = "Helvetica",
                                          size    = 12,
                                          face    = "bold"),
          axis.title.x = element_text( family = "Helvetica",
                                       size   = 8),
          axis.title.y = element_text( family = "Helvetica",
                                       size   = 8)) +
    labs(x ="Protein IDs",
         y = paste0("Average logFC(",
                    conditions.to.compare[1],
                    "/",
                    conditions.to.compare[2],
                    ")"),
         title = "Value Ordered fold-change plot") +
    scale_color_manual(values = colorblind.palette,
                       name   = "Protein is Differentially Expressed")
  
  # Create the plots path for each analysis
  plots.path <- paste("data-output", analysis.name, "plots", sep = "/")
  
  # Change to the plots directory
  setwd(here(plots.path))
  
  # Initialize the plot format
  plot.format <- ""
  
  # In user has specified the he preffed different image format, change to the appropriate format
  switch(plots.format,
         plot.format <- "eps",
         plot.format <- "ps",
         plot.format <- "tex",
         plot.format <- "pdf",
         plot.format <- "jpeg",
         plot.format <- "tiff",
         plot.format <- "png",
         plot.format <- "bmp",
         plot.format <- "svg")
  
  # Save the volcano plot with no description in high resolution 1920x1080 pixels, in the appropriate format
  ggsave(filename = paste0("plain-ratio-plot",".",plot.format),
         plot = ratio.plot.no.title,
         device = plot.format,
         width = 6.4,
         height = 3.6,
         units = "in",
         dpi = "print",
         limitsize = FALSE)
  
  
  # Save the descriptive volcano plot in high resolution 1920x1080 pixels, in the appropriate format
  ggsave(filename = paste0("descriptive-ratio-plot",".",plot.format),
         plot = ratio.plot.with.title,
         device = plot.format,
         width = 6.4,
         height = 3.6,
         units = "in",
         dpi = "print",
         limitsize = FALSE)
  
  # Change back to the src directory
  setwd(here("src"))
  
}

do.MA.plots <- function(limma.results, conditions.to.compare, analysis.name,
                        plots.format = 5) {
  #
  # Produces 2 MA plots, one with tile and one with no title
  #
  # Args:
  #   limma.results:          The data.frame with the limma results
  #   conditions.to.compare:  The condition for the descriptive plot
  #   analysis.name:          The name of the experiment e.g. "SILAC HUMAN"
  #   plots.format:           Default is 5 (jpeg format). A numeric value indicating the format of the plots. 
  #                           The numbers correspont to:  1     2     3     4      5      6     7     8     9
  #                                                     "eps" "ps"  "tex" "pdf" "jpeg" "tiff" "png" "bmp" "svg"
  # Returns:
  #   The saved plots inside the data-out/plots folder
  #
  
  # Get the conditions' names
  condition.A <- conditions.to.compare[[1]]
  condition.B <- conditions.to.compare[[2]]
  
  # Colorblind-friendly palette with grey:
  colorblind.palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # Make a MA plot tha does not contain any title or description
  ma.plot.no.title <- ggplot(data       = limma.results,
                             aes(x      = AveExpr,
                                 y      = logFC,
                                 color  = is.diffenetially.expressed)) +
    geom_point(alpha  = 0.4,
               size   = 1.75) +
    theme(legend.position="none",
          axis.title.x = element_text( family = "Helvetica",
                                       size   = 8),
          axis.title.y = element_text( family = "Helvetica",
                                       size   = 8)) +
    labs(x ="A (Average Protein Abundance)",
         y = paste0("M (log2 ", condition.A, "/", condition.B,")")) +
    scale_color_manual(values = colorblind.palette)
  
  # Make a MA plot that with a description
  ma.plot.with.title <- ggplot(data       = limma.results,
                               aes(x      = AveExpr,
                                   y      = logFC,
                                   color  = is.diffenetially.expressed)) +
    geom_point(alpha  = 0.4,
               size   = 1.75) +
    theme(legend.position="none",
          plot.title      = element_text( hjust   = 0.5,
                                          family  = "Helvetica",
                                          size    = 12,
                                          face    = "bold"),
          axis.title.x    = element_text( family  = "Helvetica",
                                          size    = 8),
          axis.title.y    = element_text( family  = "Helvetica",
                                          size    = 8))+
    labs(x     = "A (Average Protein Abundance)",
         y     = paste0("M (log2 ", condition.A, "/", condition.B,")"),
         title = "MA-plot") +
    scale_color_manual(values = colorblind.palette)
  
  # Create the plots path for each analysis
  plots.path <- paste("data-output", analysis.name, "plots", sep = "/")
  
  # Change to the plots directory
  setwd(here(plots.path))
  
  # Initialize the plot format
  plot.format <- ""
  
  # In user has specified the he preffed different image format, change to the appropriate format
  switch(plots.format,
         plot.format <- "eps",
         plot.format <- "ps",
         plot.format <- "tex",
         plot.format <- "pdf",
         plot.format <- "jpeg",
         plot.format <- "tiff",
         plot.format <- "png",
         plot.format <- "bmp",
         plot.format <- "svg")
  
  # Save the volcano plot with no description in high resolution 1920x1080 pixels, in the appropriate format
  ggsave(filename = paste0("plain-ma-plot",".",plot.format),
         plot = ma.plot.no.title,
         device = plot.format,
         width = 6.4,
         height = 3.6,
         units = "in",
         dpi = "print",
         limitsize = FALSE)
  
  
  # Save the descriptive volcano plot in high resolution 1920x1080 pixels, in the appropriate format
  ggsave(filename = paste0("descriptive-ma-plot",".",plot.format),
         plot = ma.plot.with.title,
         device = plot.format,
         width = 6.4,
         height = 3.6,
         units = "in",
         dpi = "print",
         limitsize = FALSE)
  
  # Change back to the src directory
  setwd(here("src"))
}

do.volcano.plots <- function(limma.results, conditions.to.compare, analysis.name,
                             plots.format = 5, error.correction.method = "B", fold.change.cut.off = 1.5, FDR = 0.05) {
  #
  # Generates and saves the 2 volcano plots in the appropriate format, one plain and one with descriptive info
  #
  # Args:
  #   limma.results:          The data.frame with the limma results
  #   conditions.to.compare:  The condition for the descriptive plot
  #   analysis.name:          The name of the experiment e.g. "SILAC HUMAN"
  #   plots.format:           Default is 5 (jpeg format). A numeric value indicating the format of the plots. 
  #                           The numbers correspont to:  1     2     3     4      5      6     7     8     9
  #                                                     "eps" "ps"  "tex" "pdf" "jpeg" "tiff" "png" "bmp" "svg"
  #   error.correction.method:  Default is "B". Corrects the Type 1 errors using the Bonferroni correction method or
  #                             the Benjamini-Hochberg method. Values can be "B" or "BH'
  #   fold.change.cut.off:      Default is 1.5. The least fold change  in order to consider a protein differentially expressed
  #   FDR:                      Default is 0.05. The least fold change  in order to consider a protein differentially expressed
  #
  # Returns:
  #   The saved plots inside the data-out/plots folder
  #
  
  # Store the limma results
  plot.data <- limma.results
  
  # Get the total number of proteins in the experiment
  number.of.total.proteins <- dim(limma.results)[1]
  
  # Count the differentially expressed proteins
  number.of.differentially.expressed.proteins <- length(which(plot.data$is.diffenetially.expressed == TRUE))
  
  # Colorblind-friendly palette with grey:
  colorblind.palette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  
  # Make a volcano plot tha does not contain any title or description
  volcano.plot.no.title <- ggplot(data       = plot.data,
                                  aes(x      = logFC,
                                      y      = -log10(adj.P.Val),
                                      color  = is.diffenetially.expressed)) +
    geom_point(alpha  = 0.4,
               size   = 1.75) +
    theme(legend.position="none",
          axis.title.x = element_text( family = "Helvetica",
                                       size   = 8),
          axis.title.y = element_text( family = "Helvetica",
                                       size   = 8)) +
    labs(x ="log2 fold change",
         y = "-log10 p-value") +
    scale_color_manual(values = colorblind.palette)
  
  # Now make a title for the descriptive volcano plot
  volcano.title <- paste0( "Differentially Expressed Proteins between conditions ",
                           conditions.to.compare[[1]],
                           " and ",
                           conditions.to.compare[[2]])
  
  # In case that the description is longer than 70 characters, break the title into multiple lines
  volcano.title <- paste0(strwrap(volcano.title, 70), collapse = "\n")
  
  if (error.correction.method == "B") {
    p.value.correction.method <- "Bonferroni"
  } else {
    p.value.correction.method <- "Benjamini–Hochberg"
  }
  
  # Do the same for the subtitle
  volcano.subtitle <- paste0(number.of.differentially.expressed.proteins,
                             "(",
                             round(number.of.differentially.expressed.proteins/number.of.total.proteins, 2) * 100,
                             "%)  differentially expressed proteins in a total of ",
                             number.of.total.proteins, 
                             " detected by MS",
                             " (logFC.cutoff = ",
                             fold.change.cut.off,
                             " , FDR = ",
                             FDR,
                             ", p-value correction method = ",
                             p.value.correction.method,
                             ")")
  
  # Do the same for the subtitle
  volcano.subtitle <- paste0(strwrap(volcano.subtitle, 90), collapse = "\n")
  
  # # Make a volcano plot with a title and a description
  volcano.plot.with.title <- ggplot(data     = plot.data,
                                    aes(x      = logFC,
                                        y      = -log10(adj.P.Val),
                                        color  = is.diffenetially.expressed)) +
    geom_point(alpha  = 0.4,
               size   = 1.75) +
    theme( legend.position ="none", 
           plot.title      = element_text( hjust   = 0.5,
                                           family  = "Helvetica",
                                           size    = 12,
                                           face    = "bold"),
           plot.subtitle   = element_text( hjust   = 0.5,
                                           family  = "Helvetica",
                                           size    = 10),
           axis.title.x    = element_text( family  = "Helvetica",
                                           size    = 8),
           axis.title.y    = element_text( family  = "Helvetica",
                                           size    = 8))+
    labs(x          = "log2 fold change",
         y          = "-log10 p-value",
         title      = volcano.title,
         subtitle   = volcano.subtitle) +
    scale_color_manual(values = colorblind.palette)
  
  # Create the plots path for each analysis
  plots.path <- paste("data-output", analysis.name, "plots", sep = "/")
  
  # Change to the plots directory
  setwd(here(plots.path))
  
  # Initialize the plot format
  plot.format <- ""
  
  # In user has specified the he preffed different image format, change to the appropriate format
  switch(plots.format,
         plot.format <- "eps",
         plot.format <- "ps",
         plot.format <- "tex",
         plot.format <- "pdf",
         plot.format <- "jpeg",
         plot.format <- "tiff",
         plot.format <- "png",
         plot.format <- "bmp",
         plot.format <- "svg")
  
  # Save the volcano plot with no description in high resolution 1920x1080 pixels, in the appropriate format
  ggsave(filename = paste0("plain-volcano-plot",".",plot.format),
         plot = volcano.plot.no.title,
         device = plot.format,
         width = 6.4,
         height = 3.6,
         units = "in",
         dpi = "print",
         limitsize = FALSE)
  
  
  # Save the descriptive volcano plot in high resolution 1920x1080 pixels, in the appropriate format
  ggsave(filename = paste0("descriptive-volcano-plot",".",plot.format),
         plot = volcano.plot.with.title,
         device = plot.format,
         width = 6.4,
         height = 3.6,
         units = "in",
         dpi = "print",
         limitsize = FALSE)
  
  # Change back to the src directory
  setwd(here("src"))
  
}

do.limma.analysis <- function(aggregated.data, conditions.to.compare, experimental.metadata,
                              error.correction.method = "BH", fold.change.cut.off = 1.5, FDR = 0.05) {
  #
  # Does the limma analysis and returns a dataframe with the results of limme
  #
  # Args:
  #   aggregated.data:          The protein abundancis data.table
  #   conditions.to.compare:    The condition to compare
  #   experimental.metadata:    The metadata list of the experiment (condition >  number of biological replicates/
  #                                                                               number of technical replicates/
  #                                                                               number of fractions)
  #   error.correction.method:  Default is "BH". Corrects the Type 1 errors using the Bonferroni correction method or
  #                             the Benjamini-Hochberg method. Values can be "B" or "BH'
  #   fold.change.cut.off:      Default is 1.5. The least fold change  in order to consider a protein differentially expressed
  #   FDR:                      Default is 0.05. The least fold change  in order to consider a protein differentially expressed
  #   
  # Returns:
  #   A data.frame with the results of the limma analysis
  #
  
  # Copy the data
  limma.data <- copy(aggregated.data)
  
  # Get the protein names
  protein.ids <- limma.data[, protein.ids]
  
  # And remove the protein names from the data.table
  limma.data[, protein.ids:= NULL]
  
  # Now convert the data.table to matrix
  limma.data <- as.matrix(limma.data)
  
  # And add the protein names as rownames  
  rownames(limma.data) <- protein.ids
  
  # Then get the colnames
  limma.column.names <- colnames(limma.data)
  
  # Remove the "protein.abundance." prefic and the technical replicate prefix
  limma.column.names <- gsub("^(protein.abundance )|(T[[:digit:]*]$)", "",
                             limma.column.names,
                             perl = TRUE)
  
  # Reset the column names
  colnames(limma.data) <- limma.column.names
  
  #  Get the condition names
  condition.A <- conditions.to.compare[1]
  condition.B <- conditions.to.compare[2]
  
  # In case of an Isotopic Experiment we have only one element in the list
  if (length(names(experimental.metadata)) == 1) {
    condition.A <- "Labeled Experiment"
    condition.B <- "Labeled Experiment"
  }
  
  # Get the number of biological replicates for each condition
  condition.A.number.of.biological <- experimental.metadata[[condition.A]]$number.of.biological.replicates
  condition.B.number.of.biological <- experimental.metadata[[condition.B]]$number.of.biological.replicates
  
  # Get the number of technical replicates for each condition
  condition.A.number.of.technical <- experimental.metadata[[condition.A]]$number.of.technical.replicates
  condition.B.number.of.technical <- experimental.metadata[[condition.B]]$number.of.technical.replicates
  
  # Construct the design matrix
  design.condition.A <- c(rep.int(1, condition.A.number.of.biological * condition.A.number.of.technical),
                          rep.int(0, condition.B.number.of.biological * condition.B.number.of.technical))
  
  design.condition.B <- c(rep.int(0, condition.A.number.of.biological * condition.A.number.of.technical),
                          rep.int(1, condition.B.number.of.biological * condition.B.number.of.technical))
  
  #  Get the condition names
  condition.A <- conditions.to.compare[1]
  condition.B <- conditions.to.compare[2]
  
  design <- matrix(c(design.condition.A, design.condition.B), ncol = 2)
  
  # Add colnames to the design matrix
  colnames(design) <- c(condition.A, condition.B)
  
  # In we also have technical replicates, add the appropriate column to the design
  if (condition.A.number.of.technical != 1 &
      condition.B.number.of.technical != 1) {
    # design <- cbind(design, technical.replicate = c(rep.int(0:1, condition.A.number.of.biological),
    #                                               rep.int(0:1, condition.B.number.of.biological)))

    # Prepare the limma block matrix
    limma.block <- c(rep(1:condition.A.number.of.biological,
                         c( rep.int(condition.A.number.of.technical,
                                    condition.A.number.of.biological))),
                     rep(1:condition.B.number.of.biological,
                         c( rep.int(condition.B.number.of.technical,
                                    condition.B.number.of.biological))))
    
    # Calculate the correlation between the technical replicates
    duplicate.correlation <- duplicateCorrelation(limma.data,
                                                  design,
                                                  block = limma.block)
    
    # Get the correlation
    duplicate.correlation.consensus <- duplicate.correlation$consensus.correlation
    
    # Do the linear modelling taking into account the technical replicates and their correlation
    limma.fit <- lmFit(limma.data,
                       design,
                       block = limma.block,
                       correlation = duplicate.correlation.consensus)
  
  } else {
    
    # Do the linear modelling in the case we do not have technical replicates
    limma.fit <- lmFit(limma.data, design)
  }
  
  contrast.matrix <- makeContrasts(paste0(conditions.to.compare[1], "-", conditions.to.compare[2]),
                                   levels = c(conditions.to.compare[1], conditions.to.compare[2]))
  
  limma.fit.contrast <- contrasts.fit(limma.fit, contrast.matrix)
  
  # And calculates the statistics
  limma.fit.treat <- treat(limma.fit.contrast,
                             lfc = log2(fold.change.cut.off),
                             trend = TRUE,
                             robust = TRUE)

  # Make a data.frame with the statistics for each protein, sorted by p-value, either it is significant or not
  limma.results <- topTreat(limma.fit.treat,
                            adjust.method = error.correction.method,
                            sort.by = "P",
                            n = Inf)
  
  # And a column with the protein names
  limma.results <- cbind(Protein=rownames(limma.results), limma.results)
  
  # And restore the rownames to numbers
  rownames(limma.results) <- 1:nrow(limma.results)
  
  # Depending on the the results conservativeness, use the appropriate method
  limma.results$is.diffenetially.expressed =  (abs(limma.results$logFC) > fold.change.cut.off) &
                                              limma.results$adj.P.Val  < FDR

  return (limma.results)
}

# 
# tmp.evalutate.correctness.t.test <- function(aggregated.data, proteins=c("PA3479", "PA5346", "PA0176", "PA4625", "PA4624", "PA3724", "PA3385", "PA5272", "PA3217", "PA3544", "PA5060")) {
#   results <- lapply(proteins, function(proteinID) {
#     test.data <- aggregated.data[grep(proteinID, aggregated.data[[1]]), ]
#     group1 <- unlist(test.data[1,2:7])
#     group2 <- unlist(test.data[1,8:13])
#     return (t.test(group1, group2)$"p.value")
#   })
#   results <- unlist(results)
#   significant.proteins <- proteins[results < 0.05]
#   number.of.significant <- length(significant.proteins)
#   cat(toString(number.of.significant),"proteins where significant:", significant.proteins)
# }
# 
# tmp.evalutate.correctness.limma <- function(limma.results, proteins=c("PA3479", "PA5346", "PA0176", "PA4625", "PA4624", "PA3724", "PA3385", "PA5272", "PA3217", "PA3544", "PA5060")) {
#   results <- c()
#   for (protein in proteins) {
#     id <- grep(protein, limma.results$Protein)
#     results <- c(results, limma.results$is.diffenetially.expressed[id])
#   }
#   cat(length(proteins[results]), "proteins where significant:",proteins[results],"\n")
# }
# }
