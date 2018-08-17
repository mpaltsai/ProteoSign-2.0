make.folder <- function(folder.name) {
  
  data.output.path <- here("data-output")
  
  folder.output <- paste(data.output.path, folder.name, sep = "/")
  
  # If folder already exist, remove it
  if( file.exists(folder.output) == TRUE) {
    unlink(folder.output, recursive=T, force=T)
  }
  
  # And then create the new folder
  tryCatch(
    {
      dir.create(folder.output)
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

make.data.output.folders <- function() {
  #
  # Makes the data output folders limma-output, plots 
  #
  # Args:
  #
  # Returns:
  #   Prints a message if the folder was created successfully, otherwise terminates the analysis.
  #
  # Set limma output path
  data.output.path <- here("data-output")
  
  # Set limma output folder name
  folders.to.create <- c("limma-output", "plots", "intermediate-data")
  
  invisible(lapply(folders.to.create, make.folder))
  
}

make.Venn.diagram <- function(analysis.data, conditions.to.compare, analysis.title) {
  #
  # Makes a Venn diagram between 2 conditions
  #
  # Args:
  #   analysis.data:          The data.table with the evidence data
  #   conditions.to.compare:  A vector with the 2 conditions to compare
  #   analysis.title:         The analysis title provided by the user
  #
  # Returns:
  #   A venn diagram in png format and the data needed for its construction
  #
  
  # Subset the evidence data keeping only the "Protein IDs", "Condition","biological replicate",
  # "technical replicate", "fraction" columns
  venn.table <- evidence.data[,
                              .SD,
                              .SDcols =  c("Protein IDs",
                                           "Condition",
                                          "biological replicate",
                                          "technical replicate",
                                          "fraction")]
  
  # Order the subseted data.table
  setkey(venn.table,  "Protein IDs",
         "Condition",
         "biological replicate",
         "technical replicate",
         "fraction")
  
  # Set limma output path
  data.output.path <- here("data-output")
  
  # Set limma output folder name
  limma.folder <- "limma-output"
  
  # Prepare the limma output path
  limma.output <- paste(data.output.path, limma.folder, sep = "/")
  
  # Change directory to limma output
  setwd(limma.output)
  
  # Count the occurences of each protein in each condition
  venn.table <- venn.table[, .("Occurences"=.N), 
                           by=c( "Protein IDs",
                                  "Condition")]
  
  # Order the data.table
  setkey(venn.table, Condition)
  
  # Write the data in a TSV file
  write.table(venn.table,
              file = paste0(analysis.title, 
                          "-venn-data.txt"),
              sep="\t",
              row.names = FALSE,
              quote = FALSE)
  
  # Now store the 2 conditions
  conditions <- unique(venn.table[, Condition])
  
  # Store in the variable the condition A
  condition.1.table <- venn.table[conditions[1]]
  
  # Store in the variable the condition B
  condition.2.table <- venn.table[conditions[2]]
  
  # Take the Protein IDs and make the result a vector
  condition.1.proteins <- unlist(condition.1.table[, "Protein IDs"])
  condition.2.proteins <- unlist(condition.2.table[, "Protein IDs"])
  
  # Make the set
  sets <- list(condition.1.proteins, condition.2.proteins)
  
  # Give the names of the conditions to the list's set
  names(sets) <- conditions
  
  # By default the VennDiagram package makes a log file, that we want to supperss
  futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
  
  # Now draw the venn diagram
  venn.diagram <- venn.diagram(sets,
                               "test.venn.png",
                               imagetype = "png",
                               main = paste0("Venn diagram between the conditions ",
                                            conditions[1],
                                            " and ",
                                            conditions[2],
                                            "."),
                               category = conditions,
                               scaled = FALSE,
                               alpha = c(0.7, 0.7),
                               fill = c("blue", "red"),
                               cat.default.pos = "text",
                               cat.pos = c(1,1),
                               cex =  2,
                               cat.cex = 2,
                               cat.dist = c(0.05, 0.05))
  
  # And finally reset the working directory to the scripts folder
  setwd(here("src"))  
}

save.intermediate.data.tables <- function(data, file.name, output.folder = "intermediate-data") {
  #
  # Saves the intermediate data.table into csv format
  #
  # Args:
  #   data:           The data.table
  #   file.name:      The csv name
  #   output.folder:  Default is "intermediate-data". The folder to store the csvs
  #
  # Returns:
  #   A csv with the data of the appropriate data.table
  #
  
  # Construct the output folder path
  output.folder.path <- paste("data-output", output.folder, sep = "/")
  
  # Change to the appropriate directory
  setwd(here(output.folder.path))
  
  # Substitute the dots to dashes
  dashed.file.name <- gsub("\\.", "-", file.name)
  
  # Add the csv extension to the name
  dashed.file.name <- paste0(dashed.file.name, ".csv")
  
  # Write the data.table to csv
  fwrite(data, file = dashed.file.name)
  
  # Inform the user
  cat(dashed.file.name, "saved under", here("data-output/intermediate-data"), "\n")
  
  # Change back to the src directory
  setwd(here("src"))
}

filter.out.reverse.and.contaminants <- function(analysis.data) {
  #
  # Removes Contamintants and Reverse flagged peptides from the evidence data
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
                                        `Protein IDs`,
                                        perl = TRUE) == FALSE)
  
  # Clean them from the Protein IDs starting with REV__ regarding the Reverse sequences
  data.no.contaminants.no.reverse <- subset(data.no.contaminants,
                                            grepl("^REV__",
                                            `Protein IDs`,
                                            perl = TRUE) == FALSE)
  
  return (data.no.contaminants.no.reverse)
}

use.peptides.median <- function (intensities, minimum.detections = 2) {
  #
  # Get the median intensity of each peptide
  #
  # Args:
  #   intensities:        The list of the peptide's intensities between the technical replicates and fractions
  #   minimum.detections: Default is 2. The minimum number of detections to assume that the peptide was detected
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
  if (length(intensities) < minimum.detections) {
    return (list(NA))
  } else {
    return (list(median(intensities)))  
  }
  
}

do.vsn.normalization <- function(filtered.data, conditions.to.compare, minimum.detections = 2, do.norm = TRUE) {
  #
  # Does variance stabilizing normalization (VSN) on the peptides intensities
  # 
  # Args:
  #   filtered.data:          The filtered data data.table
  #   conditions.to.compare:  The conditions to compare
  #   minimum.detections:     Default is 2. A strictness parameter regarding how many times is a peptide measured
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
  setnames(data, "description", "Biological Replicate")
  
  # Reorder the data.table
  setkey(data, `Biological Replicate`,`Protein IDs`,`Unique Sequence ID` )
  
  # Iterate over conditions
  for (condition in conditions.to.compare)  {
    
    # Replace the NaN values with NA
    set(data, which(is.nan(data[, eval(condition)])), condition, NA)
    
    # Transfrorm the data from long format to wide where we have the peptides on the y axis and the biological
    # replicates on the x axis
    wide.data <- dcast(data,
                   `Protein IDs` + `Unique Sequence ID` ~ `Biological Replicate`,
                   value.var = condition,
                   fun.aggregate = use.peptides.median)
    
    # Get only the part of the data.table with the intensities of the peptides
    vsn.matrix <- as.matrix(wide.data[, .SD, .SDcols = -c(1, 2)])
    
    # Do the normalization and suppress the warnings regardind the NA removal the the vsn package throws by default
    if (do.norm == TRUE) {
      suppressWarnings(vsn.matrix.normalized <- justvsn(vsn.matrix))  
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
                                   by=c("Protein IDs", "Unique Sequence ID"))
    }
  }
  
  return (vsn.normalized.data)
}

do.peptide.intensities.plots <- function(intensities.before.normalization, intensities.after.normalization, plots.format = 5) {
  #
  # Makes 4 boxplot, before and after the normalization, with and without titles
  #
  # Args:
  #   intensities.before.normalization: The data.table of the non-normalized peptide intensities
  #   intensities.after.normalization:  The data.table of the normalized peptide intensities
  #   plots.format:           Default is 5 (jpeg format). A numeric value indicating the format of the plots. 
  #                           The numbers correspont to:  1     2     3     4      5      6     7     8     9
  #                                                     "eps" "ps"  "tex" "pdf" "jpeg" "tiff" "png" "bmp" "svg"
  #
  # Returns:
  #   The boxplots
  #
  
  # Copy the data.tables
  raw.data <- copy(intensities.before.normalization)
  normalized.data <- copy(intensities.after.normalization)
  
  # Remove the 'Protein IDs' columns
  raw.data[, `Protein IDs`:= NULL]
  normalized.data[, `Protein IDs`:= NULL]
  
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
  
  # Change to the limma directory
  setwd(here("data-output/plots"))
  
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

do.LCMD.imputation <- function(vsn.normalized.data) {
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
  
  # Make a matrix with only the intensities
  data.to.impute <- as.matrix(vsn.normalized.data[, -c(1,2)])
  
  # Do the imputation
  imputed.matrix <- impute.QRILC(data.to.impute)[[1]]
  
  # And update the initial data.table
  for (index in 1:ncol(imputed.matrix)) {
    imputed.data[, index + 2] <- imputed.matrix[, index]  
  }
  
  return (imputed.data)  
}

do.peptides.aggregation <- function(imputed.data) {
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
  
  # Do the aggregation of the peptides' intensity that belong to the same protein
  aggregated.data <- aggregated.data[, lapply(.SD, sum), .SDcols=!"Unique Sequence ID", by="Protein IDs"]
  
  # And rename the columns
  old.column.names <- colnames(aggregated.data)[-c(1)]
  
  # In order to  know that
  new.column.names <- paste("Protein Abundance", old.column.names)
  
  # The columns now hold the protein abundances
  colnames(aggregated.data)[-c(1)] <- new.column.names
  
  return (aggregated.data)
}

do.QQ.plots <- function(aggregated.data, conditions.to.compare, plots.format = 5) {
  #
  # Produces 2 QQ plots, one with tile and one with no title
  #
  # Args:
  #   aggregated.data:        The aggregated.data data.table
  #   conditions.to.compare:  The condition for the descriptive plot
  #   plots.format:           Default is 5 (jpeg format). A numeric value indicating the format of the plots. 
  #                           The numbers correspont to:  1     2     3     4      5      6     7     8     9
  #                                                     "eps" "ps"  "tex" "pdf" "jpeg" "tiff" "png" "bmp" "svg"
  # Returns:
  #   The saved plots inside the data-out/plots folder
  #
  
  # Copy the data.table
  data <- copy(aggregated.data)
  
  # Remove the protein ids column
  data[, "Protein IDs":= NULL]
  
  # Get the conditions' names
  condition.A <- conditions.to.compare[[1]]
  condition.B <- conditions.to.compare[[2]]
  
  
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
  
  # Get the data for condition A
  condition.B.data <- as.double(as.matrix(data[, ((condition.A.number.of.biological * condition.A.number.of.technical) + 1): ncol(data)]))
  
  # Tag them
  condition.B.data <- cbind(condition.B.data, 2)
  
  # Wrap them up in a data.frame
  qq.plot.data <-  as.data.frame(rbind(condition.A.data, condition.B.data))  
  
  # Fix the column names
  colnames(qq.plot.data) <- c("Abundancies", "Condition")
  
  # Make a QQ plot tha does not contain any title or description
  qq.plot.no.title <- ggplot(data       = qq.plot.data,
                            aes(sample  = Abundancies,
                                group   = Condition,
                                color   = as.factor(Condition))) +
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
                           y = "Sample QUantiles") +
                      scale_color_manual(values = colorblind.palette,
                                         name = "Conditions",
                                         labels = conditions.to.compare)
  
  # Now make a title for the descriptive QQ plot
  qq.plot.title <- paste0( "Q-Q plot of conditions ",
                           conditions.to.compare[[1]],
                           " and ",
                           conditions.to.compare[[2]],
                           "protein abundancies")
  
  # In case that the description is longer than 70 characters, break the title into multiple lines
  qq.plot.title <- paste0(strwrap(qq.plot.title, 70), collapse = "\n")      
  
  # Make a QQ plot tha contains title or description
  qq.plot.with.title <- ggplot(data       = qq.plot.data,
                             aes(sample  = Abundancies,
                                 group   = Condition,
                                 color   = as.factor(Condition))) +
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
                           y = "Sample QUantiles",
                           title = qq.plot.title) +
                      scale_color_manual(values = colorblind.palette,
                                         name = "Conditions",
                                         labels = conditions.to.compare)
  
  # Change to the limma directory
  setwd(here("data-output/plots"))
  
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

do.value.ordered.ratio.plot <- function(limma.results, conditions.to.compare, plots.format = 5) {
  #
  # Does 2 value ordered ratio plots, one with title and one with not.
  #
  # Args:
  #   limma.results:          The data.frame with the limma results
  #   conditions.to.compare:  The condition for the descriptive plot
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
                            
  # Change to the limma directory
  setwd(here("data-output/plots"))
  
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

do.MA.plots <- function(limma.results, conditions.to.compare, plots.format = 5) {
  #
  # Produces 2 MA plots, one with tile and one with no title
  #
  # Args:
  #   limma.results:          The data.frame with the limma results
  #   conditions.to.compare:  The condition for the descriptive plot
  #   plots.format:           Default is 5 (jpeg format). A numeric value indicating the format of the plots. 
  #                           The numbers correspont to:  1     2     3     4      5      6     7     8     9
  #                                                     "eps" "ps"  "tex" "pdf" "jpeg" "tiff" "png" "bmp" "svg"
  # Returns:
  #   The saved plots inside the data-out/plots folder
  
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
  
  # Change to the limma directory
  setwd(here("data-output/plots"))
  
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

do.volcano.plots <- function(limma.results, conditions.to.compare, plots.format = 5, fold.change.cut.off = 1.5, FDR = 0.05 ) {
  #
  # Generates and saves the 2 volcano plots in the appropriate format, one plain and one with descriptive info
  #
  # Args:
  #   limma.results:          The data.frame with the limma results
  #   conditions.to.compare:  The condition for the descriptive plot
  #   plots.format:           Default is 5 (jpeg format). A numeric value indicating the format of the plots. 
  #                           The numbers correspont to:  1     2     3     4      5      6     7     8     9
  #                                                     "eps" "ps"  "tex" "pdf" "jpeg" "tiff" "png" "bmp" "svg"
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
  
  # Change directory to the limma-output folder, inside the data-output folder
  setwd(here("data-output/plots"))
  
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
                              error.correction.method = "B", fold.change.cut.off = 1.5, FDR = 0.05) {
  #
  # Does the limma analysis and returns a dataframe with the results of limme
  #
  # Args:
  #   aggregated.data:          The protein abundancis data.table
  #   conditions.to.compare:    The condition to compare
  #   experimental.metadata:    The metadata list of the experiment (condition >  number of biological replicates/
  #                                                                               number of technical replicates/
  #                                                                               number of fractions)
  #   fold.change.cut.off:      Default is 1.5. The least fold change  in order to consider a protein differentially expressed
  #   FDR:                      Default is 0.05. The least fold change  in order to consider a protein differentially expressed
  #   error.correction.method:  Default is "B". Corrects the Type 1 errors using the Bonferroni correction method or
  #                             the Benjamini-Hochberg method. Values can be "B" or "BH'
  #
  # Returns:
  #   A data.frame with the results of the limma analysis
  #
  
  # Copy the data
  limma.data <- copy(aggregated.data)
  
  # Get the protein names
  protein.ids <- limma.data[[1]]
  
  # And remove the protein names from the data.table
  limma.data[[1]] <- NULL
  
  # Now convert the data.table to matrix
  limma.data <- as.matrix(limma.data)
  
  # And add the protein names as rownames  
  rownames(limma.data) <- protein.ids
  
  # Then get the colnames
  limma.column.names <- colnames(limma.data)
  
  # Remove the protein abundance prefic and the technical replicate prefix
  limma.column.names <- gsub("^(Protein Abundance )|(T[[:digit:]*]$)", "", limma.column.names, perl = TRUE)
  
  #  Get the condition names
  condition.A <- conditions.to.compare[[1]]
  condition.B <- conditions.to.compare[[2]]
  
  # Get the number of biological replicates for each condition
  condition.A.number.of.biological <- experimental.metadata[[condition.A]]$number.of.biological.replicates
  condition.B.number.of.biological <- experimental.metadata[[condition.B]]$number.of.biological.replicates
  
  # Get the number of technical replicates for each condition
  condition.A.number.of.technical <- experimental.metadata[[condition.A]]$number.of.technical.replicates
  condition.B.number.of.technical <- experimental.metadata[[condition.B]]$number.of.technical.replicates
  
  # Construct the design matrix
  design <- cbind(Intercept = 1, condition.B = c( rep.int(0, condition.A.number.of.biological * condition.A.number.of.technical),
                                                  rep.int(1, condition.B.number.of.biological * condition.B.number.of.technical))
                  )
  # ,
  #                               technical.replicate = c(rep.int(1:condition.A.number.of.technical, condition.A.number.of.biological),
  #                                                       rep.int(1:condition.B.number.of.technical, condition.B.number.of.biological)))
  # 

  # Do the linear modelling
  limma.fit <- lmFit(limma.data, design) 
  
  # And calculates the statistics
  limma.fit <- treat(limma.fit, lfc = log2(fold.change.cut.off), trend = TRUE, robust = TRUE)
  
  # Make a data.frame with the statistics for each protein, sorted by p-value, either it is significant or not
  limma.results <- topTable(limma.fit, coef = 2, sort.by="P", n=Inf)
  
  # And a column with the protein names
  limma.results <- cbind(Protein=rownames(limma.results), limma.results)
  
  # And restore the rownames to numbers
  rownames(limma.results) <- 1:nrow(limma.results)
  
  # Depending on the the results conservativeness, use the appropriate method
  switch(error.correction.method,
         "B" =  {
                   limma.results$is.diffenetially.expressed =  (abs(limma.results$logFC) > fold.change.cut.off) &
                                                                limma.results$P.Value  < (0.05/length(protein.ids))
                },
         "BH" = {
                  limma.results$is.diffenetially.expressed =  (abs(limma.results$logFC) > fold.change.cut.off) &
                                                              limma.results$adj.P.Val  < FDR
                },
         {
           stop("Invalid limma method for error correction.\n")
         })
  
  
  return (limma.results)
 
}

tmp.evalutate.correctness.t.test <- function(aggregated.data, proteins=c("PA3479", "PA5346", "PA0176", "PA4625", "PA4624", "PA3724", "PA3385", "PA5272", "PA3217", "PA3544", "PA5060")) {
  results <- lapply(proteins, function(proteinID) {
    test.data <- aggregated.data[grep(proteinID, aggregated.data[[1]]), ]
    group1 <- unlist(test.data[1,2:7])
    group2 <- unlist(test.data[1,8:13])
    return (t.test(group1, group2)$"p.value")
  })
  results <- unlist(results)
  significant.proteins <- proteins[results < 0.05]
  number.of.significant <- length(significant.proteins)
  cat(toString(number.of.significant),"proteins where significant:", significant.proteins)
}

tmp.evalutate.correctness.limma <- function(limma.results, proteins=c("PA3479", "PA5346", "PA0176", "PA4625", "PA4624", "PA3724", "PA3385", "PA5272", "PA3217", "PA3544", "PA5060")) {
  results <- c()
  for (protein in proteins) {
    id <- grep(protein, limma.results$Protein)
    results <- c(results, limma.results$is.diffenetially.expressed[id])
  }
  cat(length(proteins[results]), "proteins where significant:",proteins[results],"\n")
}