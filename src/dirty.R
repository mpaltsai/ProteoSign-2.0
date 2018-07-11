# All functions exist in this file. If they are too many, they can be separeted in
# several functions_XXX.R files.

# library(log4r)
# library(gdata)
# library(pryr)

# Confirmed linraries
# library(data.table)
# library(feather)

# Moved to functions_build.R
# replicates.per.condition <- function(biological.replicates, technical.replicates, conditions) {
#   
#   replicates.per.condition <- list()
#   unique.conditions <- unique(conditions)
#   
#   for (condition in unique.conditions) {
#     raw.files.per.condition <- which(condition == conditions)
#     biological.replicates.per.condition <- biological.replicates[raw.files.per.condition]
#     technical.replicates.per.condition <- technical.replicates[raw.files.per.condition]
#     new.element <- list("biological" = biological.replicates.per.condition,
#                         "technical"= technical.replicates.per.condition)
#     if (length(replicates.per.condition) == 0 ) {
#       replicates.per.condition <- new.element
#     } else {
#       replicates.per.condition <- list(replicates.per.condition, new.element)
#     }
#   }
#   names(replicates.per.condition) <-  unique.conditions
#   return (replicates.per.condition)
# }

# Moved to functions_build.R
# replicates.status.per.condition <- function(replicates.per.condition) {
#   
#   replicates.status.per.condition <- list()
#   for (condition in replicates.per.condition) {
#     biological.replicates <- condition[[1]]
#     technical.replicates <- condition[[2]]
#     
#     new.status <- check.replicates(biological.replicates,
#                                    technical.replicates)
#     
#     if (length(replicates.status.per.condition) == 0) {
#       replicates.status.per.condition <- new.status
#     } else {
#       replicates.status.per.condition <- list(replicates.status.per.condition, new.status)
#     }
#   }
#   names(replicates.status.per.condition) <- names(replicates.per.condition)
#   return (replicates.status.per.condition)
# }

# Moved to functions_build.R
# check.replicates <- function(biological.replicates, technical.replicates) {
#   
#   unique.biological.replicates <-  unique(biological.replicates)
#   unique.technical.replicates <-  sort(unique(technical.replicates))
#   
#   expected.biological.replicates <- c(1:length(unique.biological.replicates))
#   expected.technical.replicates <- c(1:length(unique.technical.replicates))
#   
#   biological.replicates.are.ok <- all(unique.biological.replicates == expected.biological.replicates)
#   technical.replicates.are.ok <- all(unique.technical.replicates == expected.technical.replicates)
#   
#   replicates.status <- list( "biological" = biological.replicates.are.ok,
#                              "technical" = technical.replicates.are.ok)
#   
#   return (replicates.status)
# }
# 
# Moved to functions_build.R
# fix.replicates <- function(biological.replicates, technical.replicates) {
#   
#   unique.biological <- unique(biological.replicates)
#   biological.sample.id <- 1
#   for (biorep in unique.biological) {
#     
#     same.samples.biological <- which(biological.replicates == biorep)
#     biological.replicates[same.samples.biological] <- biological.sample.id
#     biological.sample.id <- biological.sample.id + 1
#     unique.technical <- unique(technical.replicates[same.samples.biological])
#     technical.sample.id <- 1
#     
#     for (techrep in unique.technical) {
#       same.samples.technical <- which(technical.replicates == techrep & biological.replicates == biological.sample.id - 1)
#       technical.replicates[same.samples.technical] <- technical.sample.id
#       technical.sample.id <- technical.sample.id +1
#     }
#   }
#   fixed.replicates <- list(biological.replicates, technical.replicates)
#   names(fixed.replicates) <- c("biological", "technical")
#   return (fixed.replicates)
# }
# 
# fix.replicates.per.condition <- function (replicates.per.condition, replicates.status.per.condition) {
#   
#   biological.replicate.id <- 1
#   
#   for (index in 1:length(replicates.status.per.condition)) {
#     condition <- names(replicates.status.per.condition)[index]
#     biological.status <- replicates.status.per.condition[[index]]$biological
#     technical.status  <- replicates.status.per.condition[[index]]$technical
#     
#     case <- (biological.status * 2  + technical.status * 1) + 1
#     switch(case,
#            cat("Condition", condition, "has bad biological and technical replicates.\n") ,
#            cat("Condition", condition, "has bad biological replicates.\n") ,
#            cat("Condition", condition, "has bad technical replicates.\n") ,
#            cat("Condition", condition, "has good replicates.\n"))
#     if (case != 4) {
#       cat("Fixing...\n")
#       biological.replicates <- replicates.per.condition[[index]]$biological
#       technical.replicates <- replicates.per.condition[[index]]$technical
#       fixed.replicates <- fix.replicates(biological.replicates, technical.replicates)
#       replicates.per.condition[[index]] <- fixed.replicates
#     }
#   }
#   return (replicates.per.condition)
# }
# Moved to functions_build.R
# restore.replicates <- function(replicates.per.condition) {
#   biological <- c()
#   technical <- c()
#   for (rep in replicates.per.condition) {
#     biological <- c(biological, rep$biological)
#     technical <- c(technical, rep$technical)
#   }
#   restored.replicates <- list("biological" = biological, "technical" = technical)
#   return (restored.replicates)
# }

# Moved to functions_build.R
# check.replicates.number <- function(replicate.multiplexing.is.used, replicates) {
#   number.of.replicates <- length(unique(replicates)) 
#   if (replicate.multiplexing.is.used == FALSE && number.of.replicates <= 1) {
#     return (FALSE)
#   } else {
#     return (TRUE)
#   }
# }
# Moved to functions_build.R
# make.experimental.description <- function(experimental.setup.id, biological.replicates, technical.replicates, fractions) {
#   experimental.description <- c()
#   experimental.description.biological <- paste0("b", biological.replicates)
#   experimental.description.technical <- paste0("t", technical.replicates)
#   experimental.description.fractions <- paste0("f", fractions)
#   switch(experimental.setup.id,
#          {
#            cat("We have bioreps.\n")
#            experimental.description <- experimental.description.biological
#          },
#          {
#            cat("No bioreps, no fractions? Really? \n")
#            cat("Invalid experimental description id in make.experimental.description function.\n")  
#            return (FALSE)
#          },
#          {
#            cat("We have bioreps and techreps.\n")
#            experimental.description <- paste0(experimental.description.biological,
#                                               experimental.description.technical)
#          },
#          {
#            cat("We have bioreps and fractions.\n")
#            experimental.description <- paste0(experimental.description.biological,
#                                               experimental.description.fractions)
#          },
#          {
#            cat("No bioreps? Really?\n")
#            cat("Invalid experimental description id in make.experimental.description function.\n")  
#            return (FALSE)
#          },
#          {
#            cat("We have bioreps, techreps and fractions.\n")
#            experimental.description <- paste0(experimental.description.biological,
#                                               experimental.description.technical,
#                                               experimental.description.fractions)
#          })
#   cat(experimental.description,"\n")
#   return (experimental.description)
# }
#  No need for this, used tr unix instead
# clean.file.from.quotes <- function(file.name) {
#   file.first.line <- readLines(file.name, n=1)
#   have.quotes <- grepl("\"", file.first.line)
#   cat('Reading', file.name, '...\n')
#   # Moved to load_data.R
#   # input.data <- fread(file.name, integer64 = "numeric")
#   if (have.quotes == TRUE) {
#     cat("Removing double quotes from input data file", file.name, "...\n")
#     
#     print(typeof(input.data))
#     # Too slow, do i need it?
#     #cleaned.data <- gsub("\"", "", input.data)
#     
#     # evidence.file_cleaned <- file(evidence.file, open="w")
#     # writeLines(tmpdata, con=evidence.file_cleaned)
#     # close(evidence.file_cleaned)
#   } else {
#     cat("File is empty from quotes!\n")
#     cleaned.data <- input.data
#   }
#   return (cleaned.data)
# }
# Moved to functions_build.R
# find.proteome.discoverer.protein.column <- function(evidence.data) {
#   protein.groups.column <- FALSE
#   if ('Protein.Group.Accessions' %in% colnames(evidence.data) == TRUE) {
#     protein.groups.column<-'Protein.Group.Accessions'
#   } else if ('Protein.Accessions' %in% colnames(evidence.data) == TRUE) {
#     protein.groups.column<-'Protein.Accessions'
#   } else {
#     cat("Error User: The dataset does not contain the columns 'Protein Group Accessions' or 'Protein Accessions'\n")
#   }
#   
#   return (protein.groups.column)
# }

correct.maxquant.protein.groups <- function(protein.groups.data, evidence.data) {
  ## For MaxQuant correct protein groups in the evidence file using the
  # protein groups file.
  # protein.groups.data <- fread(protein.groups.file, integer64 = "numeric")
  # print(length(protein.groups.data$`Peptide counts (all)`))
  # If there isn't a Protein.Names or Protein.names column
  # (depends on MQ version), create one from the Fasta Headers column
  protein.groups.column.names <- colnames(protein.groups.data)
  evidence.column.names <- colnames(evidence.data)
  
  number.of.Protein.Names.columns <- length(which(protein.groups.column.names == 'Protein.Names'))
  number.of.Protein.names.columns <- length(which(protein.groups.column.names == 'Protein.names'))
  
  
  if(number.of.Protein.Names.columns == 0 && number.of.Protein.names.columns == 0){
    protein.groups.data$`Protein Names` <- gsub("(^>[[:alnum:]]+[[:punct:][:blank:]])|(;>.*)|(;[[:alnum:]]+)",
                                                "",
                                                protein.groups.data$`Fasta headers`,
                                                perl = TRUE)
  }
  if (number.of.Protein.names.columns > 0)
  {
    setnames(protein.groups.data, "Protein names", "Protein.Names")
  }
  
  # Construct a table, mapping the correct protein groups IDs (and the corresponding proteins names) to the evidence IDs
  #First check if there are any blank lines and remove them:
  setkey(protein.groups.data, "Evidence IDs")
  
  protein.groups.data <- protein.groups.data[!""]
  
  if (!"Protein IDs" %in% protein.groups.column.names & "Peptide IDs" %in% protein.groups.column.names)
  {
    setnames(protein.groups.data, "Peptide IDs", "Protein IDs")
  }
  
  protein.groups.subset <- protein.groups.data[, .SD, .SDcols = c("Protein IDs", "Protein Names", "Evidence IDs")]
  
  protein.groups.subset.multiline.evidence <- protein.groups.subset[ ,list(`Evidence IDs` = unlist(strsplit(`Evidence IDs`,
                                                                                                            ";"))),
                                                                     by= list(`Protein IDs`, `Protein Names`)]
  
  setnames(protein.groups.subset.multiline.evidence,
           "Evidence IDs",
           "id")
  
  class(protein.groups.subset.multiline.evidence$id) <- 'integer'
  setkey(protein.groups.subset.multiline.evidence, id)
  
  # Get the evidence records
  setkey(evidence.data, id)
  # Remove the incorrect Protein.IDs column
  column.delete.case <- (any(grepl("Protein IDs", evidence.column.names, perl = TRUE)) * 1) +
    (any(grepl("Protein Names", evidence.column.names, perl = TRUE)) * 2) + 1
  switch(column.delete.case,
         cat("No need for column deletion...\n"),
         evidence.data[, "Protein IDs" := NULL],
         evidence.data[, "Protein Names" := NULL],
         evidence.data[, c("Protein IDs", "Protein Names") := NULL])
  
  evidence.column.names <- colnames(evidence.data)
  
  # Inner join the mapping table with the evidence table and return the data frame that we ought to have in the first place
  
  merged.evidence.table <- merge(protein.groups.subset.multiline.evidence, evidence.data, by = "id")
  
  merged.evidence.table.columns <- colnames(merged.evidence.table)
  
  if("Contaminant" %in% merged.evidence.table.columns){
    setkey(merged.evidence.table, Contaminant)
    merged.evidence.table <- merged.evidence.table[!"+"]
  }
  
  if("Reverse" %in% merged.evidence.table.columns){
    setkey(merged.evidence.table, Reverse)
    merged.evidence.table <- merged.evidence.table[!"+"]
  }
  
  setkey(merged.evidence.table, id)
  
  corrected.files <- list("protein.groups.file"= protein.groups.data, "evidence.file" = merged.evidence.table)
  return (corrected.files)
}

read.protein.groups.data_v3 <- function(protein.groups.file,
                                        evidence.file,
                                        time.point,
                                        keep.evidence.ids=F,
                                        is.proteome.discoverer.data) {
  # Moved to load_data.R
  # cat("Reading data file...\n")
  # protein.groups.column <- ""
  # 
  # evidence.data <- clean.file.from.quotes(evidence.file)
  # protein.groups.data <- NULL
  # if (is.proteome.discoverer.data == TRUE) {
  #   cat('Proteome Discoverer data...\n')
  #   protein.groups.column <- find.proteome.discoverer.protein.column(evidence.data)
  # } else {
  #   cat('Maxquant data...\n')
  #   cleaned.protein.groups.data <- clean.file.from.quotes(protein.groups.file)
  #   protein.groups.data <- cleaned.protein.groups.data
  #   protein.groups.column<-'^Proteins$'
  # }
  # 
  # protein.groups.column.position <- which(colnames(evidence.data) == protein.groups.column)
  # 
  # colnames(evidence.data)[protein.groups.column.position] <- "Protein IDs"
  
  if (is.proteome.discoverer.data == FALSE) {
    corrected.files <- correct.maxquant.protein.groups(protein.groups.data, evidence.data)
    
    protein.groups.data <- corrected.files$protein.groups.file
    evidence.data <- corrected.files$evidence.file
  }
  
  # In the case of Isobaric labeling we should reformat the table before
  # proceeding, afterwards we will treat the data as
  # if they were label-free data
  
  analysis.parameters <- list("dataset.source" = "Maxquant",
                              "is.isobaric" = FALSE,
                              "is.label.free" = FALSE,
                              "allow.label.swap" = FALSE,
                              "allow.label.rename" = FALSE,
                              "replicate.multiplexing" = FALSE)
  
  if (is.isobaric == TRUE) {
    if(is.proteome.discoverer.data == FALSE) {
      evidence$Intensity <- NULL
      varcolnames <- grep("^Reporter.intensity.[[:digit:]]", colnames(evidence.data), value = TRUE)
      evidence <- reshape(evidence, varying = varcolnames, v.names = "Intensity", timevar = "Labeling.State", times = varcolnames, direction = "long", new.row.names=sequence(prod(length(varcolnames), nrow(evidence.data))))
      conditions.labels << -sub("^X", "Reporter.intensity.", conditions.labels)
      if (AllowLabelRename == T)
      {
        Rename_Array$old_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1", Rename_Array$old_label)
        Rename_Array$new_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1", Rename_Array$new_label)
      }
      if (AllowLS == T)
      {
        Ls_array$first_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1",  Ls_array$first_label)
        Ls_array$second_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1",  Ls_array$second_label)
      }
      is.label.free<-T;
      filterL_lbl <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1",  filterL_lbl)
      if(replicate.multiplexing.is.used){
        RMtagsdata$name <- sub("^([[:digit:]]*[[:alpha:]]?)$", "Reporter.intensity.\\1",  RMtagsdata$name)
      }
    }
    else{
      evidence$Intensity <- NULL
      if (any(grepl("^Abundance..", colnames(evidence.data)))){
        colnames(evidence.data) <- sub("^Abundance..", "X", colnames(evidence.data))
      }
      varcolnames <- grep("^X[[:digit:]]*[[:alpha:]]?$", colnames(evidence.data), value = TRUE)
      evidence <- reshape(evidence, varying = varcolnames, v.names = "Intensity", timevar = "Modifications", times = varcolnames, direction = "long", new.row.names=sequence(prod(length(varcolnames), nrow(evidence.data))))
      is.label.free<-T;
      if (AllowLabelRename == T)
      {
        Rename_Array$old_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", Rename_Array$old_label)
        Rename_Array$new_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", Rename_Array$new_label)
      }
      if (AllowLS == T)
      {
        Ls_array$first_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", Ls_array$first_label)
        Ls_array$second_label <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", Ls_array$second_label)
      }
      filterL_lbl <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", filterL_lbl)
      if(replicate.multiplexing.is.used){
        RMtagsdata$name <- sub("^([[:digit:]]*[[:alpha:]]?)$", "X\\1", RMtagsdata$name)
      }
    }
  }
  # #Erase all rows in rename_array that try to rename a label to an already existing
  # mi<-which(Rename_Array$new_label %in% conditions.labels)
  # Rename_Array <- Rename_Array[-mi,]
  # Rename_Array <<- Rename_Array
  
  #levellog(paste0("read.protein.groups.data_v3: Identified proteins: ",length(unique(evidence$Protein.IDs))," (",time.point,")"))
  
  old.evidence.data <- evidence.data
  setkey(old.evidence.data, "Protein IDs")
  evidence.data <- old.evidence.data[!" "]
  
  cat("read.protein.groups.data_v3: Discarded PSM records due to unassigned protein group: ", nrow(old.evidence.data) - nrow(evidence.data))
  
  ## Make Protein.IDs human-readable
  if(is.proteome.discoverer.data == TRUE){
    protein.groups.column <- "Protein Descriptions"
    if (!"Protein Descriptions" %in% colnames(evidence.data)) {
      evidence.data[, `Protein Descriptions` := ""]
    }
  } else {
    protein.groups.column <- "Protein Names"
  }
  #### Shady part, maybe errors, bugs, accidental demon summonning
  #### Maybe duplicate of Maxquant block
  
  tmp.table <- evidence.data[, .SD, .SDcols = c("Protein IDs", protein.groups.column)]
  tmp.table$id <- 1:nrow(tmp.table)
  setkey(tmp.table, "Protein IDs")
  
  # Generate data.table with unique Protein.IDs
  tmp.table2<-tmp.table[,.("Protein ID quantity"=.N), by="Protein IDs"]
  setkey(tmp.table2, "Protein IDs")
  
  # Make a new protein description column in other data.table
  tmp.table[, tmp.description := do.call(paste, c(.SD)), .SDcols = c("Protein IDs", "Protein Names")]
  
  # set the Protein.IDs in the original data frame
  evidence.shady <- tmp.table2[tmp.table][order(id), tmp.description]
  evidence.data$`Protein IDs` <- tmp.table2[tmp.table][order(id), tmp.description]
  #merged.evidence.table <- merge(protein.groups.subset.multiline.evidence, evidence.data, by = "id")
  
  #### End of shady part, satan has left the building
  
  ## Assign defined labels (conditions), one for each PSM record
  
  if (is.proteome.discoverer.data == TRUE) {
    raw.file.column <- "Spectrum File" 
  } else {
    if ( "Raw File" %in% evidence.column.names) {
      raw.file.column <- 'Raw File'
    } else {
      raw.file.column <- 'Raw file'
    }
  }
  
  
  if(is.label.free == TRUE){
    condition.specification.column <- raw.file.column
    
    if(is.isobaric == TRUE){
      if (is.proteome.discoverer.data == TRUE) {
        condition.specification.column <- "Modifications"
      } else {
        condition.specification.column <- "Labeling State"
      }
    }
  } else {
    if (is.proteome.discoverer.data == TRUE) {
      condition.specification.column <- "Modifications"
    } else {
      condition.specification.column <- "Labeling State"
    }
  }
  
  if(replicate.multiplexing.is.used == TRUE){
    levellog("read.protein.groups.data_v3: Transforming data for Replication Multiplexing ...")
    #when RM is chosen, in this line evidence
    #has two main columns, one called Labeling.State or Modifications
    #that tells us what tag it was tagged and one called 
    #Spectrum File or Raw File that says from which raw file did 
    #the psm come from
    #in case of RM breps treps and conditions may come from 
    # either raw files or tags but in this data format creating
    #two new columns describing the structure correctly is not difficult
    #first create the column for conditions
    RMrawfilesdata <- RMrawfilesdata[!RMrawfilesdata$used == 'false',]
    RMtagsdata <- RMtagsdata[!RMtagsdata$used == 'false',]
    if(RMconditionsinrawfiles)
    {
      evidence <- merge(evidence, RMrawfilesdata[, c('name', 'cond')], by.x = raw.file.column, by.y = 'name')
    } else {
      evidence <- merge(evidence, RMtagsdata[, c('name', 'cond')], by.x = condition.specification.column, by.y = 'name')
    }
    #the conditions.labels array is already set from the front end
    #Now we will initialize the new_raw_file column that will contain pseudo-raw files describing our bioreps, techreps and fracs
    colnames(RMrawfilesdata)[3:5] <- c('new_brep', 'new_trep', 'new_frac')
    colnames(RMtagsdata)[3:5] <- c('new_brep', 'new_trep', 'new_frac')
    if(RMbrepsinrawfiles)
    {
      evidence <- merge(evidence, RMrawfilesdata[, c('name', 'new_brep')], by.x = raw.file.column, by.y = 'name')
    } else {
      evidence <- merge(evidence, RMtagsdata[, c('name', 'new_brep')], by.x = condition.specification.column, by.y = 'name')
    }
    #do the same thing for treps
    if(RMtrepsinrawfiles)
    {
      evidence <- merge(evidence, RMrawfilesdata[, c('name', 'new_trep')], by.x = raw.file.column, by.y = 'name')
    } else {
      evidence <- merge(evidence, RMtagsdata[, c('name', 'new_trep')], by.x = condition.specification.column, by.y = 'name')
    }
    if(RMbrepsinrawfiles | RMtrepsinrawfiles) {
      #do the same thing for fracs
      evidence <- merge(evidence, RMrawfilesdata[, c('name', 'new_frac')], by.x = raw.file.column, by.y = 'name')
    }
    if(!RMbrepsinrawfiles & !RMtrepsinrawfiles) {
      evidence$new_raw_file <- paste0('b', evidence$new_brep, 't', evidence$new_trep)
    } else {
      evidence$new_raw_file <- paste0('b', evidence$new_brep, 't', evidence$new_trep, 'f', evidence$new_frac)
    }
    #Now lets refresh the rep_structure array the pseudo raw files we created are descriptive and contain the the breps treps and fracs
    pseudo_raw_files <- unique(evidence$new_raw_file)
    new_rep_structure <- rep_structure
    new_rep_structure <- new_rep_structure[0,]
    for (i in 1:length(pseudo_raw_files)) {
      levels(new_rep_structure$raw_file) <- c(levels(new_rep_structure$raw_file), pseudo_raw_files[i])
      new_rep_structure[i,] <- c(as.character(pseudo_raw_files[i]), NA, NA, NA, NA)
    }
    colnames(new_rep_structure) <- c('raw_file','biorep','techrep','fraction', 'rep_desc')
    if(RMbrepsinrawfiles | RMtrepsinrawfiles) {
      for(i in 1:nrow(new_rep_structure))
      {
        new_rep_structure[i, c('raw_file', 'biorep','techrep','fraction')] <- str_match_all(new_rep_structure$raw_file[i], "b(.*?)t(.*?)f(.*)")[[1]][1,]
      }
    } else{
      for(i in 1:nrow(new_rep_structure))
      {
        new_rep_structure[i, c('raw_file', 'biorep','techrep')] <- str_match_all(new_rep_structure$raw_file[i], "b(.*?)t(.*)")[[1]][1,]
      }
      new_rep_structure[, 'fraction'] <- "1"
    }
    
    if(length(unique(new_rep_structure$techrep)) > 1){
      if(length(unique(new_rep_structure$fraction)) > 1){
        # we have techreps and fractions
        new_rep_structure$rep_desc<-paste(paste(paste('b',new_rep_structure$biorep,sep=''),'t',new_rep_structure$techrep,sep=''),'f',new_rep_structure$fraction,sep='')
      }else{
        #we have bioreps and techreps
        new_rep_structure$rep_desc<-paste(paste('b',new_rep_structure$biorep,sep=''),'t',new_rep_structure$techrep,sep='')
      }
    }else{
      if(length(unique(new_rep_structure$fraction)) > 1){
        # we have fractions but not techreps
        new_rep_structure$rep_desc<-paste(paste(paste('b',new_rep_structure$biorep,sep=''),'t',new_rep_structure$techrep,sep=''),'f',new_rep_structure$fraction,sep='')
      }else{
        # we just have bioreps
        new_rep_structure$rep_desc<-paste(paste('b',new_rep_structure$biorep,sep=''),'t',new_rep_structure$techrep,sep='')
      }
    }
    .GlobalEnv[["rep_structure"]]<-new_rep_structure
    #now erase the old columns for raw files and conditions and replace with the new ones
    evidence[, c(raw.file.column, condition.specification.column)] <- list(NULL)
    colnames(evidence.data)[colnames(evidence.data) == "new_raw_file"] <- raw.file.column
    colnames(evidence.data)[colnames(evidence.data) == "cond"] <- condition.specification.column
    if(length(unique(rep_structure$biorep)) == 1){
      levellog("Error User: Cannot accept dataset with just one biological replicate. Aborting ...")
      return(F)
    }
  }
  
  cat("read.protein.groups.data_v3: Assigning labels ...")
  # levellog("",change=1)
  evidence.data[, `Condition Label`:= ""]
  
  ### What is this?
  background_species_lbl<-NA
  
  evidence.labeling
  for(i in 1:length(conditions.labels)){
    if (is.proteome.discoverer.data == TRUE) {
      if (is.label.free == TRUE) {
        if (is.isobaric == FALSE) {
          mi<-which(grepl(conditions.labels[i], LFQ_conds[, "condition"]))
          mi2<-which(grepl(paste(LFQ_conds[mi,]$raw_file, collapse="|"), evidence[, condition.specification.column]))
          evidence[mi2,]$label_<-conditions.labels[i]
        } else {
          evidence$label_<-evidence$Modifications
        }
      } else {
        evidence$label_<-evidence$Quan.Channel
      }
    } else {
      if (is.label.free == TRUE){
        if (is.isobaric == FALSE){
          mi<-which(grepl(conditions.labels[i], LFQ_conds[, "condition"]))
          mi2<-which(grepl(paste(LFQ_conds[mi,]$raw_file, collapse="|"), evidence[, condition.specification.column]))
          evidence[mi2,]$label_<-conditions.labels[i]
        } else {
          evidence$label_<-evidence$Labeling.State
        }
      } else {
        # MQ nomenclature for labels: 0 the first label, 1 the second etc ...
        # Since the user might opted for excluding some labels
        # find the index of the included labels and parse them to the label_ column
        mi0<-which(All_MQ_Labels == conditions.labels[i])
        mi<-which(grepl((mi0[1]-1), evidence[, condition.specification.column]))
        evidence[mi,]$label_<-conditions.labels[i]
      }
    }
    levellog(paste0("read.protein.groups.data_v3: Assigned label '", conditions.labels[i],"'."))
  }
  #Rename any labels if necessary
  if (AllowLabelRename == T)
  {
    if (length(Rename_Array$old_label) != 0)
    {
      for(i in 1:length(Rename_Array$old_label))
      {
        if(Rename_Array$old_label[i] != Rename_Array$new_label[i])
        {
          #The case where old label = new label is common since if the user did not ask for a label rename (merge) for a label
          #this label is sent to be renamed in R to itself
          #in any other case rename the labels that are merged to the same label so that they become indistinguishable
          #and refresh the conditions labels by erasing the old label and adding the new if necessary
          mi<-which(evidence$label_ == Rename_Array$old_label[i])
          if (is.label.free == FALSE & is.isobaric == FALSE)
          {
            #in case of precursor ion data add you need to create a new column with the new cond name containing
            #the intensities of the old labels to the line where the label is found:
            prefix<-NA
            if (is.proteome.discoverer.data)
            {
              prefix<-""
            }
            else
            {
              prefix<-"Intensity."
            }
            evidence[mi, paste0(prefix, Rename_Array$new_label[i])] <- evidence[mi,  paste0(prefix, Rename_Array$old_label[i])]
          }
          evidence$label_[mi] <- Rename_Array$new_label[i]
          mi<-which(conditions.labels == Rename_Array$old_label[i])
          if (length(mi) > 0)
          {
            conditions.labels <- conditions.labels[-mi]
          }
          mi<-which(conditions.labels == Rename_Array$new_label[i])
          if (length(mi) == 0)
          {
            conditions.labels <- c(conditions.labels, Rename_Array$new_label[i])
          }
        }
      } 
      conditions.labels <<- conditions.labels
      nConditions<<-length(conditions.labels)
    }
    #No condition can be named "N" in ProteoSign so take care of this scarce situation:
    for(i in 1:length(conditions.labels)){
      if(conditions.labels[i] == "N")
      {
        conditions.labels[i] <- "condN"
        conditions.labels <<- conditions.labels
        evidence$label_ <- sub("^N$", "condN", evidence$label_)
      }
    }
  }
  
  #Rename again the labels in case of label swapping
  if (AllowLS == T)
  {
    for(i in 1:length(Ls_array$first_label))
    {
      mi1<-which(evidence[, raw.file.column] == Ls_array$selected_raw_file[i] & evidence$label_ == Ls_array$first_label[i])
      mi2<-which(evidence[, raw.file.column] == Ls_array$selected_raw_file[i] & evidence$label_ == Ls_array$second_label[i])
      evidence$label_[mi1] <- as.character(Ls_array$second_label[i])
      evidence$label_[mi2] <- as.character(Ls_array$first_label[i])
      if (is.label.free == FALSE & is.isobaric == FALSE)
      {
        prefix<-NA
        if (is.proteome.discoverer.data)
        {
          prefix<-""
        }
        else
        {
          prefix<-"Intensity."
        }
        #in case of precursor ion we should also swap the values between the columns "Intensity" of the first and second label
        if (length(mi1) > 0)
        {
          evidence[mi1, "temp_LS_Intensities"] <- evidence[mi1, paste0(prefix, Ls_array$second_label[i])]
          evidence[mi1, paste0(prefix, Ls_array$second_label[i])] <- evidence[mi1, paste0(prefix, Ls_array$first_label[i])]
          evidence[mi1, paste0(prefix, Ls_array$first_label[i])] <- evidence[mi1, "temp_LS_Intensities"]
        }
        if (length(mi2) > 0)
        {
          evidence[mi2, "temp_LS_Intensities"] <- evidence[mi2, paste0(prefix, Ls_array$second_label[i])]
          evidence[mi2, paste0(prefix, Ls_array$second_label[i])] <- evidence[mi2, paste0(prefix, Ls_array$first_label[i])]
          evidence[mi2, paste0(prefix, Ls_array$first_label[i])] <- evidence[mi2, "temp_LS_Intensities"]
        }
      }
    }
  }
  levellog("",change=-1)
  mi<-which(is.na(evidence$label_))
  if(is.na(background_species_lbl)){
    if(length(mi) > 0){
      evidence<-evidence[-mi,]
      levellog(paste("read.protein.groups.data_v3: Discarded PSM records due to unassigned label: ",length(mi),sep=""))
    }
  }else{
    evidence[mi,]$label_<-background_species_lbl
  }
  # Now add the experimental structure information
  evidence<-merge(evidence, .GlobalEnv[["rep_structure"]], by.x=c(raw.file.column), by.y=c('raw_file'))
  new_cond_labels <- NULL
  for (cond_i in conditions.labels)
  {
    if (!(cond_i %in% evidence$label_))
    {
      levellog(paste0("Warn User: ", cond_i, " label was not found in the selected raw files and was not used in comparisons!"))
      if(filterL_lbl == cond_i)
      {
        filterL<-F
        levellog("Warning!: the filter label was not found in the selected raw files. Filtering will not take place!")
      }
    }
    else
    {
      new_cond_labels <- c(new_cond_labels, cond_i)
    }
  }
  if(length(new_cond_labels)>1)
  {
    conditions.labels <- new_cond_labels
    conditions.labels <<- conditions.labels
    nConditions<<-length(conditions.labels)
  }else{
    levellog(paste0("Error User: Not enough labels left, aborting..."))
    return(F)
  }
  ## If we have fractionation, remake the rep_desc column and don't take into account the fraction number
  if(length(unique(.GlobalEnv[["rep_structure"]]$fraction)) > 1){
    evidence$rep_desc <- paste0('b',evidence$biorep,'t',evidence$techrep)
  }
  
  ## Generate Venn data for the identified proteins and output to a file
  levellog("read.protein.groups.data_v3: Generating ID Venn data ...")
  tmp.table<-data.table(evidence[, c('Protein.IDs', 'biorep', 'techrep', 'fraction')])
  setkey(tmp.table,  Protein.IDs, biorep, techrep, fraction)
  setwd(limma_output)
  write.table(tmp.table[, .(n=.N), by=.(Protein.IDs,rep=biorep)][,.(Protein.IDs,rep)],file=paste0(outputFigsPrefix,"_id_venn3-data_",time.point,".txt"),sep="\t",row.names=F)
  setwd("..")    
  
  # Bring Labeled or Label-free data to the following common format
  # (table headers):
  # rep_desc Protein.IDs UniqueSequences.Intensity.condition_1 ... UniqueSequences.Intensity.condition_N Intensity.condition_1 ... Intensity.condition_N
  
  levellog("read.protein.groups.data_v3: Standarizing data format ...")
  if(!is.proteome.discoverer.data){
    colnames(evidence.data)[grepl('Peptide.ID',colnames(evidence.data))]<-'Unique.Sequence.ID'
    if (!is.isobaric){
      # colnames(evidence.data)[grepl('Intensity\\..+',colnames(evidence.data))]<-conditions.labels
      colnames(evidence.data) <- sub('Intensity\\.(.+)', "\\1", colnames(evidence.data))
    }
    # else{
    # evidence[,conditions.labels]<-NA
    # for (my_cond in conditions.labels){
    #   mi<-which(grepl(my_cond, evidence$Labeling.State))
    #   evidence[mi, my_cond] <- evidence[mi, "Intensity"]
    # }
    # }
    
  }
  if (!'Unique.Sequence.ID' %in% colnames(evidence.data)){
    if ('Annotated.Sequence' %in% colnames(evidence.data)){
      colnames(evidence.data)[colnames(evidence.data) == 'Annotated.Sequence'] <- 'Unique.Sequence.ID'
      evidence$Unique.Sequence.ID <- sub(".*?\\.(.*?)\\..*", "\\1", evidence$Unique.Sequence.ID)
    }
  }
  if(is.label.free){
    if(is.proteome.discoverer.data){
      # Precursor Area is unfortunately buggy (sometimes 0/NA), so we are left with Intensity to work with
      #intensityCol <- 'Precursor.Area'
      intensityCol <- 'Intensity'
    }else{
      intensityCol <- 'Intensity'
    }
    evidence.dt<-data.table(evidence[, c('Protein.IDs', 'Unique.Sequence.ID', intensityCol,'label_', 'rep_desc')])
    setkey(evidence.dt, rep_desc, Protein.IDs, Unique.Sequence.ID, label_)
    # Get maximum PSM intensity per peptide/protein/[(rep_desc/label) = raw_file]
    suppressWarnings(evidence.dt<-evidence.dt[, .(maxI=max(get(intensityCol), na.rm = T)), by=.(rep_desc, Protein.IDs, Unique.Sequence.ID, label_)][maxI != -Inf])
  }else{
    if(is.proteome.discoverer.data){
      if (!'Quan.Usage' %in% colnames(evidence.data) & 'Peptide.Quan.Usage' %in% colnames(evidence.data)){
        colnames(evidence.data)[colnames(evidence.data) == 'Peptide.Quan.Usage'] <- 'Quan.Usage'
      }
      evidence.dt<-data.table(evidence[, c('Quan.Usage','Protein.IDs', 'Unique.Sequence.ID', conditions.labels,'rep_desc', 'label_')])
    }else{
      evidence.dt<-data.table(evidence[, c('Protein.IDs', 'Unique.Sequence.ID', conditions.labels,'rep_desc', 'label_')])
    }
    setkey(evidence.dt, rep_desc, Protein.IDs, Unique.Sequence.ID)    
  }
  ## Calculate identified peptide counts per protein for each condition/label and replicate in the following three steps
  # 1. For each condition (per sequnce, protein and replicate), set a corresponding column to TRUE if there are > 0 evidence.dt (PSMs) records, FALSE otherwise
  evidence.dt.seqCounts<-dcast.data.table(evidence.dt[, .(n=.N > 0), by=.(rep_desc, Protein.IDs, Unique.Sequence.ID, label_)], rep_desc + Protein.IDs + Unique.Sequence.ID ~ label_, fill=FALSE)
  # 2. Add a column flagging the common, between conditions/labels, sequences.
  # In case of more than two conditions/labels, the flag designates that there are at least two conditions/labels where the peptide is common
  evidence.dt.seqCounts[, 'common' := rowSums(.SD) > 1,.SDcols=conditions.labels]    
  # 3. Collapse the records for each protein (per replicate) and count the TRUEs.
  # evidence.dt[, .(n.Unique.Sequence.IDs=.N), by=.(rep_desc, Protein.IDs)]
  evidence.dt.seqCounts<-evidence.dt.seqCounts[,c(n.Unique.Sequence.IDs=.N,lapply(.SD, function(x){return(length(which(x)))})), by=.(rep_desc,Protein.IDs),.SDcols=c(conditions.labels, 'common')]
  # 4. Calculate the percentage columns
  evidence.dt.seqCounts[, paste0(conditions.labels,'p') := lapply(.SD, function(x){return((x/sum(.SD))*100)}), by=.(rep_desc,Protein.IDs),.SDcols=c(conditions.labels)]
  ## Rename the peptide counts columns
  setnames(evidence.dt.seqCounts,colnames(evidence.dt.seqCounts)[which(colnames(evidence.dt.seqCounts) %in% conditions.labels)],paste('UniqueSequences',conditions.labels,sep='.'))    
  ## Calculate the protein intensity = (sum of unique peptide intensities) for each condition/label and replicate in the following two steps
  if(is.label.free){
    # 1. Cast the data so that we have columns for each label and intensity separately
    evidence.dt<-dcast.data.table(evidence.dt, rep_desc + Protein.IDs + Unique.Sequence.ID ~ label_, fill=0)    
  }else{
    if(is.proteome.discoverer.data){
      # 1. Take the (Quan.Usage == 'Used') records and for each peptide keep only the PSM record with the highest intensity
      evidence.dt<-evidence.dt[Quan.Usage == 'Used' | Quan.Usage == 'Use', lapply(.SD, max), by=.(rep_desc, Protein.IDs, Unique.Sequence.ID), .SDcols=conditions.labels]    
    }else{
      # 2. Take the records with Intensity != NA across labels/conditions and for each peptide keep only the PSM record with the highest intensity
      evidence.dt[, sumI := rowSums(.SD, na.rm = T), .SDcols=conditions.labels]
      evidence.dt<-evidence.dt[sumI > 0, lapply(.SD, max), by=.(rep_desc, Protein.IDs, Unique.Sequence.ID), .SDcols=conditions.labels]    
      evidence.dt[, sumI := NULL]
    }
  }
  
  # Get a vector of unique peptides intensities
  tmp.I<-sort(unique(evidence.dt[,get(conditions.labels)]))
  # If the minimum intensity is zero
  if(tmp.I[1] == 0){
    # Replace 0's with minimum intensity (PD can do this automatically for us)
    minI<-tmp.I[2]
    evidence.dt[, (conditions.labels) := lapply(.SD, function(x){ t<-which(x == 0); if(length(t) > 0){x[t] <- minI}; return(x) }), .SDcols=conditions.labels]
  }else{
    minI<-tmp.I[1]
  }
  
  ## If enabled, do filter out peptides where all 'channels' except filterL_lbl channel have noise-level intensity (peptide-level filtering)
  if(filterL && filterL_lvl){
    evidence.dt[, minIcount := rowSums(.SD == minI), .SDcols=conditions.labels[! conditions.labels %in% filterL_lbl]]
    n1<-nrow(evidence.dt)
    evidence.dt<-evidence.dt[minIcount < (nConditions - 1)]
    n2<-nrow(evidence.dt)
    if(n2 < n1){
      levellog(paste0("read.protein.groups.data_v3: Filtered out ", (n1-n2)," peptides having noise-level intensity in all channels except the '", filterL_lbl,"' channel ..."));
    }
    evidence.dt[, minIcount := NULL]
  }
  
  # 2. Calculate the protein intensity (= sum of unique peptide intensities) for each condition/label and replicate
  # Also count the number of quantifiable peptides (those which do not have intensity NA)
  if(is.label.free){
    # Top three in abundance
    #evidence.dt<-evidence.dt[, lapply(.SD, function(x){x<-x[!is.na(x)]; x<-sort(x, decreasing<-T); if(length(x)<3){return(sum(x))}else{return(sum(x[1:3]))}}), by=.(rep_desc, Protein.IDs), .SDcols=conditions.labels]
    evidence.dt<-evidence.dt[, c(n=.N, nas=length(which(is.na(.SD))) ,lapply(.SD, function(x){x<-x[!is.na(x)]; x<-sort(x, decreasing<-T); if(length(x)<3){return(sum(x))}else{return(sum(x[1:3]))}})), by=.(rep_desc, Protein.IDs), .SDcols=conditions.labels]
  }else{
    # All peptides
    evidence.dt<-evidence.dt[, c(n=.N, nas=length(which(is.na(.SD))) ,lapply(.SD, sum, na.rm = T)), by=.(rep_desc, Protein.IDs), .SDcols=conditions.labels] 
  }
  ## Rename the intensity columns
  setnames(evidence.dt,colnames(evidence.dt)[which(colnames(evidence.dt) %in% conditions.labels)],paste('Intensity',conditions.labels,sep='.'))
  ## Merge with the evidence.dt.seqCounts table
  evidence.dt<-merge(evidence.dt, evidence.dt.seqCounts)
  
  # Add the experimental structure information to evidence.dt based on rep_desc (raw file at this point has no information and is dropped)
  tmp.rep_struct<-.GlobalEnv[["rep_structure"]][! duplicated(.GlobalEnv[["rep_structure"]][,c('biorep','techrep')]), !grepl('raw_file', colnames(.GlobalEnv[["rep_structure"]])) & !grepl('fraction', colnames(.GlobalEnv[["rep_structure"]]) )]
  tmp.rep_struct$rep_desc<-paste0('b',tmp.rep_struct$biorep,'t',tmp.rep_struct$techrep)
  evidence.dt<-merge(evidence.dt ,data.table(tmp.rep_struct), by='rep_desc')
  
  ## If enabled, do filter out proteins based on percentage labeling for the desired label (protein-level filtering)
  if(filterL && !filterL_lvl){
    n1<-length(unique(evidence.dt[get(paste0(filterL_lbl,"p")) == 100.0]$Protein.IDs))
    evidence.dt<-evidence.dt[get(paste0(filterL_lbl,"p")) < 100.0]
    levellog(paste0("read.protein.groups.data_v3: Filtered out ", n1," proteins which where identified solely by '", filterL_lbl, "' peptides ..."));
  }
  
  ## Get protein IDs that were quantified with a total of at least 'nRequiredLeastBioreps' different peptides accross at least 'nRequiredLeastBioreps' biological replicates.
  # E.g. 1: with 3 biological replicates, a protein that was quantified by a single peptide in 'nRequiredLeastBioreps' out of the 3 replicates will be discarded if 'nRequiredLeastBioreps' > 1 (retained otherwise).
  # E.g. 2: with 3 biological replicates, a protein that was quantified by a single peptide in 1 out of the 3 replicates will be discarded if 'nRequiredLeastBioreps' > 1 (retained otherwise).
  # E.g. 3: with 3 biological replicates, a protein that was quantified by two peptides in at one of replicates will be discarded if 'nRequiredLeastBioreps' > 1 (retained otherwise).
  # E.g. 4: with 3 biological replicates, a protein that was quantified by two peptides (in total) in 2 out of the 3 replicates will be discarded if 'nRequiredLeastBioreps' > 2 (retained otherwise).
  
  Protein.IDs.quant <- evidence.dt[, .(c1 = sum(N.x-nas)) , by=.(Protein.IDs, biorep)][, .(nQuantPeps = sum(c1), geqXnRequiredLeastBioreps = .N >= .GlobalEnv[["nRequiredLeastBioreps"]]), by=.(Protein.IDs)][nQuantPeps >= .GlobalEnv[["nRequiredLeastPeps"]] & geqXnRequiredLeastBioreps == T]$Protein.IDs
  levellog(paste0("read.protein.groups.data_v3: Filtered out ", (length(unique(evidence.dt$Protein.IDs)) - length(Protein.IDs.quant))," proteins which were not identified in at least ",nRequiredLeastBioreps," biological replicate(s) with at least a total of ",nRequiredLeastPeps," peptide(s)"));
  evidence.dt[,nQuantPeps := N.x-nas]
  evidence.dt<-evidence.dt[Protein.IDs %in% Protein.IDs.quant]
  
  ## Experimental filter based on outlier removal (grubbs method) based on the first condition specified.
  ## NOTE: It is applied when there are no technical replicates in Label-free data, where variability is expected to be very high.
  # If a protein intensity in condition i and biological replicate j is found to be an outlier based on the distribution
  # of intensities from all biological replicates, then
  # the biological replicate j is removed for that particular protein for all conditions.
  #if(is.label.free && .GlobalEnv[["n_techreps"]] < 2){
  #  evidence.dt.bad <- suppressWarnings(evidence.dt[, lapply(.SD, function(x){p.val = grubbs.test(x)$p.value; if(!is.na(p.val) && p.val < 0.05){outlier.true <- T}else{outlier.true <- F}; if(outlier.true){return(.I[outlier(x, logical=T)][1] )}else{return(as.integer(0))} }),by=.(Protein.IDs),.SDcols=paste0('Intensity.',conditions.labels[1])][,get(paste0('Intensity.',conditions.labels[1]))])
  #  evidence.dt.bad <- evidence.dt.bad[evidence.dt.bad > 0]
  #  evidence.dt<-evidence.dt[! evidence.dt.bad]
  #  levellog(paste0("read.protein.groups.data_v3: Filtered out ", length(evidence.dt.bad)," protein intensities based on outlier detection on condition '",conditions.labels[1],"'."));
  #  Protein.IDs.quant <- evidence.dt[, .(c1 = sum(n-nas)) , by=.(Protein.IDs, biorep)][, .(nQuantPeps = sum(c1), geqXnRequiredLeastBioreps = .N >= .GlobalEnv[["nRequiredLeastBioreps"]]), by=.(Protein.IDs)][nQuantPeps >= .GlobalEnv[["nRequiredLeastBioreps"]] & geqXnRequiredLeastBioreps == T]$Protein.IDs
  #  levellog(paste0("read.protein.groups.data_v3: Filtered out another ", (length(unique(evidence.dt$Protein.IDs)) - length(Protein.IDs.quant))," proteins which were not identified in at least ",nRequiredLeastBioreps," biological replicate(s) with at least a total of ",nRequiredLeastBioreps," peptide(s)"));
  #  evidence.dt[,nQuantPeps := n-nas]
  #  evidence.dt<-evidence.dt[Protein.IDs %in% Protein.IDs.quant]
  #}
  
  
  
  ## Generate Venn data for the identified proteins and output to a file
  levellog("read.protein.groups.data_v3: Generating quant Venn data ...")
  setwd(limma_output)  
  write.table(evidence.dt[, .(Protein.IDs, rep=biorep)],file=paste0(outputFigsPrefix,"_quant_venn3-data-",.GlobalEnv[["nRequiredLeastBioreps"]],"reps_",time.point,".txt"),sep="\t",row.names=F)
  setwd("..")
  
  ## Cast the table to the following format
  # Protein.IDs Intensity.[<rep_desc_X>.<label/condition_Y> ...] [<rep_desc_X>.Ratio.counts ...] [<rep_desc_X>.uniqueSequences ...] time.point [<label/condition_Y> ...] [<label/condition_Y>p ...]
  
  ## Step 1: For each 'rep_desc', add to a growing dataframe the evidence.dt data, renaming the columns accordingly
  # Also, calculate the missing columns required by the target format and drop the unnecessary columns
  setkey(evidence.dt, Protein.IDs)
  protein.groups.data<-data.frame(Protein.IDs = unique(evidence.dt)$Protein.IDs)
  setkey(evidence.dt, rep_desc)
  for(rep_desc_i in unique(evidence.dt)$rep_desc){
    rep_desc_i_protein.groups.data<-data.frame(evidence.dt[rep_desc == rep_desc_i,])
    allcols<-colnames(rep_desc_i_protein.groups.data)
    # Rename Intensity cols
    colsl<-grepl('^Intensity' ,allcols)
    colnames(rep_desc_i_protein.groups.data)[colsl]<-gsub("^Intensity(.+)$",paste("Intensity\\1",rep_desc_i,sep='.'), allcols[colsl])
    # Rename UniqueSequences cols
    colsl<-grepl('^UniqueSequences' ,allcols)
    colnames(rep_desc_i_protein.groups.data)[colsl]<-gsub("^UniqueSequences(.+)$",paste(rep_desc_i,"uniqueSequences\\1",sep='.'), allcols[colsl])
    # Add new column <rep_desc_X>.uniqueSequences
    rep_desc_i_protein.groups.data[, paste(rep_desc_i,'uniqueSequences',sep='.')]<-rowSums(rep_desc_i_protein.groups.data[, colnames(rep_desc_i_protein.groups.data)[colsl]])
    # Rename 'p' (percentage) cols
    colsl<-allcols %in% paste0(conditions.labels,'p')
    colnames(rep_desc_i_protein.groups.data)[colsl]<-gsub("^(.+)$",paste("\\1",rep_desc_i,sep='.'), allcols[colsl])
    # Rename the 'nQuantPeps' column to <rep_desc_i>.Ratio.counts
    colsl<-allcols %in% c('nQuantPeps')
    colnames(rep_desc_i_protein.groups.data)[colsl]<-paste(rep_desc_i,'Ratio.counts',sep='.')
    # merge with the growing data frame
    cc<-intersect(names(protein.groups.data), names(rep_desc_i_protein.groups.data))
    protein.groups.data<-merge(protein.groups.data, rep_desc_i_protein.groups.data[, ! colnames(rep_desc_i_protein.groups.data) %in% c('biorep', 'techrep', 'fraction', 'rep_desc', cc[! grepl('Protein.IDs', cc)] )], all.x = T)
  }
  # Step 2: Calculate the columns [<label/condition_Y> ...] containing the number of unique sequences found per condition in all replicates
  allcols<-colnames(protein.groups.data)
  for(cond_i in conditions.labels){
    colsl<-grepl(paste('uniqueSequences', cond_i,sep='\\.') ,allcols)
    protein.groups.data[, cond_i]<-rowSums(protein.groups.data[, allcols[colsl]])
  }
  # Step 3: Calculate the columns [<label/condition_Y>p ...] containing the percentage of unique sequences that were found in a specific condition in all replicates
  allcols<-colnames(protein.groups.data)
  for(cond_i in conditions.labels){
    colsl<-allcols %in% conditions.labels & ! allcols %in% cond_i
    protein.groups.data[, paste0(cond_i,'p')]<-(protein.groups.data[, cond_i]/rowSums(protein.groups.data[, c(cond_i, allcols[colsl])]))*100
  }
  # Step 4: Add time-point column
  protein.groups.data$time.point<-time.point
  # Step 5: Remove unnecessary columns (uniqueSequences per rep_desc and percentage unique peptides per rep_desc)
  allcols<-colnames(protein.groups.data)
  protein.groups.data<-protein.groups.data[,-which(grepl('uniqueSequences\\.', allcols) | grepl('p\\.b[0-9]+t[0-9]+$', allcols) | grepl('^common$', allcols))]
  ##
  levellog(paste0("read.protein.groups.data_v3: Quantifiable proteins: ", nrow(protein.groups.data)," (",time.point,")"))
  levellog("",change=-1)
  ## 
  return(protein.groups.data)  
}


perform_analysis <- function() {
  # Moved to initialize.R
  # ProteinQuantitation <- TRUE
  
  # rep_structure<-read.table(experimental_structure_file,col.names=c("raw_file","biorep","techrep","fraction"), sep="\t")
  # Moved to initialize.R
  # replicate.multiplexing.is.used <- FALSE
  
  # Added to load_data_R
  # experimental.structure <- read.csv("~/experimental-structure.csv")
  # 
  # experimental.structure <- experimental.structure[
  #   order(experimental.structure$conditions,
  #         experimental.structure$biorep,
  #         experimental.structure$techrep,
  #         experimental.structure$fractions),
  #   ]
  # 
  # rownames(experimental.structure) <- c(1:length(experimental.structure$raw.file))
  # 
  # biological.replicates.list <- experimental.structure$biorep
  # technical.replicates.list <- experimental.structure$techrep
  # experimental.fraction.list <- experimental.structure$fraction
  # experimental.conditions.list <- experimental.structure$condition
  # 
  # biological.replicates.number.status <- check.number.of.replicates(replicate.multiplexing.is.used, biological.replicates.list)
  # 
  # if (biological.replicates.number.status  == FALSE) {
  #   cat("Error User: Cannot accept dataset with just one biological replicate. Aborting ...\n")
  #   return (FALSE)
  # }
  # 
  # replicates.per.condition <- replicates.per.condition(biological.replicates.list,
  #                                                      technical.replicates.list,
  #                                                      experimental.conditions.list)
  # 
  # replicates.status.per.condition <- replicates.status.per.condition(replicates.per.condition)
  # 
  # 
  # fixed.replicates.per.condition <- fix.replicates.per.condition(replicates.per.condition,
  #                                                                replicates.status.per.condition)
  # 
  # replicates.per.condition <- fixed.replicates.per.condition
  # 
  # restored.replicates <- restore.replicates(replicates.per.condition)
  # 
  # biological.replicates.list <- restored.replicates$biological
  # 
  # technical.replicates.list <- restored.replicates$technical
  # 
  # experimental.structure$biorep <- biological.replicates.list
  # 
  # experimental.structure$techrep <- technical.replicates.list
  # 
  
  # TODO Maybe later for condition fix
  fixed.labels <- conditions.labels

  bad.labels.index <- which(grepl("^(?!cond).*$", fixed.labels, perl = TRUE))

  fixed.labels[bad.labels.index] <- paste0("cond", "N")

  # TODOEND
  
  # Made it my own way
  # original_rep_structure$rep_desc<-paste(paste(paste("b",original_rep_structure$biorep,sep=""),"t",original_rep_structure$techrep,sep=""))
  # 
  # Moved to build.R 
  # experimental.setup.id <-  1 * check.number.of.replicates(FALSE, biological.replicates.list) +
  #   2 * check.number.of.replicates(FALSE, technical.replicates.list) +
  #   3 * check.number.of.replicates(FALSE, experimental.fraction.list)
  # 
  # experimental.description <- make.experimental.description(experimental.setup.id,
  #                                                           biological.replicates.list,
  #                                                           technical.replicates.list,
  #                                                           experimental.fraction.list)
  # 
  # experimental.structure$description <- experimental.description
  # 
  # Removed
  # clean.workspace()
  
  # ### Huh?
  # .GlobalEnv[["rep_structure"]]<-rep_structure
  # .GlobalEnv[["LFQ_conds"]]<-LFQ_conds
  # .GlobalEnv[["original_rep_structure"]]<-original_rep_structure
  # .GlobalEnv[["n_bioreps"]]<-max(rep_structure$biorep)
  # .GlobalEnv[["n_techreps"]] <- min(
  #   ddply(rep_structure[,c("biorep","techrep")],
  #         c("biorep"),
  #         function(x){
  #           return(max(x$techrep))
  #         })$V1)
  
  # ###
  
  # Should be setted on initialize.R  
  # if (ProteinQuantitation == TRUE) {
  #   quantitation.status <- "Protein"
  # } else {
  #   quantitation.status <- "Peptide"
  # }
  # Moved to analyze.R/functions_analyze.R
  # if( file.exists(limma.output.folder) == TRUE) {
  #   unlink(limma.output.folder, recursive=T, force=T)
  # }
  # 
  # dir.create(limma.output.folder)
  # 
  # Moved to load_data.R
  # evidence.file <- "Dropbox/Review/case studies/evidence.txt"
  # protein.groups.file <- "Dropbox/Review/case studies/cleanProteinGroups.txt"
  
  # Moved to initialize.R
  # is.proteome.discoverer.data <- FALSE
  
  cat("Calling read.protein.groups.data ...\n")
  protein.groups <- read.protein.groups.data_v3(protein.groups.file,
                                                evidence.file,
                                                time.points
                                                keep.evidences.ids=T,
                                                is.proteome.discoverer.data)
  
  #Restore the original rep descriptions to add to the graph
  if (replicate.multiplexing.is.used == FALSE)
  {
    newcolumns <- names(protein_groups)
    oldcolumns = newcolumns
    for(my_column in newcolumns){
      for(my_repdesc in .GlobalEnv[["rep_structure"]]$rep_desc){
        if (grepl(my_repdesc, my_column)){
          temp_name <- .GlobalEnv[["original_rep_structure"]]$rep_desc[match(my_repdesc, .GlobalEnv[["rep_structure"]]$rep_desc)]
          newcolumns[match(my_column, newcolumns)] <- sub(my_repdesc, temp_name, my_column)
        }
      }
    }
    colnames(protein_groups) <- newcolumns
  }
  
  setwd(limma.output.folder)
  temp_protein.groups.data <- protein_groups
  write.table(temp_protein.groups.data[, -which(names(temp_protein.groups.data) %in% c("N.x","N.y"))],file=paste(outputFigsPrefix,"_proteinGroupsDF.txt",sep=""),row.names=F,sep="\t")
  setwd("..")
  if (!replicate.multiplexing.is.used)
  {
    colnames(protein_groups) <- oldcolumns
  }
  #Create the expdesign table:
  expdesign<-c()
  #Rename conditions labels from condN back to N
  for(i in 1:length(conditions.labels)){
    if(conditions.labels[i] == "condN")
    {
      conditions.labels[i] <- "N"
      conditions.labels <<- conditions.labels
      colnames(protein_groups) <- sub("condN", "N", colnames(protein_groups))
    }
  }
  if(filterL_lbl == "condN")
  {
    filterL_lbl <<- "N"
  }
  #Rename condN back to N in lfq_conds and protein_groups
  mi <- which(LFQ_conds$condition == "condN")
  if(length(mi)>0)
  {
    levels(LFQ_conds$condition) <- c(levels(LFQ_conds$condition), "N")
    LFQ_conds$condition[which(LFQ_conds$condition == "condN")] <- "N"
    LFQ_conds$condition <- factor(LFQ_conds$condition)
    LFQ_conds <<- LFQ_conds
    colnames(protein_groups) <- sub("condN", "N", colnames(protein_groups))
  }
  for(cond_i in conditions.labels){
    expdesign<-rbind(expdesign,cbind(paste(sub("Intensity\\.","",sort(colnames(protein_groups)[grep(paste("Intensity.",cond_i,".b",sep=""),colnames(protein_groups))]))),cond_i))  
  }
  
  colnames(expdesign)<-c("Sample","Category")
  if(!replicate.multiplexing.is.used){
    #the following lines also deal with restoring the original breps and treps numbers
    #temp vector has only the information of the replication (e.g. b1t1 or b1t1f1 if we have fractionation)
    temp_vector <- sub("(.*)\\.","", expdesign[,1])
    temp_vector <- original_rep_structure$rep_desc[match(temp_vector, sub("f.*", "", rep_structure$rep_desc))]
    #Make sure that expdesign (column Sample) contains data in te right format by merging expdesign and tmp_vector:
    tmp_counter <- 0
    for (expdesign_i in expdesign[,1]){
      expdesign[tmp_counter + 1,1] <- sub("(.*)\\..*",paste0("\\1.", temp_vector[tmp_counter + 1]), expdesign_i)
      tmp_counter <- tmp_counter + 1
    }
  }
  #Remove the fractionation information: (if any)
  expdesign[,1] <- sub("(.*\\..*)f.*", "\\1", expdesign[,1], perl = TRUE)
  write.table(expdesign,file="curr_exp_design.txt",row.names=F,quote=F,sep = "\t")
  exp_design_protein.groups.file<<-"curr_exp_design.txt"
  
  levellog("Performing the analysis ...")
  
  do_limma_analysis(prepare_working_protein.groups.data(protein_groups),time.point,exp_design_protein.groups.file,exportFormat="pdf",outputFigsPrefix=outputFigsPrefix)
  
  levellog("Data analysis finished.")
  levellog("",change=-1)
  return(T)
}

##########################  CLEAN  #########################################
logger <- NULL

initialize.logger <- function(){
  logger <<- create.logger(level = "DEBUG")
}

do.logging <- function(message, level){
  memory <- humanReadable(mem_used(), standard = "SI")
  message <- paste(message," ", memory, ".", sep = "")
  levellog(logger, level, message)
}

# Used the global.variable for cleaning purposes
# clean.workspace <- function() {
#   variables <- c("biological.replicates.number.status",
#                  "fixed.replicates.per.condition",
#                  "replicates.status.per.condition",
#                  "restored.replicates",
#                  "experimental.setup.id")
#   
#   functions <- c("check.number.of.replicates",
#                  "check.replicates",
#                  "fix.replicates",
#                  "fix.replicates.per.condition",
#                  "make.experimental.description",
#                  "restore.replicates")
#   r.objects <- c(variables, functions)
#   suppressWarnings(remove(list = r.objects, envir = globalenv()))
# }

