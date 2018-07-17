
addLabel<-function(lblname, lbl.Modifications, LabelFree = T){
  
  #levellog("", change=1)
  #If label name is a number some routines won't work, it has to be converted to some acceptable variable name
  lblname<-make.names(lblname)
  labeltxt <- "label";
  if(!LabelFree){
    lbl.Modifications<-gsub("\\(","\\\\(",lbl.Modifications)
    lbl.Modifications<-gsub("\\)","\\\\)",lbl.Modifications)
    unmod_idx<-which(lbl.Modifications == "")
    if(length(unmod_idx) > 0){
      rest_idx<-which(lbl.Modifications[-unmod_idx] != "")
      if(length(rest_idx) > 0){
        lbl.Modifications<-c(lbl.Modifications[unmod_idx],paste(lbl.Modifications[-unmod_idx],"\\)",sep=""))
      }
    }else{
      lbl.Modifications<-paste(lbl.Modifications,"\\)",sep="")
    }
  }else{
    labeltxt <- "condition";
  }
  cat("condition labels", conditions.labels, "<<<\n")
  lblname_i<-which(grepl(paste("^",lblname,"$",sep=""),conditions.labels))
  print(lblname_i)
  if(length(lblname_i) != 0){
    cat(paste("addLabel: Error adding ",
              labeltxt," '",
              lblname,
              "': An existing ",
              labeltxt,
              " with name '",
              lblname,
              "' (specification: ",
              paste(unlist(conditions.labels.Modifications[lblname_i]),
                    collapse=", "),
              ") already exists. Please try a different name.",
              sep=""))
    return(FALSE)
  }
  
  conditions.labels<<-c(conditions.labels, lblname)
  j<-length(conditions.labels.Modifications)+1
  conditions.labels.Modifications[[j]]<<-lbl.Modifications
  nConditions<<-length(conditions.labels)
  #cat("", change=-1)
}

clearLabels<-function(){
  #levellog("", change=1)
  conditions.labels<<-c()
  conditions.labels.Modifications<<-list()
  nConditions<<-length(conditions.labels)
  #levellog("", change=-1)
}

conditions.to.raw.files <- list()
