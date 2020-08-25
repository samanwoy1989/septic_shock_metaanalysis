# A function for reading data
# Read the bead-level data for Set1
# Normalize and store the expression set object
# Variance based filtering to remove duplicated probes
# Author: Saroj Kant Mohapatra
# March 2015

readData <- function(targets.file = "Metadata/TARGETS1.txt",
                     eset.file = "Data/sepsis_1_20150717.rda",
                     data.file = "TableControl_sample probe profile.txt",
                     ctrl.data.file = "TableControl_control probe profile.txt",
                     data.path="Set1",
                     study.name = "Sepsis_SCB_Batch_1") {

  if(data.path=="Set3") {
    source("Rcode/read_set3.R")
  } else {
    targets <- readTargets(file=targets.file, row.names="ID")
    if(file.exists(eset.file)) {
      load(eset.file)
    } else {
      x1 <- read.ilmn(files=data.file,
                      ctrlfiles=ctrl.data.file, 
                      path=data.path, ctrlpath=data.path)
      if(data.path=="nibmg/Set2") {
        x1 <- x1[, colnames(x1)!="9553976020_K"]       #x1 <- x1[,-11] # remove sample K
        #targets <- targets[-11,]
      }
      if(data.path=="nibmg/Set7") {
        x1 <- x1[, colnames(x1)!="C15"]       #x1 <- x1[,-11] # remove sample C15
      }
      # Normalize
      y1 <- neqc(x1)
      # Add Entrez gene ids to each probe to set1
      allgsyms <- y1$genes[,"SYMBOL"]
      allgsyms[allgsyms==""] <- "NEGATIVE"
      allegids <- mget(allgsyms, revmap(org.Hs.egSYMBOL), ifnotfound=NA)
      allegids <- sapply(allegids, function(x) x[[1]][1])
      uniegids <- sort(unique(allegids))
      # For each egid, get a single row id with
      # duplicated probes are resolved by selecting the high-variance probe id
      ids.vf <- sapply(uniegids, function(eg) {
        ids <- which(y1$genes[,"SYMBOL"]%in%names(which(allegids==eg)))
        ifelse(length(ids)>1,
               ids[which.max(apply(as.matrix(y1$E[ids,]), 1, var, na.rm=T))], ids)
      })
      exps <- y1$E[ids.vf,]
      rownames(exps) <- uniegids
      
      common.ids <- intersect(colnames(exps), rownames(targets))
      targets <- targets[common.ids,]
      exps <- exps[,common.ids]
      rm(common.ids)
      pData <- data.frame(targets)
      rownames(pData) <- targets[,"Sample_ID"]
      if(all(pData$ID==colnames(exps))) {
        colnames(exps) <- rownames(pData)
      } else {
        cat("Error!! Column names in exps matrix not matching pData.\n")
      }
      
      phenoData <- new("AnnotatedDataFrame", data=pData)
      experimentData <- new("MIAME", name="Saroj K Mohapatra",
                            lab="Sepsis Genomics Lab",
                            contact="skm1@nibmg.ac.in",
                            title=study.name,
                            other=list(notes="Created from Illumina probe profile files "))
      annotation <- "org.Hs.eg.db"
      eset <- ExpressionSet(assayData=exps, phenoData=phenoData,
                            experimentData=experimentData,
                            annotation=annotation)
      save(eset, file=eset.file)
    } 
  }
  eset
}
