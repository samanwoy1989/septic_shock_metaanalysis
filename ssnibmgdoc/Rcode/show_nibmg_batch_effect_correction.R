# show_nibmg_batch_effect_correction.R
# Show that batch effect has been removed from nibmg data

  # Function definition for reading NIBMG data
  source("nibmg/Rcode/readData.R")

  #Reading the data files
  params.lst <- read.table(file="nibmg/Metadata/params_list.txt", header=T, 
                     stringsAsFactors=FALSE, sep="\t")
  eset.list <- list()
  for(i in c(1:nrow(params.lst))) { # Avoid number 3
    print(i)

    targets.file=paste("nibmg/Metadata/TARGETS", i, ".txt", sep="")
    #datestr <- strsplit(as.character(Sys.time()), split=" ")[[1]][1]
    #datestr <- gsub("-","",datestr)
    datestr <- "20150720"
    eset.file <- paste("./nibmg/Data/sepsis_", i, "_", datestr, ".rda", sep="")
    data.path <- paste("nibmg/Set", i, sep="")

    if(i==3) {
      source("nibmg/Rcode/read_set3.R")
      eset.list[[i]] <- eset
      rm(eset)
      gc()
    } else {
    
      eset.list[[i]] <- readData(targets.file=targets.file,
           eset.file=eset.file,
           data.file=params.lst$data.files[i],
           ctrl.data.file=params.lst$ctrl.data.files[i],
           data.path=data.path,
           study.name=params.lst$study.name[i])
    }
  }

  # Combine the expression data for common genes into one eset
  # also add a batch id to each array
  egs <- featureNames(eset.list[[1]])
  for(i in 2:length(eset.list)) {
    egs <- intersect(egs, featureNames(eset.list[[i]]))
  }

  for(i in 1:length(eset.list)) {
    eset <- eset.list[[i]]
    eset$batchid <- rep(i, ncol(eset))
    eset.list[[i]] <- eset[egs,]
    rm(eset)
    gc()
  }

  eset <- eset.list[[1]]
  gexp <- exprs(eset)
  ID <- eset$ID
  Sample_ID <- eset$Sample_ID
  Group <- eset$Group
  PTID <- eset$PTID
  Age <- eset$Age
  Sex <- eset$Sex
  Diagnosis <- eset$Diagnosis
  batchid <- eset$batchid
  pData <- data.frame(cbind(ID,Sample_ID,Group,PTID, Age, Sex, 
  	Diagnosis, batchid))
  rownames(pData) <- pData$Sample_ID

  for(i in 2:length(eset.list)) {
    eset <- eset.list[[i]]
    gexp <- cbind(gexp, exprs(eset))
    ID <- eset$ID
    Sample_ID <- eset$Sample_ID
    Group <- eset$Group
    PTID <- eset$PTID
    Age <- eset$Age
    Sex <- eset$Sex
    Diagnosis <- eset$Diagnosis
    batchid <- eset$batchid
    pData.new <- 
  	  data.frame(cbind(ID,Sample_ID,Group,PTID, Age, Sex, Diagnosis, batchid))
    rownames(pData.new) <- pData.new$Sample_ID
    pData <- rbind(pData, pData.new)
  }

  if(!all(colnames(gexp)==rownames(pData))) {
    cat("Error!! Column names in expression matrix not matching pData.\n")
  }
  phenoData <- new("AnnotatedDataFrame", data=pData)
  experimentData <- new("MIAME", name="Saroj K Mohapatra",
                      lab="Sepsis Genomics Lab",
                      contact="skm1@nibmg.ac.in",
                      title="Sepsis_SCB_combined_March_2015",
                      other=list(notes="Data are from six batches. Batch 1,2,4,5,6: Created from neqc-normalized Illumina data; Batch 3: Affymetrix RMA-normalized data. Duplicated probes have been removed. Feature name corresponds to human Entrez gene id."))
  annotation <- "org.Hs.eg.db"
  eset <- ExpressionSet(assayData=gexp, phenoData=phenoData,
                      experimentData=experimentData,
                      annotation=annotation)

  eset.before <- eset
  

  # Remove batch effect
  cat("Performing batch effect correction ...")
  pheno <- pData(eset)
  edata <- exprs(eset)
  batch <- as.integer(eset$batchid)
  mod <- model.matrix(~as.factor(Group), data=pheno)
  combat_edata = ComBat(dat=edata, batch=batch, mod=mod, 
	par.prior=TRUE, prior.plots=FALSE)
  exprs(eset) <- combat_edata
  cat(" done!\n")

  # remove a chip of questionable quality
  #eset <- eset[, !(sampleNames(eset)%in%c("SKM1", "CMC2","CMC3"))]
  #eset <- eset[, !eset$PTID%in%c("ks1","ks2")] # remove non-SCB samples
  eset$Group <- factor(eset$Group)
  
  eset.after <- eset
  
  # Boxplots before and after batch effect correction
  png("Result/nibmg_batch_effect_correction.png", width=10)
  par(mfrow=c(1,2))
  boxplot(exprs(eset.before), las=2, names=1:ncol(eset), 
  	main="Before correction")  	 
  boxplot(exprs(eset.after), las=2, names=1:ncol(eset), 
  	main="After correction")
  dev.off()
  #system("evince Result/nibmg_batch_effect_correction.pdf")

