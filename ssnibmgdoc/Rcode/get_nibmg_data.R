# get_nibmg_data.R
# Jan 14 2016	Included batch effect correction
# Mar 16, 2016  

nibmg.file <- "nibmg20160317.rda"
if(file.exists(file=nibmg.file)) {
  cat("Loading nibmg gene expression data ...")
  load(file=nibmg.file)
  cat(" done!\n")
} else {
  require(limma)

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
    if(i==7) {
      datestr <- "20160129"
    } else {
      datestr <- "20150720"
    }
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

  pdata <- pData(eset)
  diagnos <- as.character(pdata$Diagnosis)
  diagnos[is.na(diagnos)] <- "Sepsis"
  diagnos[!diagnos%in%c("Sepsis","Severe Sepsis","Septic Shock")] <- "Control"
  ptids <- as.character(pdata$PTID)
  ptids[ptids=="ctrl_sm"] <- "C1"
  not.starting.C <- substr(ptids, 1, 1)!="C"
  ptids[not.starting.C] <- paste("S", ptids[not.starting.C], sep="")
  pData(eset) <- data.frame("Diagnosis"=factor(diagnos),
  			"PTID" = factor(ptids),
  			pdata[, c("Group","batchid")])
  

  # Remove batch effect
  esetnibmg.b <- eset
  cat("Performing batch effect correction ...")
  pheno <- pData(eset)
  edata <- exprs(eset)
  batch <- as.integer(eset$batchid)
  mod <- model.matrix(~as.factor(Group), data=pheno)
  combat_edata = ComBat(dat=edata, batch=batch, mod=mod, 
	par.prior=TRUE, prior.plots=FALSE)
  exprs(eset) <- combat_edata
  cat(" done!\n")
  esetnibmg.a <- eset

  # Keep only control and septic shock day 1 (group == 0)
  #is.of.nibmg <- !(eset$PTID%in%c("CMC2", "CMC3", "ks1", "ks2"))
  which.con <- which(eset$Group=="Ctrl")
  ssptids <- as.character(unique(eset$PTID[which(eset$Diagnosis=="Septic Shock")]))
  which.ss <- match(ssptids, eset$PTID)
  sel <- c(which.con, which.ss)
  esetn <- esetnibmg.a[,sel]
  esetn.b <- esetnibmg.b[,sel]
  
  sampleNames(esetn) <- as.character(esetn$PTID)
  pData(esetn) <- data.frame("Group"=esetn$Diagnosis, "PTID"=esetn$PTID,
  	"batchid"=esetn$batchid)
  esetn$Group <- factor(esetn$Group)
  esetn$PTID <- factor(esetn$PTID)
  esetn$batchid <- factor(esetn$batchid)

  sampleNames(esetn.b) <- as.character(esetn.b$PTID)
  pData(esetn.b) <- data.frame("Group"=esetn.b$Diagnosis, "PTID"=esetn.b$PTID,
  	"batchid"=esetn$batchid)
  esetn.b$Group <- factor(esetn.b$Group)
  esetn.b$PTID <- factor(esetn.b$PTID)
  esetn.b$batchid <- factor(esetn.b$batchid)
  
  # Save data
  cat("Saving data to file ...")
  nibmg.date <- date()
  nibmg.info <- sessionInfo()
  save(esetn, esetn.b, nibmg.date, nibmg.info, file=nibmg.file)
  cat(" done!\n")
}
