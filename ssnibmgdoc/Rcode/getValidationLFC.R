# Get the log-fold changes for the validation data sets
# 1 - GENE
# 2 - ENTREZ_GENE_ID
# 3 - GENE_SYMBOL

egs <- egs25
gsyms <- as.character(unlist(mget(egs, org.Hs.egSYMBOL)))
lfcV <- data.frame(matrix(NA, nrow=length(egs), ncol=3))
colnames(lfcV) <- names(gsets)
rownames(lfcV) <- egs

for(i in 1:length(egs)) {
  eg <- as.character(egs[i])
  sym <- as.character(gsyms[i])
  print(eg)
  # study 1
  study.id <- names(gsets)[1]
  gset <- gsets[[study.id]]
  selcols <- gset$dx%in%c("Control", "Septic_Shock")
  sel <- which(eg==featureData(gset)$GENE)
  probe <- featureNames(gset)[sel]
  if(length(probe)!=0) {
    if(length(probe)>1) {
      probe <- names(which.max(apply(exprs(gset[probe,]), 1, var)))
    }
  }
  gxp <- exprs(gset[probe, selcols])
  colnames(gxp) <- gset$dx[selcols]
  xCon <- mean(gxp[colnames(gxp)=="Control"], na.rm=T)
  xSS  <- mean(gxp[colnames(gxp)=="Septic_Shock"], na.rm=T)
  lfcV[eg, study.id] <- xSS-xCon
  
  # study 2
  study.id <- names(gsets)[2]
  gset <- gsets[[study.id]]
  selcols <- gset$dx%in%c("Control", "Septic_Shock")
  sel <- which(eg==featureData(gset)$ENTREZ_GENE_ID)
  probe <- featureNames(gset)[sel]
  if(length(probe)!=0) {
    if(length(probe)>1) {
      probe <- names(which.max(apply(exprs(gset[probe,]), 1, var)))
    }
  }
  gxp <- exprs(gset[probe, selcols])
  colnames(gxp) <- gset$dx[selcols]
  xCon <- mean(gxp[colnames(gxp)=="Control"], na.rm=T)
  xSS  <- mean(gxp[colnames(gxp)=="Septic_Shock"], na.rm=T)
  lfcV[eg, study.id] <- xSS-xCon

  # study 3
  study.id <- names(gsets)[3]
  gset <- gsets[[study.id]]
  selcols <- gset$dx%in%c("Control", "Septic_Shock")
  sel <- which(sym==featureData(gset)$GENE_SYMBOL)
  probe <- featureNames(gset)[sel]
  if(length(probe)!=0) {
    if(length(probe)>1) {
      probe <- names(which.max(apply(exprs(gset[probe,]), 1, var)))
    }
  }
  gxp <- exprs(gset[probe, selcols])
  colnames(gxp) <- gset$dx[selcols]
  xCon <- mean(gxp[colnames(gxp)=="Control"], na.rm=T)
  xSS  <- mean(gxp[colnames(gxp)=="Septic_Shock"], na.rm=T)
  lfcV[eg, study.id] <- xSS-xCon
}


