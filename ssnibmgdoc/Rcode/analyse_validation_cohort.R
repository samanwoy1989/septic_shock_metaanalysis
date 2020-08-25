eset <- gset
pdata <- pData(eset)

dx.str <- sapply(strsplit(as.character(pdata[, "characteristics_ch1.2"]),
                          split=": "), "[", 2)                    
dx <- dx.str
dx[dx.str=="NA"] <- "Control"
dx[dx.str!="NA"] <- "Septic_Shock"    

ptids <- sapply(strsplit(as.character(pdata$title), split="_"),"[", 2)

day.str <- sapply(strsplit(as.character(pdata[, "characteristics_ch1.3"]),
                           split=": "), "[", 2)
day <- sapply(strsplit(day.str, split=" hr"), "[", 1)
day[day=="0"]  <- 0
day[day=="24"] <- 1
day[day=="48"] <- 2            

pData(eset) <- data.frame(pdata, ptids, dx.str, dx, day)  	
vset <- eset

egs25<-egs.toshow



is.pathgene <- which(featureData(vset)$"ENTREZ_GENE_ID"%in%egs25)
egs25<-unique(featureData(vset)$"ENTREZ_GENE_ID"[is.pathgene])



gexp <-exprs(vset)
vars <- apply(gexp, 1, var)
probes25 <-c()
for(i in 1:length(egs25)) {
  eg <- egs25[i]
  probes25[i] <- names(which.max(vars[(featureData(vset)$ENTREZ_GENE_ID==eg)]))
}
names(probes25) <- egs25 

dat<-gexp[probes25,]
