#####################################################################
# Septic Shock data													#
# Signaling Pathway Impact Analysis									#
#####################################################################
# SPIA of all studies
spia.file <- "spia.res_2016.05.02.rda"
if(file.exists(spia.file)) {
  load(file=spia.file)
} else {
  spia.res <- list()
  print.noquote("Running SPIA... Process will take Time...")
  for(i in 1:length(ssids)) {
    study.id <- ssids[i]
    print.noquote(study.id)
    eset <- ss.eset[[study.id]]
    rtt <- rowttests(eset, "Group")
    #rtt <- all.con.rtt[[study.id]]
    sel <- which(rtt$p.value<0.05)
    egs.sel <- rownames(rtt)[sel]
    lfc.sel <- rtt[egs.sel, "dm"]
    names(lfc.sel) <- egs.sel
    universeGeneIds <- rownames(rtt)

    # SPIA
    res<-spia(de=lfc.sel, all=universeGeneIds, organism="hsa", nB=2000,
	  plots=F, beta=NULL, combine="fisher", verbose=FALSE)
    spia.res[[i]] <- res
  }
  spia.date <- date()
  save(spia.res, spia.date, file=spia.file)
}

res <- spia.res[[1]]
resall <-  data.frame(ssids[1], res$ID,res$pG,res$Name)
colnames(resall)<-c("Study", "KEGG_id","pG_score","KEGG_name")
for(i in 2:length(ssids)) {
  res <- spia.res[[i]]
  resedited<-data.frame(ssids[i], res$ID,res$pG,res$Name)
  colnames(resedited)<-c("Study", "KEGG_id","pG_score","KEGG_name")
  resall <- rbind(resall, resedited, deparse.level=0)
}

# Calculate Fisher's product of p-values for all pathways
# Check if the pG_score can be considered a p-value
keggs <- as.character(unique(resall$KEGG_id))
spiaFp <- vector(mode="numeric", length=length(keggs))
names(spiaFp) <- keggs
for(id in keggs) {
  ps <- resall[resall$KEGG_id==id, "pG_score"]
  spiaFp[id] <- Fisher.test(ps)["p.value"]
}

# Get the top 10 perturbed pathways
top10ids <- names(sort(spiaFp, decreasing=F)[1:10])
#print("The pathways discovered by SPIA of septic shock studies\n")
#print(resall[match(top10ids,resall[,2]),c(2,4)])

#    KEGG_id                            KEGG_name
# 1    05134                        Legionellosis
# 2    05152                         Tuberculosis
# 7    05132                 Salmonella infection
# 5    04010               MAPK signaling pathway
# 3    04380           Osteoclast differentiation
# 20   05150      Staphylococcus aureus infection
# 10   05140                        Leishmaniasis
# 22   04620 Toll-like receptor signaling pathway
# 13   05322         Systemic lupus erythematosus
# 9    05145                        Toxoplasmosis

