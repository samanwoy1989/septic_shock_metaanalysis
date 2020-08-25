require("resample")
simulp.file <- "res_per.2016.04.19.rda"
if(file.exists(file=simulp.file)) {
  load(file=simulp.file)
} else {

  egs.hsa04380 <- genes.by.pathway$hsa04380 # hsa04380 genes
  egs.sel <- intersect(egs.hsa04380, featureNames(esetn)) # should be in NIBMG data
  egs.sel <- intersect(egs.sel, allegs) # present also in metaanalysis data

  grpids <- as.character(esetn$Group)
  dat <- data.frame(grpids, t(exprs(esetn[egs.sel,])))
  colnames(dat)[1] <- "Group"
  res_per <- permutationTest(data=dat, alternative="greater", 
	statistic=function(mydat) {
      is.con <- mydat$Group=="Control"
      is.ss <- mydat$Group=="Septic Shock"
      x <- t(mydat[, colnames(dat)!="Group"])
      pvalsn <- apply(x, 1, function(x) {
		t.test(x[is.con], x[is.ss], alternative="less")$p.value
	  })
	  
      padjn <- pvalsn
      xcon <- rowMeans(x[, is.con])
      xss <- rowMeans(x[, is.ss])
      which.lowp <- which(padjn<0.05)
      length(which.lowp)/length(padjn)
    },
	resampleColumns = "Group", R=99999)
  save(res_per, file=simulp.file)
}
pstr <- paste("p = ", formatC(res_per$stat$PValue, digits=3), sep="")

