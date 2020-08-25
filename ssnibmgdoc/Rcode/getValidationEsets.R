# Keep all the three validation esets in one list
# 2016.06.27 -  Added hours post-shock onset for GSE57065

valfn <- "Data/validation.rda"
if(file.exists(file=valfn)) {
  cat("Loading validation data ...")
  load(valfn)
  cat(" done!\n")
} else {
  # There are three gse ids
  gseids.val <- c("GSE48080","GSE57065","GSE67652")
  vsets <- list()
  for(i in 1:length(gseids.val)) {
    gseid <- gseids.val[i]
    print(gseid)
    fn <- paste("Data/", gseid, ".rda", sep="")
    if(file.exists(fn)) {
      cat("Loading data from file ... ")
      load(fn)
      cat("done!\n")
    } else {
      cat("Reading the smatrix file and converting to a .rda ...")
      gset <- getGEO(filename=
  	    paste("Data/", gseid, "_series_matrix.txt.gz", sep=""),
  	    GSEMatrix=TRUE)
      save(gset, file=fn)
      cat("done!\n")
    }
    vsets[[i]] <- gset
  }
  names(vsets) <- gseids.val
  save(vsets, file=valfn)
}



