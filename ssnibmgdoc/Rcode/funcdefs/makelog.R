#################################################
# A function for transforming the		#
# expressionset to log base 2, minimum set to 0	#
# First check if the maxval of data > 20	#
# Then shift the values to make min=1		#
# Then apply log2				#
# otherwise do nothing				#
#################################################
makelog <- function(eset) {
  gexp <- exprs(eset)
  if(max(gexp, na.rm=T)>20) {
    cat("Log transforming gene expression data ...")
    gexp <- gexp+1-min(gexp, na.rm=T)
    gexp[!is.na(gexp)] <- log2(gexp[!is.na(gexp)])
    exprs(eset) <- gexp
    cat("\n")
  } else {
    cat("Looks like this expression set is already log-transformed; doing nothing.\n")
  }
  eset
}

