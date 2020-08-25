#####################################################################
# Septic Shock data													#
# Gene-level meta-analysis											#
#####################################################################
# Compare Septic_Shock group with Zcontrol, run t-test
# Compute Z-stat from t-stat
# Save the results into a list							

# The correlation matrix between any pair of studies can be obtained 
# using the formula in the screenshot provided by SB.
# The R code implementing this formula is given below:

#Input: case overlap and control overlap matrices: nmat11, nmat10
# nmat10 = cross overlap matrix is optional
corr.mat.logit <- function(nmat11, nmat00, nmat10=NULL)
{
	k <- nrow(nmat11)
	if(is.null(nmat10)) nmat10 <- matrix(0, k, k)
        nmat01 <- t(nmat10)
	
	nvec0 <- diag(nmat00)
	nvec1 <- diag(nmat11)
	nvec <- nvec0 + nvec1
	
	neff <- (nvec0 * nvec1)/nvec
	rneff <- sqrt(neff)
	
	rmat <- (nmat11 * outer(1/nvec1, 1/nvec1) + nmat00 * outer(1/nvec0, 1/nvec0)
			 -nmat10 * outer(1/nvec1, 1/nvec0) - nmat01 * outer(1/nvec0, 1/nvec1)) * sqrt(outer(neff, neff))
	rmat
}


# read the overlap matrix; calculate pair-wise correlation
fcon <- "Metadata/overlapcontrol.txt"
fcase <- "Metadata/overlapcase.txt"
nmat11 <- as.matrix(read.table(file=fcase, header=T, sep="\t"))
nmat00 <- as.matrix(read.table(file=fcon, header=T, sep="\t"))
rmat <- corr.mat.logit(nmat11, nmat00)


# List all row-wise t-test result  ; for six studies
# Basic idea is to calculate z-statistic for each gene
# t-statistic is used in the formula
# log-fold change is calculated here and stored; used later

n6 <- vector(mode="numeric", length=6) # The denominator for a formula below
z6 <- data.frame(matrix(0, nrow=length(allegs), ncol=6)) # z-statistic
tstat6 <- data.frame(matrix(0, nrow=length(allegs), ncol=6)) # t-statistic
lfc6 <- data.frame(matrix(0, nrow=length(allegs), ncol=6)) # log-fold change
rownames(z6) <- allegs
rownames(tstat6) <- allegs 
rownames(lfc6) <- allegs
  
for(i in 1:length(ssids)) {
  # get the study id, platform info
  sid <- ssids[i]

  print(paste("Now analysing study ", sid))
  gset <- ss.eset[[i]]
  gset <- gset[allegs,]
  gset <- gset[, gset$Group%in%c("Zcontrol","Septic_Shock")]

  # log-fold change and t-statistic
  rtt <- rowttests(gset, "Group")
  lfc6[i] <- rtt$dm
  tstat <- rtt$statistic
  tstat6[i] <- tstat
  colnames(lfc6)[i] <- sid
  colnames(tstat6)[i] <- sid
    
  # Calculate the number 'n' for the formula 
  # z~ = (sqrt(n1)*z1 + sqrt(n2)*z2 + ...)/sqrt(n1+n2+...)
  ncon <- sum(gset$Group=="Zcontrol")
  nss <- sum(gset$Group=="Septic_Shock")
  n6[i] <- (ncon*nss)/(ncon+nss) # harmonic mean? 
  names(n6)[i] <- sid
  
  # convert t-stat to z-stat
  degf <- length(gset$Group)-2
  ptres <- pt(tstat, df=degf)
  z6[[i]] <- qnorm(ptres)
  colnames(z6)[i] <- sid
}

# Get the mean log-fold change
lfc <- apply(lfc6, 1, mean, na.rm=T)

# Calculate combined zstat
rneff <- sqrt(n6)

# Note by SB
# Numerator is as before:
#   numr <- sum(rneff * z)
# The current denominator 
#   denr <- sqrt(sum(neff))
# should be replaced by 
#   denr <- sqrt((matrix(rneff, nrow=1) %*% rmat %*% matrix(rneff, ncol=1))
# z.meta <- numr/denr

#denr <- sqrt(sum(n6)) # old denominator
#rmat <- diag(6)
denr <- sqrt((matrix(rneff, nrow=1) %*% rmat %*% matrix(rneff, ncol=1)))
z.meta <- apply(z6, 1, function(x) sum(x*rneff)/denr)
pvals <- 2*pnorm(abs(z.meta), lower.tail=F)
padj <- p.adjust(pvals, method="BH")

