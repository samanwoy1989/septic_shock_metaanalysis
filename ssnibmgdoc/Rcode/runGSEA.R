#########################################################################
# Septic Shock data														#
# Gene Set Enrichment Analysis											#
#########################################################################
# Create the incidence matrix
pathway.names <- names(genes.by.pathway)
pathway.genes <- unique(unlist(genes.by.pathway))
numPathways <- length(pathway.names)
numGenes <- length(pathway.genes)
Am <- matrix(0, nrow=numPathways, ncol=numGenes)
rownames(Am) <- pathway.names
colnames(Am) <- pathway.genes
for(i in 1:length(genes.by.pathway)) {
  Am[i,genes.by.pathway[[i]]] <- 1
}

# Reduce the incidence matrix by removing all gene sets that have 10 genes
# or fewer
selectedRows = (rowSums(Am)>10)
Am2 = Am[selectedRows, ]

############################
# Permutation-based
############################

#gseap.fn <- "Result/gseap_01Aug2014_nperm.10000.rda"
#gseap.fn <- "gseap_28April2016_nperm.10000.rda"
gseap.fn <- "gseap_nperm.10000_2016.05.02.rda"
if(file.exists(file=gseap.fn)) {
  load(file=gseap.fn)
} else {
  gseap <- list()
  for(i in 1:length(ss.eset)) {
    cat(names(ss.eset)[i])
    nsF <- ss.eset[[i]]
    set.seed(123)
    NPERM <- 10000
    Am2new <- Am2[, intersect(featureNames(nsF), colnames(Am2))]
    gseap[[i]] <- gseattperm(nsF, nsF$Group, Am2new, NPERM)
    names(gseap)[i] <- names(ss.eset)[i]
    cat("\n")
  }
  save(gseap, file=gseap.fn)
}

# Get the lowp and highp values
gseap.dn <- sapply(gseap, function(x) x[,1])
gseap.up <- sapply(gseap, function(x) x[,2])

# Calculate Fisher product
# Up-regulated genes
gseaFtest.up <- apply(gseap.up, 1, Fisher.test) # get the Fisher-product for all paths
gseaFp.up  <- gseaFtest.up["p.value",] # Extract the p-value for each pathway
gseaXsq.up <- gseaFtest.up["Xsq",] # Extract the Xsq for each pathway

# Down-regulated genes
gseaFtest.dn <- apply(gseap.dn, 1, Fisher.test) # get the Fisher-product for all paths
gseaFp.dn  <- gseaFtest.dn["p.value",] # Extract the p-value for each pathway
gseaXsq.dn <- gseaFtest.dn["Xsq",] # Extract the Xsq for each pathway

# correct for multiple testing
gseaFp.up <- p.adjust(gseaFp.up, method="bonferroni")
gseaFp.dn <- p.adjust(gseaFp.dn, method="bonferroni")


