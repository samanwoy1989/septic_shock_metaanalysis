#####################################################################
# Over-representation analysis
#   Get the pathways enriched with the up-regulated genes
#   There are no down-regulated genes selected for equivalent 
#   thresholds of p-values and log-fold change
# Over-representation analysis  						
# approach based on KEGGREST
# Some help taken from Pauls's contribution on the bioc mailing list
# https://stat.ethz.ch/pipermail/bioconductor/2013-September/054902.html
###########################################################################
print.noquote("Over-representation analysis: using KEGGREST")
# Up-regulated pathways from gene-level meta-analysis
# KEGG approach
up_egs.v2 <- names(which(padj<0.01 & lfc>=1))
oraGene2Kegg.up <-
  t(sapply(genes.by.pathway, hyperg.test, up_egs.v2, allegs))
o <- order(as.numeric(oraGene2Kegg.up[,"p"]))  
oraGene2Kegg.up <- data.frame(oraGene2Kegg.up[o,])
oraGene2Kegg.up$Name <-  
  as.character(pathways.list[paste("path:",rownames(oraGene2Kegg.up),sep="")])
print(oraGene2Kegg.up[oraGene2Kegg.up[,1]<0.001,])
#                    p     odds  expected
# hsa04610 5.865805e-06 11.51028 0.6921457
# hsa04380 5.673352e-05 6.590726  1.324105
# hsa05202 0.0004608154 4.767788  1.795566
#                                                                   Name
# hsa04610     Complement and coagulation cascades - Homo sapiens (human)
# hsa04380              Osteoclast differentiation - Homo sapiens (human)
# hsa05202 Transcriptional misregulation in cancer - Homo sapiens (human)
dn_egs.v2 <- names(which(padj<0.01 & lfc<=(-1.0)))
oraGene2Kegg.dn <-
  t(sapply(genes.by.pathway, hyperg.test, dn_egs.v2, allegs))
o <- order(as.numeric(oraGene2Kegg.dn[,"p"]))  
oraGene2Kegg.dn <- data.frame(oraGene2Kegg.dn[o,])
oraGene2Kegg.dn$Name <-  
  as.character(pathways.list[paste("path:",rownames(oraGene2Kegg.dn),sep="")])
oraGene2Kegg.dn[oraGene2Kegg.dn[,1]<0.001,]
# [1] p        odds     expected Name    
# <0 rows> (or 0-length row.names)

