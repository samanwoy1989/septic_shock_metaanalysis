### R code from vignette source 'nibmgss_vignette.Rnw'

###################################################
### code chunk number 1: nibmgss_vignette.Rnw:25-26
###################################################
options(digits=2, width=80, continue=" ")


###################################################
### code chunk number 2: nibmgss_vignette.Rnw:47-48
###################################################
datapath <- "./"


###################################################
### code chunk number 3: nibmgss_vignette.Rnw:52-53
###################################################
source("Rcode/prelim.R")


###################################################
### code chunk number 4: nibmgss_vignette.Rnw:57-60
###################################################
library("ssnibmg") 
data("ss.eset")
allegs <- featureNames(ss.eset[[1]])


###################################################
### code chunk number 5: nibmgss_vignette.Rnw:64-65
###################################################
head(exprs(ss.eset[[1]])[1:5,1:5])


###################################################
### code chunk number 6: nibmgss_vignette.Rnw:68-69
###################################################
print(head(exprs(ss.eset[[1]])[1:5,1:5]))


###################################################
### code chunk number 7: nibmgss_vignette.Rnw:77-78
###################################################
source("Rcode/runGenelevelMetaanalysis.R")


###################################################
### code chunk number 8: nibmgss_vignette.Rnw:92-94
###################################################
data("genes.by.pathway")
data("pathways.list")


###################################################
### code chunk number 9: nibmgss_vignette.Rnw:118-119
###################################################
source("Rcode/runORA.R")


###################################################
### code chunk number 10: nibmgss_vignette.Rnw:123-124
###################################################
keggids.ora <- rownames(oraGene2Kegg.up[oraGene2Kegg.up[,1]<0.001,])


###################################################
### code chunk number 11: nibmgss_vignette.Rnw:129-130
###################################################
source("Rcode/runGSEA.R")


###################################################
### code chunk number 12: nibmgss_vignette.Rnw:134-135
###################################################
keggids.gsea <- names(which(gseaFp.up==0.0))


###################################################
### code chunk number 13: nibmgss_vignette.Rnw:139-140
###################################################
source("Rcode/runSPIA.R")


###################################################
### code chunk number 14: nibmgss_vignette.Rnw:144-145
###################################################
keggids.spia <- paste("hsa",top10ids, sep="")


###################################################
### code chunk number 15: nibmgss_vignette.Rnw:150-151 (eval = FALSE)
###################################################
## intersect(intersect(keggids.ora, keggids.gsea), keggids.spia)


###################################################
### code chunk number 16: nibmgss_vignette.Rnw:154-155
###################################################
print(intersect(intersect(keggids.ora, keggids.gsea), keggids.spia))


###################################################
### code chunk number 17: nibmgss_vignette.Rnw:162-164
###################################################
data("esetn.b")
data("esetn")


###################################################
### code chunk number 18: batchplot (eval = FALSE)
###################################################
## par(mfrow=c(1,2))
## boxplot(exprs(esetn.b), main="Before correction", las=2)
## boxplot(exprs(esetn), main="After correction", las=2)


###################################################
### code chunk number 19: nibmgss_vignette.Rnw:185-186
###################################################
source("Rcode/getPvalResample.R")


###################################################
### code chunk number 20: scplot
###################################################
source("Rcode/drawScatterplot.R")


###################################################
### code chunk number 21: nibmgss_vignette.Rnw:210-212
###################################################
which.lowp <- which(padjn<0.05)
print(length(which.lowp))


###################################################
### code chunk number 22: nibmgss_vignette.Rnw:217-220
###################################################
which.lowp.lowexp <- which(padjn<0.05 & xcon>log2(100) & xss-xcon>1)
egs.toshow <- names(which.lowp.lowexp)
print(length(egs.toshow))


###################################################
### code chunk number 23: nibmgss_vignette.Rnw:225-226
###################################################
print(as.character(unlist(mget(egs.toshow, org.Hs.egSYMBOL))))


###################################################
### code chunk number 24: boxplot
###################################################
source("Rcode/drawBoxplot25genes.R")


###################################################
### code chunk number 25: nibmgss_vignette.Rnw:246-248
###################################################
data(gset)
source("Rcode/analyse_validation_cohort.R")


###################################################
### code chunk number 26: pca
###################################################
pc<-prcomp(data.frame(t(dat)),scale=TRUE)
mycol<-rep(c("red"),length(dx))
mycol[which(dx=="Control")]<-c("green")
pca3d(pc, components = 1:3, radius=3, col = mycol,title = NULL, new = FALSE,axes.color = "grey", bg = "lightcyan")


###################################################
### code chunk number 27: survival (eval = FALSE)
###################################################
## source("Rcode/runSurvivalAnalysis.R")


###################################################
### code chunk number 28: proxplot
###################################################
prox.score <- read.delim(file="Metadata/egs25proximalityscore.txt", header=T, sep="\t")
o <- order(prox.score[,4])
egs25 <- as.character(prox.score[o,1])
lfc1 <- lfc[egs25]
keggsym25 <- as.character(prox.score[o,2])
mycol <- as.character(prox.score[o,5])
prox <- as.character(prox.score[o,"Position"])
par(mar=c(12,5.5,2,2))
plotdat <- split(lfc1, prox)
b <- boxplot(plotdat[c("Membrane","Intermediate","Nucleus")],
outline=F, col=c("red","pink","brown"), ylim=c(0,1.3),
axes=F, ylab="Relative Gene Expression (log2 scale)")
 box(col="gray")
 axis(2)
 mtext(text=c("Membrane","Intermediate","Nucleus"),
side=1, line=0.5, at=1:3, col="darkblue",
font=2, las=2, cex=2)


###################################################
### code chunk number 29: sessionInfo
###################################################
sessionInfo()


