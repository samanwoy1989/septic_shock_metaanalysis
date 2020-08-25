
################################################################################
# JPM suggests check correlation with severity, mortality, NE requirement, etc.
#				GSE13904	Day of SS # not used
#				GSE26378	Age, Survival
#				GSE26440	Age, Survival, Sub-group
#				GSE4607		survival, day of SS, PRISM, steroid, organism,
#							infection site, Sub-group
#				GSE8121		Day of SS
#				GSE9692		Gender, survival, organism 
################################################################################

# get expression data for the 25 significant genes
gexp <- exprs(ss.eset[[1]][egs.toshow,])
for(i in 2:length(ss.eset)) {
  gexp <- cbind(gexp, exprs(ss.eset[[i]][egs.toshow,]))
}

# read phenodata table
pdata <- read.table("Metadata/phenodata.csv", header=T, sep="\t", row.names=1)
#pdata <- read.table("./Metadata/phenodata.csv", header=T, sep="\t", row.names=1)# SM
# create subgroup, outcome, prism vectors
for(i in 1:length(ss.eset)) {
  print(i)
  gseid <- names(ss.eset)[i]
  eset <- ss.eset[[gseid]]



  # get the outcome information
  sel <- as.character(pdata[gseid, "survival.column"])
  if(is.na(sel)) {
    psurv <- rep(NA, ncol(eset))
  } else {
    surv <- pData(eset)[, sel]
    sel <- as.integer(pdata[gseid, "survival.index"])
    psurv <- as.character(sapply(as.character(surv), function(x) 
      strsplit(x, split=": ")[[1]][sel]))
    psurv[eset$Group=="Zcontrol"] <- "Control"
    psurv <- toupper(psurv)
  }

  # get the PRISM information
  sel <- as.character(pdata[gseid, "PRISM.column"])
  if(is.na(sel)) {
    prismdat <- rep(NA, ncol(eset))
  } else {
    prismdata <- pData(eset)[, sel]
    sel <- as.integer(pdata[gseid, "PRISM.index"])
    prismdat <- as.character(sapply(as.character(prismdata), function(x) 
      strsplit(x, split=": ")[[1]][sel]))
    prismdat[eset$Group=="Zcontrol"] <- "Control"
    prismdat <- toupper(prismdat)
  }


  # add to the subgroup and outcome vector
  if(i == 1) {   
    outcome <- as.character(psurv)
    names(outcome) <- sampleNames(eset)
    
    prism <- as.character(prismdat)
    names(prism) <- sampleNames(eset)
    
  } else {
    
    outcome.new <- as.character(psurv)
    names(outcome.new) <- sampleNames(eset)  
    outcome <- c(outcome, outcome.new)
    
    prism.new <- as.character(prismdat)
    names(prism.new) <- sampleNames(eset)
    prism <- c(prism, prism.new)
  }
  outcome[outcome%in%c("NON SURVIVOR", "NONSURVIVOR")] <- "NS"
  outcome[outcome=="SURVIVOR"] <- "S"
  outcome[outcome=="CONTROL"] <- "CON"
  
  prism[prism=="CONTROL"] <- "0"
  prism[prism=="N/A"] <- NA  
}

################################################################################
# box plot of data versus outcome
################################################################################
y <- outcome
xo <- apply(gexp, 1, split, outcome)
names(xo) <- rownames(gexp)
# to record median_nonsurv_x - median_surv_x; useful for ordering the boxplots
dsurv <- sapply(xo, function(x) median(x$NS)-median(x$S))
o <- order(dsurv)
xo <- xo[o]

xo <- lapply(xo, function(x) {# normalize the control median to zero
  xbar <- median(x$CON)
  for(i in 1:length(x)) {
    x[[i]] <- x[[i]]-xbar
  }
  x})
xout <- list()
j <- 1
pvals <- rep(NA, length(xo))
names(pvals) <- names(xo)
for(i in 1:length(xo)) {
  xout[[j]] <- xo[[i]]$CON
  xout[[j+1]] <- xo[[i]]$NS
  xout[[j+2]] <- xo[[i]]$S
  pvals[i] <- t.test(xo[[i]]$NS, xo[[i]]$S, alternative="less")$p.value
  xout[[j+3]] <- NA
  
  eg <- names(xo)[i]
  gsym <- as.character(get(eg, org.Hs.egSYMBOL))
  names(xout)[j:(j+3)] <- c("",paste(eg,gsym,sep=" | "),"","")
  j <- j+4
}
par(mar=c(8,4,4,2))
b <- boxplot(xout, outline=F,col=c("white","red","blue", "white"),
	axes=F, border=c("white","pink","lightblue", "white"),
	main="Effect of outcome on gene expression", 
	ylab="Gene expression (log-intensity)")
legend(x=30, y=3, legend=c("Non-survivor","Survivor"), bty="n",
	fill=c("red","blue"))
abline(h=0, lty=2, col="gray50")
mtext(text=names(xout), side=1, las=2, line=0.5, at=1:length(xout))
xpos <- seq(2, 100, by=4)
ypos <- -1 # b$stats[5,xpos]+0.3
pstr <- c("","*")[as.integer(pvals<0.1)+1]
text(x=xpos, y=ypos, labels=pstr, cex=2, col="blue")
axis(2); box()

#dev.copy(device=pdf, file="/boxplot_gexp_surv.pdf", width=12); dev.off()


