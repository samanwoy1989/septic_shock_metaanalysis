o <- order(sapply(egs.toshow, function(eg)
	median(exprs(esetn[eg, is.ss])) - median(exprs(esetn[eg, is.con])) ))
egs.plot <- egs.toshow[o]
gsyms.plot <- gsyms.toshow[o]
plotdat <- list()
lcounter <- 1
for(i in 1:length(egs.toshow)) {
  eg <- egs.plot[i]
  sym <- gsyms.plot[i]
  condat <- exprs(esetn[eg, is.con])
  ssdat <- exprs(esetn[eg, is.ss])
  baseline <- median(condat)
  condat <- condat-baseline
  ssdat <- ssdat-baseline
  zero1 <- NA
  plotdat[[lcounter]] <- condat
  plotdat[[lcounter+1]] <- ssdat
  plotdat[[lcounter+2]] <- zero1
  lcounter <- lcounter+3
}
b <- boxplot(plotdat, col=c("darkgreen","red", "white"),
	pch=".", axes=F)
namepos <- seq(1, length(plotdat), by=3)+0.5
mtext(text=gsyms.plot, side=1, at=namepos, las=2, line=0.2, 
	font=2, col="darkblue")
abline(h=0, lty=2)
axis(2); box()
title(ylab="Log-gene expression (normalized)", line=2.5)
legend(legend=c("Control","SS"), fill=c("darkgreen","red"), bty="y", 
	x=30, y=4)

