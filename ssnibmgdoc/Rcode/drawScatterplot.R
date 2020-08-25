egs.hsa04380 <- genes.by.pathway$hsa04380 
egs.sel <- intersect(egs.hsa04380, featureNames(esetn)) 
egs.sel <- intersect(egs.sel, allegs) 
is.con <- esetn$Group=="Control"
is.ss <- esetn$Group=="Septic Shock"
pvalsn <- apply(exprs(esetn[egs.sel,]), 1, function(x) {
		t.test(x[is.con], x[is.ss], alternative="two.sided")$p.value
	})
padjn <- p.adjust(pvalsn, method="fdr")
xcon <- rowMeans(exprs(esetn[egs.sel, is.con]))
xss  <- rowMeans(exprs(esetn[egs.sel, is.ss]))
which.lowp <- which(padjn<0.05 & xcon>log2(100) & xss-xcon>1)
egs.toshow <- names(which.lowp)
gsyms.toshow <- as.character(unlist(mget(egs.toshow, org.Hs.egSYMBOL)))
egs.noshow <- setdiff(egs.sel, egs.toshow)
plot(x=xcon, y=xss, xlab="Control (Log 2 scale)", 
	ylab="Septic Shock (Log 2 scale)",
	xlim=c(2,max(c(xcon,xss))+0.5),
	ylim=c(2,max(c(xcon,xss))+0.5),
	pch=16, col="darkblue", type="n")
abline(0,1, col="gray", lty=1, lwd=2)
points(x=xcon[egs.noshow], y=xss[egs.noshow], pch=16, col="gray80")
text(x=jitter(xcon[egs.toshow], factor=100), y=jitter(xss[egs.toshow], factor=100),
	labels=gsyms.toshow,
	cex=0.5, font=4, col="red")
#text(x=5, y=13, cex=1.2, font=3)
