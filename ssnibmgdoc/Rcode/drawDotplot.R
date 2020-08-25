#####################################################################
# Dot plot
# Show that there are more up-regulated genes compared to down-regulated
#####################################################################
# Dot-plot of log-fold changes
par(mar=c(6,5,4,2))
plot(x=1:6, y=1:6, type="n", xlim=c(0.5,6.5), ylim=c(-2,2), axes=F, xlab="",
	ylab="Log-fold change in gene expression")
for(i in 1:ncol(lfc6)) {
 y <- lfc6[,i]
 cols <- rep("green", nrow(lfc6))
 cols[y>=0] <- "red"
 x <- seq(i-0.25, i+0.25, length.out=nrow(lfc6))
 points(x=x, y=y, col=cols, pch=16, cex=0.5)
}

axis(1, at=1:6, labels=colnames(lfc6), tick=F, las=2)
axis(2, cex.axis=1); box()
abline(h=c(0,-1,1), lwd=c(2,3,3), lty=c(1,2,2), col=c("gray","darkblue","darkblue"))
#dev.copy(device=png, file="Result/Figs_Tables/Fig_2_dotplot_logfc_gene_expression.png")
#dev.off()

