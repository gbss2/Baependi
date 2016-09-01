#==========================================
#
# (/u0254) Copyleft 2016, by GBSS and Contributors.
#
# 
# -----------------
# reapPlot.r
# -----------------
# GNU GPL 2016, by GBSS and Contributors.
#
# Original Author: Giordano Bruno Soares-Souza
# Contributor(s): 
# Updated by:
#
# Command line: function(x) reapPlot(<REAP_pairs_relatedness file>)
# Dependencies: R, library gtools
# Description: Plot Reap IBD 0 vs Kinship Coefficient Estimates
# Future Developments:
# 	1) Color gradient instead of arbitrary color limits
#
################################################################################

args = commandArgs(trailingOnly=TRUE)
reap = read.table(args[1],head=T,as.is=T)

plotFile = gsub(".txt", "_plot.png", args[1])
png(plotFile)
color = ifelse(reap$KINCOEF >= 0.4,'black',ifelse(reap$KINCOEF > 0.175 & reap$KINCOEF < 0.4 & reap$IBD0_PROB < 0.1,'red',ifelse(reap$KINCOEF >= 0.175 & reap$KINCOEF < 0.4 & reap$IBD0_PROB >= 0.1,'blue',ifelse(reap$KINCOEF <= 0.175 & reap$KINCOEF > 0.075 & reap$IBD0_PROB >= 0.15,'green',ifelse(reap$KINCOEF <= 0.075 & reap$IBD0_PROB >= 0.15,'purple','yellow')))))
plot(reap$IBD0_PROB,reap$KINCOEF,ylab="REAP Kinship - Coefficient Estimates",xlab="REAP IBD 0 Estimates",col=color,pch=19)
legend("topright",pch=19,legend=c("Twins","Parent-Offspring","Full Siblings","Second Degree","Third Degree or less","?"),col=c("black","red","blue","green","purple","yellow"))
dev.off()

