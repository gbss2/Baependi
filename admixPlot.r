#==========================================
# GBSS Project
#==========================================
#
# (/u0254) Copyleft 2015, by GBSS and Contributors.
#
# 
# -----------------
# admixPlot.r
# -----------------
# GNU GPL 2015, by GBSS and Contributors.
#
# Original Author: Giordano Bruno Soares-Souza
# Contributor(s): 
# Updated by:
#
# Command line: Rscript admixPlot.r <k> <data>
# Dependencies: R
# Description: Admixture barplot Script
# Future development:
#	1) Read K from file name
# 	2) Write populations labels
#
################################################################################

if(!require(gridExtra)) { install.packages("gridExtra"); require(gridExtra)}


args = commandArgs(trailingOnly=TRUE)
k = args[1]
tblS1=read.table(args[2],head=F, as.is=T)
#popfile = paste0(strsplit(args[2],"[.]")[1],".pop")
#popfile = gsub("3.Q", "pop", args[2])
popfile = gsub(paste0(k,".Q"), "fampop", args[2])
#pdf(paste0(strsplit(args[2],"[.]")[1],".pdf"))
pdf(gsub(".Q",".pdf",args[2]),width=12,height=6)
pop=read.table(popfile,head=F, as.is=T)
colnames(pop) =c('FID','IID','MAT','PAT','SEX','PHENO',"POP","CON")
#colnames(pop) = c("pop")
data2 = cbind(pop[,c(8,7)],tblS1,pop[,2])
data2 = data2[do.call(order,as.list(data2)),]
write.table(data2[,c(ncol(data2),1:(ncol(data2)-1))],gsub(paste0(k,".Q"), paste0(k,"_data2.txt"), args[2]),col.names=T,quote=F,row.names=F)
barplot(t(as.matrix(data2[,3:(ncol(data2)-1)])), col=rainbow(k), border=NA,xaxt="n",axes=F)
plot.new()
admxtable = aggregate(data2[,3:(ncol(data2)-1)],by=list(data2[,1]),FUN=function(x) mean(as.numeric(as.character(x))))
tt = ttheme_default(colhead=list(fg_params=list(parse=TRUE),bg_params=list(fill=c("white",rainbow(k)))))
grid.table(admxtable,theme=tt)
dev.off()

