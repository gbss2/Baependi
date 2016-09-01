#==========================================
# GBSS Project
#==========================================
#
# (/u0254) Copyleft 2016, by GBSS and Contributors.
#
# 
# -----------------
# relAdmix.r
# -----------------
# GNU GPL 2016, by GBSS and Contributors.
#
# Original Author: Giordano Bruno Soares-Souza
# Contributor(s): 
# Updated by:
#
# Command line: Rscript relAdmix(<fam,Q file,pops[1..n]>)
# Dependencies: R, library gtools
# Description: Compare admixture proportions between parents and offspring
# Future development:
# 	1) Write mean and/or range on ploting area
#	2) Improve summary output
#
################################################################################

if(!require(gtools)) { install.packages("gtools"); require(gtools)}
args = commandArgs(trailingOnly=TRUE)
#objnameX=deparse(substitute(args[1]))
fam = read.table(args[1],head=F,as.is=T)
qFile = read.table(args[2],head=F,as.is=T)
pops = c(args[3:(ncol(qFile)+2)])

plotFile = gsub("fam", "relAdmix.pdf", args[1])
tblFile = gsub("fam", "_relAdmix.txt", args[1])
sumFile = gsub("fam", "_summary_relAdmix.txt", args[1])
                                      
print(dim(fam))
colnames(fam)=c("FID","IID","PAT","MAT","SEX","PHENO")

x = cbind(fam,qFile)
x = subset(x,FID!=IID)
y = as.data.frame(matrix(NA,dim(x)[1],(4+(6*ncol(qFile)))))
colnames(y)[1:(4+(6*ncol(qFile)))] = c("FID","IID","PAT","MAT",rep(pops,6))
#mom = as.data.frame(matrix(NA,dim(x)[1],(3+ncol(qFile)))
#both = as.data.frame(matrix(NA,dim(x)[1],(3+ncol(qFile)))
x[is.na(x)] <- 0

	for (i in 1 : dim(x)[1]){
		idx = y[i,2] = x[i,2]
#		print(idx)
		y[i,1:(4+ncol(qFile))] = x[i,c(1:4,7:ncol(x))]
		idP = y[i,3]
		idM = y[i,4]
		y[i,(ncol(qFile)+5):(2*ncol(qFile)+4)] = if(idP %in% x[,2]){x[which(x[,2]==idP),7:ncol(x)]} else {rep(NA,ncol(qFile))}
		y[i,(2*ncol(qFile)+5):(3*ncol(qFile)+4)] = if(idM %in% x[,2]){x[which(x[,2]==idM),7:ncol(x)]} else {rep(NA,ncol(qFile))}
		
		if(idP %in% x[,2]){x[which(x[,2]==idP),7:ncol(x)]} else {rep(NA,ncol(qFile))}
		
#		y[i,(3*ncol(qFile)+5):(4*ncol(qFile)+4)] = (y[i,(ncol(qFile)+5):(2*ncol(qFile)+4)]/2) - y[i,5:(4+ncol(tblFile))]
#		y[i,(4*ncol(qFile)+5):(5*ncol(qFile)+4)] = (y[i,(2*ncol(qFile)+5):(3*ncol(qFile)+4)]/2) - y[i,5:(4+ncol(tblFile))]

		y[i,(3*ncol(qFile)+5):(4*ncol(qFile)+4)] = ifelse(y[i,5:(4+ncol(qFile))] > (y[i,(ncol(qFile)+5):(2*ncol(qFile)+4)]/2)+0.5 | y[i,5:(4+ncol(qFile))] < y[i,(ncol(qFile)+5):(2*ncol(qFile)+4)]/2, y[i,5:(4+ncol(qFile))]/2 - y[i,(ncol(qFile)+5):(2*ncol(qFile)+4)]/2,NA)
		y[i,(4*ncol(qFile)+5):(5*ncol(qFile)+4)] = ifelse(y[i,5:(4+ncol(qFile))] > (y[i,(2*ncol(qFile)+5):(3*ncol(qFile)+4)]/2)+0.5 | y[i,5:(4+ncol(qFile))] < y[i,(2*ncol(qFile)+5):(3*ncol(qFile)+4)]/2, y[i,5:(4+ncol(qFile))]/2 - y[i,(2*ncol(qFile)+5):(3*ncol(qFile)+4)]/2,NA)
		y[i,(5*ncol(qFile)+5):(6*ncol(qFile)+4)] = y[i,5:(4+ncol(qFile))] - ((y[i,(ncol(qFile)+5):(2*ncol(qFile)+4)]/2) + (y[i,(2*ncol(qFile)+5):(3*ncol(qFile)+4)]/2))

		}
#		print(dim(y))
		pdf(plotFile,width=10,height=4)
		par(mfrow=c(1,pops),mar=c(3,4,1,2))
#		bplot = boxplot(y[,(3*ncol(qFile)+5):(4*ncol(qFile)+4)],ylab=expression(delta ~ "Father"))
		boxplot(y[,(3*ncol(qFile)+5):(4*ncol(qFile)+4)],ylab=expression(delta ~ "Father"))
#		text(x=bplot,y=-2,label=paste0("\nMean = ",mean(y[,(3*ncol(qFile)+5):(4*ncol(qFile)+4)])),pos=1 cex.lab=0.8)
#		bplot = boxplot(y[,(4*ncol(qFile)+5):(5*ncol(qFile)+4)],ylab=expression(delta ~ "Mother"))
		boxplot(y[,(4*ncol(qFile)+5):(5*ncol(qFile)+4)],ylab=expression(delta ~ "Mother"))
#		text(x=bplot,y=-2,label=paste0("\nMean = ",mean(y[,(4*ncol(qFile)+5):(5*ncol(qFile)+4)])),pos=1 cex.lab=0.8)
#		bplot = boxplot(y[,(5*ncol(qFile)+5):(6*ncol(qFile)+4)],ylab=expression(delta ~ "Parents"))
		boxplot(y[,(5*ncol(qFile)+5):(6*ncol(qFile)+4)],ylab=expression(delta ~ "Parents"))
#		text(x=bplot,y=-2,label=paste0("\nMean = ",mean(y[,(5*ncol(qFile)+5):(6*ncol(qFile)+4))),pos=1 cex.lab=0.8)
		dev.off()
		
		yRow = as.data.frame(matrix(NA,1,(4+(6*ncol(qFile)))))
#		yRow[,1:4] = c("FID","IID","PAT","MAT")
#		yRow[,5:6*ncol(qFile)] = mean(y[,5:(6*ncol(qFile))],na.rm=T)
#		y = rbind(y,yRow)
		write.table(y,file=tblFile,col.names=T,row.names=F,quote=F)
		capture.output(summary(y[,5:(6*ncol(qFile)+4)]),file=sumFile)
