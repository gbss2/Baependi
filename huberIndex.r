#==========================================
# GBSS Project
#==========================================
#
# (/u0254) Copyleft 2015, by GBSS and Contributors.
#
# 
# -----------------
# huberIndex.r
# -----------------
# GNU GPL 2015, by GBSS and Contributors.
#
# Original Author: Giordano Bruno Soares-Souza
# Contributor(s): 
# Updated by:
#
# Command line: function(x) huberIndex(<data>)
# Dependencies: R, library MASS
# Description: 
# 
# 
#
################################################################################

huberIndex <- function(x){

	if(!require(MASS)) { install.packages("MASS"); require(MASS)}

	objname=deparse(substitute(x))
	y = as.data.frame(matrix(0,dim(x)[1],39))
	colnames(y) = c("IID","FID","PAT","MAT","SEX","AGE","ETHNICITY","PHENO","GENO","hTE","hBDE","hBEE","hBE","hE","hTM","hBDM","hBEM","hBM","hM","hTR","hBDR","hBER","hBR","hR","hTG","hBDG","hBEG","hBG","hG","hTB","hBDB","hBEB","hBB","hB","V6","ENDERECO","EUR","AFR","NAT")
	x = data.matrix(x)
	y[,c(1:9,35:39)]= x[,c(1:9,55:59)]
for (j in 1 : dim(x)[1]){
	y[j,10] = try ((as.numeric(huber(x[j,c(10,12,14)])[1])),TRUE)
	y[j,11] = try ((as.numeric(huber(x[j,c(16,18,20)])[1])),TRUE)
	y[j,12] = try ((as.numeric(huber(x[j,c(22,24,26)])[1])),TRUE)
	y[j,13] = try ((as.numeric(huber(x[j,c(16,18,20,22,24,26)])[1])),TRUE)
	y[j,14] = try ((as.numeric(huber(x[j,c(10,12,14,16,18,20,22,24,26)])[1])),TRUE)

	y[j,15] = try ((as.numeric(huber(x[j,c(11,13,15)])[1])),TRUE)
	y[j,16] = try ((as.numeric(huber(x[j,c(17,19,21)])[1])),TRUE)
	y[j,17] = try ((as.numeric(huber(x[j,c(23,25,27)])[1])),TRUE)
	y[j,18] = try ((as.numeric(huber(x[j,c(23,25,27,17,19,21)])[1])),TRUE)
	y[j,19] = try ((as.numeric(huber(x[j,c(11,13,15,17,19,21,23,25,27)])[1])),TRUE)

	y[j,20] = try ((as.numeric(huber(x[j,c(28,31,34)])[1])),TRUE)
	y[j,21] = try ((as.numeric(huber(x[j,c(37,40,43)])[1])),TRUE)
	y[j,22] = try ((as.numeric(huber(x[j,c(46,49,52)])[1])),TRUE)
	y[j,23] = try ((as.numeric(huber(x[j,c(46,49,52,37,40,43)])[1])),TRUE)
	y[j,24] = try ((as.numeric(huber(x[j,c(28,31,34,37,40,43,46,49,52)])[1])),TRUE)
	
	y[j,25] = try ((as.numeric(huber(x[j,c(29,32,35)])[1])),TRUE)
	y[j,26] = try ((as.numeric(huber(x[j,c(38,41,44)])[1])),TRUE)
	y[j,27] = try ((as.numeric(huber(x[j,c(47,50,53)])[1])),TRUE)
	y[j,28] = try ((as.numeric(huber(x[j,c(47,50,53,38,41,44)])[1])),TRUE)
	y[j,29] = try ((as.numeric(huber(x[j,c(29,32,35,38,41,44,47,50,53)])[1])),TRUE)

	y[j,30] = try ((as.numeric(huber(x[j,c(30,33,36)])[1])),TRUE)
	y[j,31] = try ((as.numeric(huber(x[j,c(39,42,45)])[1])),TRUE)
	y[j,32] = try ((as.numeric(huber(x[j,c(48,51,54)])[1])),TRUE)
	y[j,33] = try ((as.numeric(huber(x[j,c(48,51,54,39,42,45)])[1])),TRUE)
	y[j,34] = try ((as.numeric(huber(x[j,c(30,33,36,39,42,45,48,51,54)])[1])),TRUE)
	
	y[j,] <- ifelse(grepl("Error",y[j,]),NA,y[j,])
	}
	
	write.table(y,file=paste0(objname,"huberIndex.out"),col.names=T,row.names=F,quote=F,sep="\t")
	
}