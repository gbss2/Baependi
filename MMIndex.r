#==========================================
# GBSS Project
#==========================================
#
# (/u0254) Copyleft 2015, by GBSS and Contributors.
#
# 
# -----------------
# MMIndex.r
# -----------------
# GNU GPL 2015, by GBSS and Contributors.
#
# Original Author: Giordano Bruno Soares-Souza
# Contributor(s): 
# Updated by:
#
# Command line: function(x) MMIndex(<data>)
# Dependencies: R, library MASS
# Description: 
# 
# 
#
################################################################################

MMindex <- function(x){

	if(!require(MASS)) { install.packages("MASS"); require(MASS)}

	objname=deparse(substitute(x))
	y = as.data.frame(matrix(0,dim(x)[1],39))
	colnames(y) = c("IID","FID","PAT","MAT","SEX","AGE","ETHNICITY","PHENO","GENO","mTE","mBDE","mBEE","mBE","mE","mTM","mBDM","mBEM","mBM","mM","mTR","mBDR","mBER","mBR","mR","mTG","mBDG","mBEG","mBG","mG","mTB","mBDB","mBEB","mBB","mB","V6","ENDERECO","EUR","AFR","NAT")
	x = data.matrix(x)
	y[,c(1:9,35:39)]= x[,c(1:9,55:59)]
for (j in 1 : dim(x)[1]){
	y[j,10] = try (sqrt(as.numeric(huber(x[j,c(10,12,14)])[1])),TRUE) #mTE
	y[j,11] = try (sqrt(as.numeric(huber(x[j,c(16,18,20)])[1])),TRUE) #mBDE
	y[j,12] = try (sqrt(as.numeric(huber(x[j,c(22,24,26)])[1])),TRUE) #mBEE
	y[j,13] = try (sqrt(as.numeric(huber(x[j,c(10,12,14,16,18,20)])[1])),TRUE) #mBE
	y[j,14] = try (sqrt(as.numeric(huber(x[j,c(10,12,14,16,18,20,22,24,26)])[1])),TRUE) #mE

	y[j,15] = try (sqrt(as.numeric(huber(x[j,c(11,13,15)])[1])),TRUE) #mTM
	y[j,16] = try (sqrt(as.numeric(huber(x[j,c(17,19,21)])[1])),TRUE) #mBDM
	y[j,17] = try (sqrt(as.numeric(huber(x[j,c(23,25,27)])[1])),TRUE) #mBEM
	y[j,18] = try (sqrt(as.numeric(huber(x[j,c(11,13,15,17,19,21)])[1])),TRUE) #mBM
	y[j,19] = try (sqrt(as.numeric(huber(x[j,c(11,13,15,17,19,21,23,25,27)])[1])),TRUE) #mM

	y[j,20] = try (sqrt(as.numeric(huber(x[j,c(28,31,34)])[1])),TRUE) #mTR
	y[j,21] = try (sqrt(as.numeric(huber(x[j,c(37,40,43)])[1])),TRUE) #mBDR
	y[j,22] = try (sqrt(as.numeric(huber(x[j,c(46,49,52)])[1])),TRUE) #mBER
	y[j,23] = try (sqrt(as.numeric(huber(x[j,c(28,31,34,37,40,43)])[1])),TRUE) #mBR
	y[j,24] = try (sqrt(as.numeric(huber(x[j,c(28,31,34,37,40,43,46,49,52)])[1])),TRUE) #mR
	
	y[j,25] = try (sqrt(as.numeric(huber(x[j,c(29,32,35)])[1])),TRUE) #mTG
	y[j,26] = try (sqrt(as.numeric(huber(x[j,c(38,41,44)])[1])),TRUE) #mBDG
	y[j,27] = try (sqrt(as.numeric(huber(x[j,c(47,50,53)])[1])),TRUE) #mBEG
	y[j,28] = try (sqrt(as.numeric(huber(x[j,c(29,32,35,38,41,44)])[1])),TRUE)  #mBG
	y[j,29] = try (sqrt(as.numeric(huber(x[j,c(29,32,35,38,41,44,47,50,53)])[1])),TRUE) #mG

	y[j,30] = try (sqrt(as.numeric(huber(x[j,c(30,33,36)])[1])),TRUE) #mTB
	y[j,31] = try (sqrt(as.numeric(huber(x[j,c(39,42,45)])[1])),TRUE) #mBDB
	y[j,32] = try (sqrt(as.numeric(huber(x[j,c(48,51,54)])[1])),TRUE) #mBEB
	y[j,33] = try (sqrt(as.numeric(huber(x[j,c(47,50,53,39,42,45)])[1])),TRUE) #mBB
	y[j,34] = try (sqrt(as.numeric(huber(x[j,c(47,50,53,39,42,45,48,51,54)])[1])),TRUE) #mB
	
	y[j,] <- ifelse(grepl("Error",y[j,]),NA,y[j,])
	}
	
	write.table(y,file=paste0(objname,"MMIndex.out"),col.names=T,row.names=F,quote=F,sep="\t")
	
}