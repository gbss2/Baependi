#==========================================
# GBSS Project
#==========================================
#
# (/u0254) Copyleft 2016, by GBSS and Contributors.
#
# 
# -----------------
# fampop.r
# -----------------
# GNU GPL 2016, by GBSS and Contributors.
#
# Original Author: Giordano Bruno Soares-Souza
# Contributor(s): 
# Updated by:
#
# Command line: Rscript fampop.r <pop_file> <fam>
# Dependencies: R
# Description: Creates pop file input for Admixture supervised analysis and 
# an auxiliary file called <input_filename>.fampop that extends fam file with
# extra pop columns: FID IID PAT MAT SEX PHENO POP_1(subpop) POP_2(large pop).
# The script needs an extra file cointaing information about pop affiliation in
# the following format with 3 fields -> IID POP_1(subpop) POP_2(pop or continent)
# Future Developments:
#
################################################################################

args = commandArgs(trailingOnly=TRUE)

popscon = read.table(args[1],head=F,as.is=T,sep="\t")

fam = read.table(args[2],head=F,as.is=T)

colnames(popscon) =c('FID','IID','POP','CON')
colnames(fam) =c('FID','IID','PAT','MAT','SEX','PHENO')

fampop = merge(fam,popscon,by="IID",all.x=T,sort=F)

fampopfile = gsub("fam", "fampop", args[2])
popfile = gsub("fam", "pop", args[2])

write.table(fampop[,c(2,1,3,4,5,6,8,9)],fampopfile,quote=F,row.names=F,col.names=F,sep="\t")
write.table(fampop[,c(9)],popfile,quote=F,row.names=F,col.names=F,sep="\t")
