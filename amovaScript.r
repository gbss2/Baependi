#==========================================
# LDGH Project
#==========================================
#
# (/u0254) Copyleft 2012, by LDGH and Contributors.
#
# 
# -----------------
# (1-3)levels Hierfstat script
# -----------------
# GNU GPL 2012, by LDGH and Contributors.
#
# Original Author: Wagner C S Magalhães
# Contributor(s): Giordano Bruno Soares-Souza
# Updated by: Giordano Bruno Soares-Souza
#
# Command line: R: source("hierfstat_script.R") || shell: nohup R CMD BATCH --save ./hierfstat_script.R &
# Dependencies: R packages: hierfstat, R.utils
# Description: Performs hierarquical (1-3) AMOVA statistics for loci
#
################################################################################

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  writeLines("Input provided")
}

if(!require(hierfstat)){install.packages("hierfstat", repos="http://cran.rstudio.com/");require(hierfstat)}
if(!require(foreach)){install.packages("foreach", repos="http://cran.rstudio.com/");require(foreach)}
if(!require(doMC)){install.packages("doMC", repos="http://cran.rstudio.com/");require(doMC)}
if(!require(data.table)){install.packages("data.table", repos="http://cran.rstudio.com/");require(data.table)}
if(!require(R.utils)){install.packages("R.utils", repos="http://cran.rstudio.com/");require(R.utils)}
if(!require(iterators)){install.packages("iterators", repos="http://cran.rstudio.com/");require(iterators)}

registerDoMC(2)


######## SCRIPT HEADER AND FILE CALL
#writeLines("\n\n####################################################################\n####################################################################\n\n(1-3)Levels Hierfstat Script\n\n(\u0254) Copyleft 2012, by LDGH and Contributors. GNU GPL 2012\n\n####################################################################\n\nPs: Notice that two Enters are required after you type the required informations!\n\nEnter the directory path where the Hierfstat input files are:")
#files.dir = scan(file="", what="character")
######## ERROR LOOKUP - MORE THAN ONE PATH SPECIFIED
#if (length(files.dir) > 1) { stop("I'm crashing! Did you insertesd more than one path?")}
######## SHOW ENTERED PATH
#writeLines("You specified this file path:")
#print(files.dir)
######## SHOW FILES IN SPECIFIED PATH
#writeLines("With those files:")
#files.dir = c("/home/giordano/Hierfstat4/")
#all.the.files = list.files(path = files.dir,full.names=TRUE)
#all.out.files = list.files(path = files.dir,full.names=FALSE)
all.the.files = args[1]
all.out.files = paste0("Results_",args[1])
######## ERROR LOOKUP - NO FILES IN THIS PATH OR INEXISTENT PATH
#if (length(all.the.files) < 1) { stop("I'm crashing! Did you made a mistake? Verify the specified path and the files that are supposed to exist there!")}
######## PRINT FILE NAMES
#for (row_data in 1: length(all.the.files)) {
#print(all.the.files[row_data], quote=FALSE)
#}
######## ENTER THE NUMBER OF HIERARCHICAL LEVELS (1-3)
#writeLines("Enter the desired number of hierarchical levels (1-3):")
#hier_lvls = as.numeric(scan(file="", what="character"))
hier_lvls = 2
######## ERROR LOOKUP - NOT ACCEPTED VALUES
#if (length(hier_lvls) > 1) { stop("I'm crashing! Did you entered more than one value?")}
#if (hier_lvls > 3) { stop("I'm crashing! Did you entered a value greater than 3?")}
#if (hier_lvls < 1) { stop("I'm crashing! Did you entered a value lesser than 1?")}
######## PERFORMS AMOVA FOR EACH FILE


for (row_data in 1: dim(as.matrix(all.the.files))[1])
{                                                    
	writeLines("Performing analysis for this file:")
	print(all.the.files[row_data], quote=FALSE)
	print(Sys.time())
	data1 = fread(as.matrix(all.the.files)[row_data,1],sep = "auto", header=TRUE, stringsAsFactors=FALSE)
	writeLines("Upload completed")
	print(Sys.time())
	levels1 = as.numeric(unlist(data1[[2]]))
	data1 = data1[,3:ncol(data1),with=FALSE]

	writeLines("Starting computation")
	print(Sys.time())
              
    ########## AMOVA FOR ONE HIERARCHICAL LEVEL
	hier = 1	
	if (hier == 1) {
	colClasses=c(POLYM="character",FST="character",FIT="character",FIS="character")
	fstats_table = dataFrame(colClasses,nrow=0) 
	tableFst = foreach(i=1:length(data1), iter=icount(),.errorhandling="pass", .combine="rbind") %dopar% {
	fst_values1 <- try (varcomp(data.frame(cbind(as.numeric(unlist(levels1))),as.numeric(unlist(data1[[iter]])))), TRUE)
			if (is.numeric(try (fst_values1$F, TRUE))) 
			{
				fstats_table[1,1] = try (colnames(data1)[iter], TRUE)
				fstats_table[1,2] = fst_values1$F[1,1]
				fstats_table[1,3] = fst_values1$F[1,2]
				fstats_table[1,4] = fst_values1$F[2,2]
			}	else
			{
				fstats_table[1,1] = try (colnames(data1)[iter], TRUE)
				fstats_table[1,2] = NA
				fstats_table[1,3] = NA
				fstats_table[1,4] = NA
			}
#		print(fstats_table)			
		fstats_table[is.na(fstats_table)]<- 0
		return(fstats_table)
		}
    	}
	writeLines("Computation done")
	print(Sys.time())
	write.table(tableFst, file=all.out.files, append = F, col.names=TRUE, row.names=FALSE, quote=FALSE)
	writeLines("Writing completed")
	print(Sys.time())

}
