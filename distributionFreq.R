
  #==========================================
  # R file
  #==========================================
  #
  # (/u0254) Copyleft 2016
  #
  # -----------------
  # distributionFreq.R
  # -----------------
  # GNU GPL 2016
  #
  # Original Author: 
  # Contributor(s): 
  # Updated by:
  #
  # Command line: 
  # Dependencies: R
  # Description: 
  # 
  #############################################

# 1) Upload Fam files

    baepFampopNAT3 = read.table('K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.fampop',head=F,as.is=T)
    baepFampopNAT5 = read.table('K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.fampop',head=F,as.is=T)
    baepFampopEAS3 = read.table('K3_EAS/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_pruned_qc_me_EUR_EAS_AFR_noDrift.fampop',head=F,as.is=T)
    baepFampopEAS5 = read.table('K5_EAS/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_pruned_qc_me_EUR_EAS_AFR_noDrift.fampop',head=F,as.is=T)

# 2) Upload Admixture Freq files

    propNAT3 = read.table('K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.3.Q',head=F,as.is=T)
    propNAT5 = read.table('K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.5.Q',head=F,as.is=T)
    propEAS3 = read.table('K3_EAS/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_pruned_qc_me_EUR_EAS_AFR_noDrift.3.Q',head=F,as.is=T)
    propEAS5 = read.table('K5_EAS/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_pruned_qc_me_EUR_EAS_AFR_noDrift.5.Q',head=F,as.is=T)

# 3) Data Merge

    colnames(baepFampopNAT5) = c("FID","IID","PID","MID","SEX","PHENO","POP_N5","CON_N5")
    colnames(baepFampopNAT3) = c("FID","IID","PID","MID","SEX","PHENO","POP_N3","CON_N3")
    colnames(baepFampopEAS5) = c("FID","IID","PID","MID","SEX","PHENO","POP_E5","CON_E5")
    colnames(baepFampopEAS3) = c("FID","IID","PID","MID","SEX","PHENO","POP_E3","CON_E3")
    
    colnames(propNAT3) = c("EUR_N3","AFR_N3","NAT_N3")
    colnames(propNAT5) = c("NEUR_N5","SEUR_N5","WAFR_N5","NAT_N5","EAFR_N5")
    colnames(propEAS3) = c("EUR_E3","EAS_E3","AFR_E3")
    colnames(propEAS5) = c("NEUR_E5","EAS_E5","SEUR_E5","WAFR_E5","EAFR_E5")
    
    nat3 = cbind(baepFampopNAT3[1:1838,],propNAT3[1:1838,])
    nat5 = cbind(baepFampopNAT5[1:1838,],propNAT5[1:1838,])
    eas3 = cbind(baepFampopEAS3[1:1838,],propEAS3[1:1838,])
    eas5 = cbind(baepFampopEAS5[1:1838,],propEAS5[1:1838,])
    
    admix = cbind(nat3,nat5[,7:13],eas3[,7:11],eas5[,7:13])

# 4) Plot Frequencies distribution

# Dataset 1 - NAT

    pdf('ancestry_qc_NAT.pdf')
    par(mfrow=c(2,4))
    hist(admix$EUR_N3,xlim=c(0,1),xlab="",main="NAT K3 - EUR")
    hist(admix$NAT_N3,xlim=c(0,1),xlab="",main="NAT K3 - NAT")
    hist(admix$AFR_N3,xlim=c(0,1),xlab="",main="NAT K3 - AFR")
    hist(admix$NEUR_N5,xlim=c(0,1),xlab="",main="NAT K5 - NEUR")
    hist(admix$NAT_N5,xlim=c(0,1),xlab="",main="NAT K5 - NAT")
    hist(admix$SEUR_N5,xlim=c(0,1),xlab="",main="NAT K5 - SEUR")
    hist(admix$WAFR_N5,xlim=c(0,1),xlab="",main="NAT K5 - WAFR")
    hist(admix$EAFR_N5,xlim=c(0,1),xlab="",main="NAT K5 - EAFR")
    mtext("Distribution of Ancestries - NAT dataset", outer = TRUE, cex = 1.5)
    dev.off()

# Dataset 2 - EAS

    pdf('ancestry_qc_EAS.pdf')
    par(mfrow=c(2,4))
    hist(admix$EUR_E3,xlim=c(0,1),xlab="",main="EAS K3 - EUR")
    hist(admix$EAS_E3,xlim=c(0,1),xlab="",main="EAS K3 - EAS")
    hist(admix$AFR_E3,xlim=c(0,1),xlab="",main="EAS K3 - AFR")
    hist(admix$NEUR_E5,xlim=c(0,1),xlab="",main="EAS K5 - NEUR")
    hist(admix$EAS_E5,xlim=c(0,1),xlab="",main="EAS K5 - EAS")
    hist(admix$SEUR_E5,xlim=c(0,1),xlab="",main="EAS K5 - SEUR")
    hist(admix$WAFR_E5,xlim=c(0,1),xlab="",main="EAS K5 - WAFR")
    hist(admix$EAFR_E5,xlim=c(0,1),xlab="",main="EAS K5 - EAFR")
    mtext("Distribution of Ancestries - EAS dataset", outer = TRUE, cex = 1.5)
    dev.off()

# 5) Plot Frequencies correlation

# Dataset 1 - NAT

pdf('ancestry_corr_NAT.pdf',width=12,height=9)
pdf('ancestry_corr_NAT.pdf')
par(mfrow=c(3,2))
plot(admix$EUR_N3,(admix$NAT_N3+admix$AFR_N3),xlab="NAT K3 - EUR",ylab="NAT K3 - NAT+AFR")
mtext(paste0("rho = ",cor(admix$EUR_N3,(admix$NAT_N3+admix$AFR_N3))),cex= 0.7, side=3, adj=0.9, line=-1.5)
plot((admix$SEUR_N5+admix$NEUR_N5),(admix$NAT_N5+admix$WAFR_N5+admix$EAFR_N5),xlab="NAT K5 - EUR (NEUR+SEUR)",ylab="NAT K5 - NAT+AFR (WAFR+EAFR)")
mtext(paste0("rho = ",cor((admix$SEUR_N5+admix$NEUR_N5),(admix$NAT_N5+admix$WAFR_N5+admix$EAFR_N5))),cex= 0.7, side=3, adj=0.9, line=-1.5)
plot(admix$EUR_N3,(admix$NEUR_N5+admix$SEUR_N5),xlab="NAT K3 - EUR",ylab="NAT K5 - NEUR+SEUR")
mtext(paste0("rho = ",cor(admix$EUR_N3,(admix$NEUR_N5+admix$SEUR_N5))),cex= 0.7, side=1, adj=0.9, line=-1.5)
plot(admix$EUR_N3,(admix$NEUR_N5+admix$SEUR_N5),xlab="NAT K3 - EUR",ylab="NAT K5 - NEUR+SEUR",main="EUR > 0.6",ylim=c(0.6,1),xlim=c(0.6,1))
mtext(paste0("rho = ",cor(admix[(admix$NEUR_N5+admix$SEUR_N5)>0.9,]$EUR_N3,(admix[(admix$NEUR_N5+admix$SEUR_N5)>0.9,]$NEUR_N5+admix[(admix$NEUR_N5+admix$SEUR_N5)>0.9,]$SEUR_N5))),cex= 0.7, side=1, adj=0.9, line=-1.5)
plot(admix$AFR_N3,(admix$WAFR_N5+admix$EAFR_N5),xlab="NAT K3 - AFR",ylab="NAT K5 - WAFR+EAFR")
mtext(paste0("rho = ",cor(admix$AFR_N3,(admix$WAFR_N5+admix$EAFR_N5))),cex= 0.7, side=1, adj=0.9, line=-1.5)
plot(admix$AFR_N3,(admix$WAFR_N5+admix$EAFR_N5),xlab="NAT K3 - AFR",ylab="NAT K5 - WAFR+EAFR",main="AFR < 0.1",ylim=c(0,0.1),xlim=c(0,0.1))
mtext(paste0("rho = ",cor(admix[(admix$WAFR_N5+admix$EAFR_N5)<0.1,]$AFR_N3,(admix[(admix$WAFR_N5+admix$EAFR_N5)<0.1,]$WAFR_N5+admix[(admix$WAFR_N5+admix$EAFR_N5)<0.1,]$EAFR_N5))),cex= 0.7, side=1, adj=0.9, line=-1.5)
dev.off()

# Dataset 2 - EAS

pdf('ancestry_corr_EAS.pdf')
par(mfrow=c(3,2))
plot(admix$EUR_E3,(admix$EAS_E3+admix$AFR_E3),xlab="EAS K3 - EUR",ylab="EAS K3 - EAS+AFR")
mtext(paste0("rho = ",cor(admix$EUR_E3,(admix$EAS_E3+admix$AFR_E3))),cex= 0.7, side=3, adj=0.9, line=-1.5)
plot((admix$SEUR_E5+admix$NEUR_E5),(admix$EAS_E5+admix$WAFR_E5+admix$EAFR_E5),xlab="EAS K5 - EUR (NEUR+SEUR)",ylab="EAS K5 - EAS+AFR (WAFR+EAFR)")
mtext(paste0("rho = ",cor((admix$SEUR_E5+admix$NEUR_E5),(admix$EAS_E5+admix$WAFR_E5+admix$EAFR_E5))),cex= 0.7, side=3, adj=0.9, line=-1.5)
plot(admix$EUR_E3,(admix$NEUR_E5+admix$SEUR_E5),xlab="EAS K3 - EUR",ylab="EAS K5 - NEUR+SEUR")
mtext(paste0("rho = ",cor(admix$EUR_E3,(admix$NEUR_E5+admix$SEUR_E5))),cex= 0.7, side=1, adj=0.9, line=-1.5)
plot(admix$EUR_E3,(admix$NEUR_E5+admix$SEUR_E5),xlab="EAS K3 - EUR",ylab="EAS K5 - NEUR+SEUR",main="EUR > 0.6",ylim=c(0.6,1),xlim=c(0.6,1))
mtext(paste0("rho = ",cor(admix[(admix$NEUR_E5+admix$SEUR_E5)>0.9,]$EUR_E3,(admix[(admix$NEUR_E5+admix$SEUR_E5)>0.9,]$NEUR_E5+admix[(admix$NEUR_E5+admix$SEUR_E5)>0.9,]$SEUR_E5))),cex= 0.7, side=1, adj=0.9, line=-1.5)
plot(admix$AFR_E3,(admix$WAFR_E5+admix$EAFR_E5),xlab="EAS K3 - AFR",ylab="EAS K5 - WAFR+EAFR")
mtext(paste0("rho = ",cor(admix$AFR_E3,(admix$WAFR_E5+admix$EAFR_E5))),cex= 0.7, side=1, adj=0.9, line=-1.5)
plot(admix$AFR_E3,(admix$WAFR_E5+admix$EAFR_E5),xlab="EAS K3 - AFR",ylab="EAS K5 - WAFR+EAFR",main="AFR < 0.1",ylim=c(0,0.1),xlim=c(0,0.1))
mtext(paste0("rho = ",cor(admix[(admix$WAFR_E5+admix$EAFR_E5)<0.1,]$AFR_E3,(admix[(admix$WAFR_E5+admix$EAFR_E5)<0.1,]$WAFR_E5+admix[(admix$WAFR_E5+admix$EAFR_E5)<0.1,]$EAFR_E5))),cex= 0.7, side=1, adj=0.9, line=-1.5)
dev.off()
