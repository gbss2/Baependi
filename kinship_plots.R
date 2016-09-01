#==========================================
#
# (/u0254) Copyleft 2016, by GBSS and Contributors.
#
# 
# -----------------
# kinship.r
# -----------------
# GNU GPL 2016, by GBSS and Contributors.
#
# Original Author: Giordano Bruno Soares-Souza
# Contributor(s): 
# Updated by:
#
# Command line: Rscript kinship(<plink.genome,REAP_pairs_relatedness file>)
# Dependencies: R, library gtools
# Description: Compare IBS and kinship between different datasets
# Future developments:
# 	1) Generalize script
#
################################################################################

plinkEAS = read.table("/scratch/capeverd/shared/Baependi/new/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_pruned_qc_me_EUR_EAS_AFR_noDrift.genome.tab",head=T,as.is=T)
plinkNAT = read.table("/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.genome.tab",head=T,as.is=T)
reapEAS3 = read.table("/scratch/capeverd/gbss2/baependi/admixture/K3_EAS/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_pruned_qc_me_EUR_EAS_AFR_noDrift_REAP_pairs_relatedness_3.txt",head=T,as.is=T)
reapNAT3 = read.table("/scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_REAP_pairs_relatedness_3.txt",head=T,as.is=T)
reapEAS5 = read.table("/scratch/capeverd/gbss2/baependi/admixture/K5_EAS/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_pruned_qc_me_EUR_EAS_AFR_noDrift_REAP_pairs_relatedness_5.txt",head=T,as.is=T)
reapNAT5 = read.table("/scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_REAP_pairs_relatedness_5.txt",head=T,as.is=T)

colnames(plinkEAS)=c("FID1","IID1","FID2","IID2","RTe","EZe","Z0e","Z1e","Z2e","PI_HATe","PHEe","DSTe","PPCe","RATIOe")
colnames(plinkNAT)=c("FID1","IID1","FID2","IID2","RTn","EZn","Z0n","Z1n","Z2n","PI_HATn","PHEn","DSTn","PPCn","RATIOn")

colnames(reapEAS3)=c("FID1","IID1","FID2","IID2","N_SNPe3","IBD0_PROBe3","IBD1_PROBe3","IBD2_PROBe3","KINCOEFe3")
colnames(reapEAS5)=c("FID1","IID1","FID2","IID2","N_SNPe5","IBD0_PROBe5","IBD1_PROBe5","IBD2_PROBe5","KINCOEFe5")
colnames(reapNAT3)=c("FID1","IID1","FID2","IID2","N_SNPn3","IBD0_PROBn3","IBD1_PROBn3","IBD2_PROBn3","KINCOEFn3")
colnames(reapNAT5)=c("FID1","IID1","FID2","IID2","N_SNPn5","IBD0_PROBn5","IBD1_PROBn5","IBD2_PROBn5","KINCOEFn5")

t1 = merge(plinkNAT,plinkEAS)
t2 = merge(t1,reapNAT3)
t3 = merge(t2,reapNAT5)
t4 = merge(t3,reapEAS3)
kinship = merge(t4,reapEAS5)

rm(t1)
rm(t2)
rm(t3)
rm(t4)

png("kinship_plots%02d.png")
qqplot(kinship$PI_HATn/2,kinship$PI_HATe/2,xlab = "PI_HAT/2 Dataset NAT", ylab = "PI_HAT/2 Dataset EAS")
qqplot(kinship$PI_HATn/2,kinship$KINCOEFn3,xlab = "PI_HAT/2 Natives", ylab = "KINCOEF Natives K3")
qqplot(kinship$PI_HATn/2,kinship$KINCOEFn5,xlab = "PI_HAT/2 Natives", ylab = "KINCOEF Natives K5")
qqplot(kinship$PI_HATe/2,kinship$KINCOEFe3,xlab = "PI_HAT/2 EAS", ylab = "KINCOEF EAS K3")
qqplot(kinship$PI_HATe/2,kinship$KINCOEFe5,xlab = "PI_HAT/2 EAS", ylab = "KINCOEF EAS K5")
qqplot(kinship$Z0n,kinship$Z0e,xlab = "IBD 0 Dataset NAT (Plink)", ylab = "IBD 0 Dataset EAS (Plink)")
qqplot(kinship$Z0n,kinship$IBD0_PROBn3,xlab = "IBD 0 Natives (Plink)", ylab = "IBD 0 NAT K3 (REAP)")
qqplot(kinship$Z0n,kinship$IBD0_PROBn5,xlab = "IBD 0 Natives (Plink)", ylab = "IBD 0 NAT K5 (REAP)")
qqplot(kinship$Z0e,kinship$IBD0_PROBe3,xlab = "IBD 0 EAS (REAP)", ylab = "IBD 0 EAS K5 (REAP)")
qqplot(kinship$Z0e,kinship$IBD0_PROBe5,xlab = "IBD 0 EAS (REAP)", ylab = "IBD 0 EAS K5 (REAP)")

qqplot(kinship$PI_HATn/2,kinship$KINCOEFn3,xlab = "PI_HAT/2", ylab = "KINCOEF",pch=1)
points(sort(kinship$PI_HATn/2), sort(kinship$KINCOEFn5), col = "red",pch=1)
points(sort(kinship$PI_HATn/2), sort(kinship$KINCOEFe3), col = "green",pch=1)
points(sort(kinship$PI_HATn/2), sort(kinship$KINCOEFe5), col = "blue",pch=1)
legend("bottomright", legend = c("NAT K3", "NAT K5", "EAS K3", "EAS K5"), pch = 1, col = c("black", "red","green","blue"))

qqplot(kinship$Z0n,kinship$IBD0_PROBn3,xlab = "IBD 0 Plink", ylab = "IBD 0 REAP",pch=1)
points(sort(kinship$Z0n), sort(kinship$IBD0_PROBn5), col = "red",pch=1)
points(sort(kinship$Z0n), sort(kinship$IBD0_PROBe3), col = "green",pch=1)
points(sort(kinship$Z0n), sort(kinship$IBD0_PROBe5), col = "blue",pch=1)
legend("bottomright", legend = c("NAT K3", "NAT K5", "EAS K3", "EAS K5"), pch = 1, col = c("black", "red","green","blue"))

dev.off()




