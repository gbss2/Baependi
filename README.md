# Baependi Analysis Pipeline

<b>Baependi releases:</b>

* Release 1 = 1117 ind
* Release 2 = 430 ind
* Release 3 = 721 ind
  
<b>Baependi Release 1+3:</b>
* 1838 individuals
*     polymorphisms

<b>Datasets for population genomics analysis:</b>

<i>a) Dataset 1 - Baependi + 1000Genomes + HGDP</i>
  * 2906 individuals
  * 62437 polymorphisms
  * Pops included: <b>Baependi</b>; <b>EUR</b>ope (Bergamo, CEU, French, GBR, IBS, TSI, Tuscan); <b>AFR</b>ican (ESN, GWD, LWK, Mandenka, YRI); <b>NAT</b>ive-american (Karitiana, Maya, Piapoco and Curripaco, Pima, Surui) 
  
<i>b) Dataset 2 - Baependi + 1000Genomes</i>
  *  individuals
  *  polymorphisms
  * Pops included: <b>Baependi</b>; <b>EUR</b>ope (Bergamo, CEU, French, GBR, IBS, TSI, Tuscan); <b>AFR</b>ican (ESN, GWD, LWK, Mandenka, YRI); <b>EAS</b>t Asian (CHB, CHD)
  
<b>Dataset 1 (Baependi + 1000 Genomes + HGDP) preparation</b>

<b>Note</b>: SNPs were pruned and selected from those included in HGDP dataset. Previous attempts resulted in less markers than this approach (eg. the merge between dataset 2 (Baependi and 1000 Genomes resulted in 19000 markers).

1) Select SNPs to extract from complete datasets


    cat Baependi_release3_filter_1_merged_chr*.rs > Baependi_release3_filter_1_merged_chrs.rs

    cat Baependi_release3_filter_05_merged_chr*.rs > Baependi_release3_filter_05_merged_chrs.rs

    cut -f2 HGDP_971Report_Forward_strand.bim > HGDP_971Report_Forward_strand.rs

    comm -1 -2 <(sort Baependi_release3_filter_1_merged_chrs.rs) <(sort HGDP_971Report_Forward_strand.rs) > intersect_baependi_hgdp.txt

    comm -1 -2 <(sort Baependi_release3_filter_05_merged_chrs.rs) <(sort HGDP_971Report_Forward_strand.rs) > intersect_05_baependi_hgdp.txt

    comm -1 -2  <(sort Baependi_release3_all_merged_chrs.rs )  <(sort HGDP_971Report_Forward_strand.rs) > intersect_all_baependi_hgdp.txt


2) Extract SNPs from Baependi imputed dataset (based in the intersection between HGDP and Baependi)

    for chr in $(seq 1 22)
    do

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --gen /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3_chr${chr}.impute2.gz --sample /alice/scratch/capeverd/shared/Baependi/Baependi_release3_batch1_3.sample --biallelic-only strict --oxford-pheno-name plink_pheno --oxford-single-chr ${chr} --keep-allele-order --extract /alice/scratch/capeverd/shared/Baependi/intersect_all_baependi_hgdp.txt --freq --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3_chr${chr}.impute2.hgdp_rs

    awk '{if($5 == 0 || $5 == 1) print $2 }' /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3_chr${chr}.impute2.hgdp_rs.frq > /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3_chr${chr}.impute2.snps2remove

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3_chr${chr}.impute2.hgdp_rs --keep-allele-order --exclude /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3_chr${chr}.impute2.snps2remove --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3_chr${chr}.impute2_hgdp_rs

    done

    gzip /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3_chr${chr}.impute2.hgdp_rs.*

3) Merge Bapendi Chrs files

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3_chr1.impute2_hgdp_rs --merge-list /alice/scratch/capeverd/shared/Baependi/baependi_merge2.txt --keep-allele-order --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs.temp 

4) Update fam file (sex and parents)

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs.temp --keep-allele-order --update-parents /home/g/gbss2/baependi_update_parents.txt --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs.temp2

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs.temp2 --keep-allele-order --update-sex /home/g/gbss2/baependi_update_sex.txt --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs

    rm /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs.temp2.*

    rm /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs.temp.* 


5) Baependi statistics and pruning

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs --hardy --keep-allele-order --freq --missing --mendel --het --write-snplist --indep-pairwise 50 10 0.1 --out /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs

6) Checking Mendelian errors (phase 1)

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs --set-me-missing --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs_me 

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs_me --hardy --keep-allele-order --freq --missing --mendel --het --write-snplist  --out /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs_me

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs_me --geno 0.05 --mind 0.1 --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs_me_qc

    cat Baependi_release3_filtered_batch1_3.impute2_hgdp_rs.lmendel |  sed -r 's/^\s+//g' | sed -r 's/\s+/\t/g' | awk '{ if($3!="N" && $3>0) print $0}'  | sort -k3 -n  | tail

    cat Baependi_release3_filtered_batch1_3.impute2_hgdp_rs.imendel |  sed -r 's/^\s+//g' | sed -r 's/\s+/\t/g' | awk '{ if($3!="N" && $3>0) print $0}'  | sort -k3 -n  | tail

    cat Baependi_release3_filtered_batch1_3.impute2_hgdp_rs.fmendel |  sed -r 's/^\s+//g' | sed -r 's/\s+/\t/g' | awk '{ if($5!="N" && $3>0) print $0}'  | sort -k5 -n  | tail

7) 1000 Genomes dataset

    gunzip /scratch/capeverd/gbss2/temp/gbss2/Documents/1000Genomes/FullData/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz

    vcftools --vcf /scratch/capeverd/gbss2/temp/gbss2/Documents/1000Genomes/FullData/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf --temp /scratch/capeverd/tmp/  --plink-tped --out ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --tfile /alice/scratch/capeverd/shared/Baependi/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_hgdp_rs --biallelic-only strict --keep-allele-order --make-bed --out /alice/scratch/capeverd/shared/Baependi/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_hgdp_rs

    grep -v '^#' /scratch/capeverd/gbss2/temp/gbss2/Documents/1000Genomes/FullData/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf | cut -f 3 | sort | uniq -d > /alice/scratch/capeverd/shared/Baependi/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.dups

    /scratch/capeverd/gbss2/soft/plink_1.9/plink/plink --bfile /scratch/capeverd/gbss2/temp/gbss2/Documents/1000Genomes/FullData/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf \
    --biallelic-only strict --keep-allele-order \
    --exclude /alice/scratch/capeverd/shared/Baependi/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.dups \
    --make-bed --out /alice/scratch/capeverd/shared/Baependi/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes

8) Merge Baependi and 1000 Genomes dataset (first step)

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/ALL.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes_hgdp_rs --bmerge /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs_pruned --extract /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs.prune.in --keep-allele-order --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs_pruned

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs --extract /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs.prune.in --keep-allele-order --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs_pruned

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs_pruned --flip /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs-merge.missnp  --keep-allele-order --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs_strand

    (# flip was not necessary) 

9) Merge Baependi + 1000 Genomes and HGDP

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs_pruned --bmerge /alice/scratch/capeverd/shared/Baependi/HGDP_971Report_Forward_strand --extract /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs.prune.in --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs_HGDP_pruned

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/HGDP_971Report_Forward_strand --flip /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs_HGDP_pruned-merge.missnp  --keep-allele-order --make-bed --out /alice/scratch/capeverd/shared/Baependi/HGDP_971Report_Forward_strand_v2

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs_pruned --bmerge /alice/scratch/capeverd/shared/Baependi/HGDP_971Report_Forward_strand_v2 --extract /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs.prune.in --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs_HGDP_pruned

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs_pruned --exclude /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs_HGDP_pruned-merge.missnp --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs_pruned_temp

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/HGDP_971Report_Forward_strand_v2  --exclude /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs_HGDP_pruned-merge.missnp --make-bed --out /alice/scratch/capeverd/shared/Baependi/HGDP_971Report_Forward_strand_v2_temp

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs_pruned_temp --bmerge /alice/scratch/capeverd/shared/Baependi/HGDP_971Report_Forward_strand_v2_temp --extract /alice/scratch/capeverd/shared/Baependi/Baependi_release3_filtered_batch1_3.impute2_hgdp_rs.prune.in --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs_HGDP_pruned

10) Quality control  (geno = 0.02 and mind = 0.1)

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs_HGDP_pruned --keep-allele-order --geno 0.02 --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs_HGDP_pruned_qc1 

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_hgdp_rs_HGDP_pruned_qc1 --keep-allele-order --mind 0.1 --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc

11) Check AT/CG flips

    grep -c 'C\sG\|G\sC\|A\sT\|T\sA' /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me.bim > Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me.ATGC.txt

    (No AT or CG SNPs in the dataset)

12) Final dataset statistics 

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc --hardy --keep-allele-order --freq --missing --mendel --het --write-snplist --out /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc

13) Set Mendelian errors as missing

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc --set-me-missing --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me 

14) Filter by continents

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me --exclude /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me.ATGC.txt  --within /alice/scratch/capeverd/shared/Baependi/ids_all_pops_con_hgdp_1000G_cv_baep_nat_CON.clust --keep-cluster-names EAFR WAFR EUR NAT ADMX --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR

15) Filter by populations and maf = 0.01

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR --within /alice/scratch/capeverd/shared/Baependi/ids_all_pops_con_hgdp_1000G_cv_baep_nat_POP.clust --maf 0.01 --remove-cluster-names ACB ASW CLM MXL PEL PUR MKK Bantu FrenchBasque Adygei Russian San BiakaPygmy MbutiPygmy Sardinian FIN Orcadian --make-bed --out /alice/scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift

<i>Final dataset = 62437 variants and 2906 people pass filters and QC<\i>

16) Admixture input (K = 3)

    Rscript /alice/scratch/capeverd/shared/Baependi/fampop.R /alice/scratch/capeverd/shared/Baependi/ids_all_pops_con_hgdp_1000G_cv_baep_nat.txt /scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.fam

    sed -i 's/ADMX/-/g' /scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.pop

    cp /scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.* /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/

17) Admixture input (K = 5)

    Rscript /alice/scratch/capeverd/shared/Baependi/fampop.R /alice/scratch/capeverd/shared/Baependi/ids_all_pops_con_hgdp_1000G_cv_baep_nat_v2.txt /scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.fam

    sed -i 's/ADMX/-/g' /scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.pop

    cp /scratch/capeverd/shared/Baependi/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.* /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/

18) Admixture Run K3

    /scratch/capeverd/gbss2/soft/admixture_linux-1.23/admixture -j6 --supervised --cv /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.bed 3 > /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.3.log

    Rscript /home/g/gbss2/admixPlot.r 3 Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.3.Q

19) Admixture Run K5

    /scratch/capeverd/gbss2/soft/admixture_linux-1.23/admixture -j6 --supervised --cv /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.bed 5 > /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.5.log

    Rscript /home/g/gbss2/admixPlot.r 5 Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.5.Q

20) REAP input K3

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift --recode12 --output-missing-genotype 0 --transpose --out /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift

    cut -d' ' -f1,2 /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.fam > /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_12cols.txt

    paste /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_12cols.txt  /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.3.Q > /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.proportions

21) REAP input K5

    /scratch/capeverd/gbss2/soft/plink_1.9/plink --bfile /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift --recode12 --output-missing-genotype 0 --transpose --out /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift

    cut -d' ' -f1,2 /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.fam > /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_12cols.txt

    paste /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_12cols.txt  /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.5.Q > /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.proportions

22) Run REAP K3

    /scratch/capeverd/gbss2/soft/REAP/./REAP -g /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.tped -p /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.tfam -a /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.proportions -f /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.3.P -r 1 -k 3 -t -0.1 > /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_reap.log

    mv /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/REAP_Individual_Index.txt /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_REAP_Individual_Index_3.txt
    mv /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/REAP_Inbreed.txt /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_REAP_Inbreed_3.txt
    mv /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/REAP_Kincoef_matrix.txt /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_REAP_Kincoef_matrix_3.txt
    mv /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/REAP_IBD0_matrix.txt /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_REAP_IBD0_matrix_3.txt
    mv /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/REAP_pairs_relatedness.txt /scratch/capeverd/gbss2/baependi/admixture/K3_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_REAP_pairs_relatedness_3.txt

    awk '{if($9 > 0.4) print $0}' Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_REAP_pairs_relatedness_3.txt

23) Run REAP K5

    /scratch/capeverd/gbss2/soft/REAP/./REAP -g /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.tped -p /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.tfam -a /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.proportions -f /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift.5.P -r 1 -k 5 -t -0.1 > /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_reap.log

    mv /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/REAP_Individual_Index.txt /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_REAP_Individual_Index_5.txt
    mv /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/REAP_Inbreed.txt /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_REAP_Inbreed_5.txt
    mv /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/REAP_Kincoef_matrix.txt /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_REAP_Kincoef_matrix_5.txt
    mv /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/REAP_IBD0_matrix.txt /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_REAP_IBD0_matrix_5.txt
    mv /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/REAP_pairs_relatedness.txt /scratch/capeverd/gbss2/baependi/admixture/K5_NAT/Baependi_r3_batch1_3_impute2_ALL_phase3_v5a.20130502_HGDP_pruned_qc_me_EUR_NAT_AFR_noDrift_REAP_pairs_relatedness_5.txt
