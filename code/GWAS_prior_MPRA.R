

###########################################   Set of all SNPs  #######################################################################

bimfiles = list.files("/n/groups/price/kushal/extras/BIMS_hg38")
bimpooltabb = c()
for(numchr in 1:22){
  bimtabb = read.table(paste0("/n/groups/price/kushal/extras/BIMS_hg38", "/", "1000G.EUR.QC.", numchr, ".bim"))
  bimpooltabb = rbind(bimpooltabb, bimtabb[, c(1, 4, 2)])
  cat("We are at chr:", numchr, "\n")
}
colnames(bimpooltabb) = c("CHR", "BP", "SNP")

###########################################   Evidence 1 (GWAS finemapping)  #######################################################################

diseases = list.files("/n/groups/price/kushal/ENCODE/data/GWAS_traits/Jamboree_GWAS_Traits/191010_UKBB_SUSIE/CredibleSets/")
diseases = diseases[c(55, 92, 72, 82)]
SNPset = list()
disease_tabb = c()
for(numd in 1:length(diseases)){
  vv = read.table(paste0("/n/groups/price/kushal/ENCODE/data/GWAS_traits/Jamboree_GWAS_Traits/191010_UKBB_SUSIE/CredibleSets/",
                         diseases[numd], "/", "variant.list.txt"), header = T)
  rsids = intersect(vv$variant[vv$pip > 0.90], bimpooltabb$SNP)
  SNPset[[numd]] = rsids
  disease_tabb = rbind(disease_tabb, cbind(rsids,
                       paste0("chr", bimpooltabb$CHR[match(rsids, bimpooltabb$SNP)]),
                       bimpooltabb$BP[match(rsids, bimpooltabb$SNP)],
                       diseases[numd]))
  cat("We are processing disease:", numd, "\n")
}
names(SNPset) = diseases

disease_tabb = as.data.frame(disease_tabb)

varscore = rep(1, nrow(disease_tabb))

###########################################   Evidence 2 (cCRE)  #######################################################################

library(data.table)
ccre_features = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/cCRE/cCRE_v04_13Jun2022.txt"))
ccre_snps = ccre_features$SNP[which(ccre_features$PLS==1 | ccre_features$pELS == 1 | ccre_features$dELS == 1)]
varscore[which(disease_tabb[,1] %in% ccre_snps == T)]=2


###########################################   Evidence 3 (cV2F)  #######################################################################

library(data.table)
cv2f_features = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/cV2F/cV2F_4Dec2021.txt"))
cv2f_snps = cv2f_features$SNP[which(cv2f_features$cV2F.bin == 1)]
varscore[which(disease_tabb[,1] %in% cv2f_snps == T)]=3


###########################################   Evidence 4 (MPRA)  #######################################################################

library(data.table)
mpra_gwas_features = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/MPRA_GWAS_allelic_effect_features_Aug032022.txt"))
mpra_snps = mpra_gwas_features$SNP[mpra_gwas_features$all == 1]
varscore[which(disease_tabb[,1] %in% mpra_snps == T)]=4

##########################################  Final   #####################################################

disease_tabb$VarScore = varscore
colnames(disease_tabb) = c("SNP", "CHR", "BP", "Trait", "VarScore")

write.table(disease_tabb, file = "/n/groups/price/kushal/ENCODE/data/Deliverables/VariantPro/variant_scores_traits_PIP90.txt",
            col.names = T, row.names = F, sep = "\t", quote=F)






disease_tabb2 = disease_tabb[which(disease_tabb$Trait != "Lym"), ]




tt = disease_tabb2[which(disease_tabb2$VarScore==4), ]
tt2 = tt[which(cv2f_features$cV2F.ALLV2F[match(tt$SNP, cv2f_features$SNP)] > 0.98), ]

length(which(cv2f_features$cV2F.ALLV2F[match(tt$SNP, cv2f_features$SNP)] > 0.99))

