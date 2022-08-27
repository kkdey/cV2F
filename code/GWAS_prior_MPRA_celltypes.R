
tissue = "K562"
mpra_match = "K562"
traits_match = "RBC"
ccre_match = "K562"


tissue = "GM12878"
mpra_match = "GM12878"
traits_match = c("WBC", "Lym")
ccre_match = "GM12878"

tissue = "HepG2"
mpra_match = "HEPG2"
traits_match = c("TBil")
ccre_match = "HepG2"

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
disease_idx = match(traits_match, diseases)
diseases = diseases[disease_idx]
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
ccre_features = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/cCRE/cCRE_celltypes_v04_13Jun2022.txt"))
ccre_snps = ccre_features$SNP[which(ccre_features[paste0(ccre_match, "_PLS")]==1 |
                                    ccre_features[paste0(ccre_match, "_pELS")] == 1 |
                                    ccre_features[paste0(ccre_match, "_dELS")] == 1)]
varscore[which(disease_tabb[,1] %in% ccre_snps == T)]=2

###########################################   Evidence 3 (cV2F)  #######################################################################

library(data.table)
cv2f_features = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/cV2F/cV2F_4Dec2021.txt"))
cv2f_snps = cv2f_features$SNP[which(cv2f_features$cV2F.bin == 1)]
varscore[which(disease_tabb[,1] %in% cv2f_snps == T)]=3


###########################################   Evidence 4 (MPRA)  #######################################################################

library(data.table)
mpra_gwas_features = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/MPRA_GWAS_allelic_effect_features_Aug032022.txt"))
mpra_snps = mpra_gwas_features$SNP[mpra_gwas_features[paste0(mpra_match)] == 1]
varscore[which(disease_tabb[,1] %in% mpra_snps == T)]=4

##########################################  Final   #####################################################

disease_tabb$VarScore = varscore
colnames(disease_tabb) = c("SNP", "CHR", "BP", "Trait", "VarScore")

write.table(disease_tabb, file = paste0("/n/groups/price/kushal/ENCODE/data/Deliverables/VariantPro/variant_scores_traits_PIP90_", tissue, ".txt"),
            col.names = T, row.names = F, sep = "\t", quote=F)


length(which(disease_tabb$VarScore == 4))


xx1 = read.table(paste0("/n/groups/price/kushal/ENCODE/data/Deliverables/VariantPro/variant_scores_traits_PIP90", ".txt"), header=T)
xx2 = read.table(paste0("/n/groups/price/kushal/ENCODE/data/Deliverables/VariantPro/variant_scores_traits_PIP90_K562", ".txt"), header=T)
xx3 = read.table(paste0("/n/groups/price/kushal/ENCODE/data/Deliverables/VariantPro/variant_scores_traits_PIP90_GM12878", ".txt"), header=T)
xx4 = read.table(paste0("/n/groups/price/kushal/ENCODE/data/Deliverables/VariantPro/variant_scores_traits_PIP90_HepG2", ".txt"), header=T)

oo1 = unique(xx1$SNP[xx1$VarScore == 4])

oo2 = unique(c(unique(xx2$SNP[xx2$VarScore == 4]),
unique(xx3$SNP[xx3$VarScore == 4]),
unique(xx4$SNP[xx4$VarScore == 4])))

