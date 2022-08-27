
library(data.table)
library(R.utils)
library(optparse)
library(xgboost)
library(pROC)
library(PRROC)


options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

option_list <- list(
  make_option("--positive_set", type="character", default = "data/positive_set.txt",
              help="rsIDs of SNPs in the positive set (celltype relevant trait)"),
  make_option("--negative_set", type="character", default = "data/negative_set.txt",
              help="rsIDs of SNPs in the negative set (non-celltype relevant traits)"),
  make_option("--feature_tabb", type="character", default = "data/feature_table.txt",
              help="Data frame of feature tables related to a cell type"),
  make_option("--output_name", type="character", default = "data/out",
              help="Name of the output file")

)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

positive_set = opt$positive_set
negative_set = opt$negative_set
feature_tabb = opt$feature_tabb
output_name = opt$output_name

positive_set = "/n/groups/price/kushal/ENCODE/data/GWAS_traits/Finemap_SNPsets/TBil_PIP0_1.txt"
negative_set = "/n/groups/price/kushal/ENCODE/data/GWAS_traits/Finemap_SNPsets/Notfinemapped_PIP0_001.txt"
feature_file = "/n/groups/price/kushal/ENCODE/data/Deliverables/Feature_Tables/K562_features_Aug032022.txt"

baseline_tabb = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/BaselineLDv2.2_features_Aug032022.txt"))
maf_tabb = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/MAF_features_Aug032022.txt"))

atabb = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/Allele_effects_all_features_v04_6Nov2021.txt"))

pre_positive_snps = read.table(positive_set, header=F)[,1]
pre_negative_snps = read.table(negative_set, header=F)[,1]

feature_tabb = data.frame(fread(feature_file))
pre_negative_snps = setdiff(intersect(pre_negative_snps, feature_tabb$SNP), pre_positive_snps)

positive_snps = intersect(pre_positive_snps, feature_tabb$SNP)
negative_snps = c()
for(nn in 1:length(pre_positive_snps)){
  bp = feature_tabb$BP[which(feature_tabb$SNP == pre_positive_snps[nn])]
  chr = feature_tabb$CHR[which(feature_tabb$SNP == pre_positive_snps[nn])]
  maf= maf_tabb$MAF[which(feature_tabb$SNP == pre_positive_snps[nn])]
  idx1 = which(feature_tabb$BP > bp - 50000 & feature_tabb$BP < bp - 1000 &
              feature_tabb$CHR == chr)
  idx2 = which(feature_tabb$BP < bp + 50000 & feature_tabb$BP > bp + 1000 &
              feature_tabb$CHR == chr)
  idx = c(idx1, idx2)
  potential_negatives = intersect(feature_tabb$SNP[idx], pre_negative_snps)
  mafs_potential_negatives = maf_tabb$MAF[match(potential_negatives, maf_tabb$SNP)]
  snp_temp = potential_negatives[order(abs(mafs_potential_negatives - maf), decreasing = F)[1]]
  negative_snps = c(negative_snps, snp_temp)
}
negative_snps = negative_snps[!is.na(negative_snps)]


pos_feature_tabb = feature_tabb[match(positive_snps, feature_tabb$SNP), -(1:4)]
neg_feature_tabb = feature_tabb[match(negative_snps, feature_tabb$SNP), -(1:4)]
combined_feature_tabb = rbind(pos_feature_tabb, neg_feature_tabb)
labels = c(rep(1, nrow(pos_feature_tabb)), rep(0, nrow(neg_feature_tabb)))


aucvec = c()
prcvec = c()
for(num_iter in 1:10){

  training_idx = sample(1:length(labels), floor(0.5*length(labels)), replace=F)
  test_idx = setdiff(1:length(labels), training_idx)

  training_labels = labels[training_idx]
  test_labels = labels[test_idx]

  training_data = combined_feature_tabb[training_idx, ]
  test_data = combined_feature_tabb[test_idx, ]

  NUMITER_XGBoost=5000
  bstSparse <-  xgboost(data = as.matrix(training_data),
                        label = training_labels,
                        nrounds = NUMITER_XGBoost,
                        objective = "binary:logistic",
                        eval_metric = "auc")
  predict_labels = predict(bstSparse, as.matrix(test_data))
  roc_test <- roc(test_labels, predict_labels, algorithm = 2)
  probs1 = predict_labels[test_labels == 1]
  probs2 = predict_labels[test_labels == 0]
  pr <- pr.curve(scores.class0 = probs1, scores.class1 = probs2, curve = T)
  aucvec = c(aucvec, auc(roc_test))
  prcvec = c(prcvec, pr$auc.integral)
  cat("We are at iter:", num_iter, "\n")
}

dff = data.frame("AUROC.mean" = mean(aucvec),
                 "AUROC.sd" = sd(aucvec),
                 "AUPRC.mean" = mean(prcvec),
                 "AUPRC.sd" = sd(prcvec))

save(dff, file = paste0(opt$output_name))




