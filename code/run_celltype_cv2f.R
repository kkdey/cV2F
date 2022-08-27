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
  make_option("--feature_file", type="character", default = "data/feature_table.txt",
              help="Data frame of feature tables related to a cell type"),
  make_option("--bimpath", type="character", help="Path and prefix of the bimfile of the BIMFILE"),
  make_option("--mafpath", type="character", help="Path and prefix of the frequency file for all SNPs"),
  make_option("--ldblockspath", type="character", help="Path and prefix of the LD blocks file"),
  make_option("--output_cv2f", type="character", default = "data/out",
              help=" Path of the output celltype cV2F file"),
  make_option("--output_metrics", type="character", default = "data/out",
              help=" Path of the output AUROC/AUPRC metrics file")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

# opt = list()
# opt$positive_set =  "/n/groups/price/kushal/ENCODE/data/GWAS_traits/Finemap_SNPsets/TBil_PIP0_1.txt"
# opt$negative_set =  "/n/groups/price/kushal/ENCODE/data/GWAS_traits/Finemap_SNPsets/Notfinemapped_PIP0_001.txt"
# opt$feature_file = "/n/groups/price/kushal/ENCODE/data/Deliverables/Feature_Tables/HepG2_features_Aug032022.txt"
# opt$mafpath = "/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/MAF_features_Aug032022.txt"
# opt$bimpath = "/n/groups/price/kushal/extras/BIMS/1000G.EUR.QC."
# opt$ldblockspath = "/n/groups/price/kushal/LAVA/loci_files/LAVA_LDblocks_published.txt"


pre_positive_snps = read.table(opt$positive_set, header=F)[,1]
pre_negative_snps = read.table(opt$negative_set, header=F)[,1]
feature_tabb = data.frame(fread(opt$feature_file))
maf_tabb = data.frame(fread(opt$mafpath))
LDblocks_tabb = data.frame(fread(opt$ldblockspath))

##########################################################################################################

cat("Define positive and LD and MAF matched negative set of variants")
cat("\n")

positive_snps = intersect(pre_positive_snps, feature_tabb$SNP)

pre_negative_snps = setdiff(intersect(pre_negative_snps, feature_tabb$SNP), pre_positive_snps)
negative_snps = c()
for(nn in 1:length(pre_positive_snps)){
  bp = feature_tabb$BP[which(feature_tabb$SNP == pre_positive_snps[nn])]
  chr = feature_tabb$CHR[which(feature_tabb$SNP == pre_positive_snps[nn])]
  maf= maf_tabb$MAF[which(feature_tabb$SNP == pre_positive_snps[nn])]
  idx1 = which(LDblocks_tabb$START < bp  & LDblocks_tabb$STOP > bp  &
                LDblocks_tabb$CHR == chr)[1]
  idx2 = which(feature_tabb$CHR == LDblocks_tabb$CHR[idx1] &
    feature_tabb$BP > LDblocks_tabb$START[idx1] &
    feature_tabb$BP < LDblocks_tabb$STOP[idx1])

  if(length(idx2) < 3){
    idx2 = which(feature_tabb$BP < bp + 100000 & feature_tabb$BP > bp - 100000 &
                   feature_tabb$CHR == chr)
  }
  potential_negatives = intersect(feature_tabb$SNP[idx2], pre_negative_snps)
  mafs_potential_negatives = maf_tabb$MAF[match(potential_negatives, maf_tabb$SNP)]
  snp_temp = potential_negatives[order(abs(mafs_potential_negatives - maf), decreasing = F)[1:2]]
  negative_snps = c(negative_snps, snp_temp)
}
negative_snps = negative_snps[!is.na(negative_snps)]

##########################################################################################################

cat("Build feature data and labels matrix")
cat("\n")

pos_feature_tabb = feature_tabb[match(positive_snps, feature_tabb$SNP), -(1:4)]
neg_feature_tabb = feature_tabb[match(negative_snps, feature_tabb$SNP), -(1:4)]
combined_feature_tabb = rbind(pos_feature_tabb, neg_feature_tabb)
labels = c(rep(1, nrow(pos_feature_tabb)), rep(0, nrow(neg_feature_tabb)))


##########################################################################################################

cat("Compute performance metrics to assess classofication accuracy")
cat("\n")

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
save(dff, file = paste0(opt$output_metrics))


##########################################################################################################

cat("Create cell-type level cV2F")
cat("\n")

cv2f_list = c()
for(numchr in 1:22){
  bimtabb = data.frame(fread(paste0(opt$bimpath, numchr, ".bim")))
  temp_positive_snps = setdiff(positive_snps, bimtabb[,2])
  temp_negative_snps = setdiff(negative_snps, bimtabb[,2])

  temp_pos_feature_tabb = feature_tabb[match(temp_positive_snps, feature_tabb$SNP), -(1:4)]
  temp_neg_feature_tabb = feature_tabb[match(temp_negative_snps, feature_tabb$SNP), -(1:4)]
  temp_combined_feature_tabb = rbind(temp_pos_feature_tabb, temp_neg_feature_tabb)
  temp_labels = c(rep(1, nrow(temp_pos_feature_tabb)), rep(0, nrow(temp_neg_feature_tabb)))

  training_data = temp_combined_feature_tabb
  training_labels = temp_labels
  test_data = feature_tabb[match(bimtabb[,2], feature_tabb$SNP), -(1:4)]
  NUMITER_XGBoost=5000
  bstSparse <-  xgboost(data = as.matrix(training_data),
                        label = training_labels,
                        nrounds = NUMITER_XGBoost,
                        objective = "binary:logistic",
                        eval_metric = "auc")
  predict_labels = predict(bstSparse, as.matrix(test_data))
  cv2f_list[[numchr]] = cbind.data.frame(bimtabb[,1], bimtabb[,4], bimtabb[,2], bimtabb[,3], predict_labels)
  cat("We are at chr:", numchr, "\n")
}

cv2f_df = do.call(rbind, cv2f_list)
colnames(cv2f_df) = c("CHR", "BP", "SNP", "CM", "cV2F")
fwrite(cv2f_df, file = paste0(opt$output_cv2f), row.names=F, col.names=T, sep = "\t", quote=F)



