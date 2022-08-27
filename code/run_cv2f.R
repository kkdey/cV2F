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
  make_option("--numPCs", default=250, type="character", help="The number of annotation PCs to use"),
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
# opt$feature_file = "/n/groups/price/kushal/ENCODE/data/Deliverables/Feature_Tables/Allfeatures_cCRE_V2F_perchr."
# opt$mafpath = "/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/MAF_features_Aug032022.txt"
# opt$bimpath = "/n/groups/price/kushal/extras/BIMS_hg38/1000G.EUR.QC."
# opt$ldblockspath = "/n/groups/price/kushal/LAVA/loci_files/LAVA_LDblocks_published.txt"
# opt$numPCs = as.integer(opt$numPCs)

pre_positive_snps = read.table(opt$positive_set, header=F)[,1]
pre_negative_snps = read.table(opt$negative_set, header=F)[,1]
maf_tabb = data.frame(fread(opt$mafpath))
LDblocks_tabb = data.frame(fread(opt$ldblockspath))

##########################################################################################################

cat("Define positive and LD and MAF matched negative set of variants")
cat("\n")

bimpooltabb = c()
for(numchr in 1:22){
  bimdf = data.frame(fread(paste0(opt$bimpath, numchr, ".bim")))
  bimpooltabb = rbind(bimpooltabb, bimdf)
  cat("Processing SNPs from chr:", numchr, "\n")
}

colnames(bimpooltabb) = c("CHR", "SNP", "CM", "BP", "A1", "A2")

positive_snps = intersect(pre_positive_snps, maf_tabb$SNP)

pre_negative_snps = setdiff(intersect(pre_negative_snps, maf_tabb$SNP), pre_positive_snps)
negative_snps = c()
for(nn in 1:length(positive_snps)){
  bp = bimpooltabb$BP[which(bimpooltabb$SNP == positive_snps[nn])]
  chr = bimpooltabb$CHR[which(bimpooltabb$SNP == positive_snps[nn])]
  maf= maf_tabb$MAF[which(bimpooltabb$SNP == positive_snps[nn])]
  idx1 = which(LDblocks_tabb$START < bp  & LDblocks_tabb$STOP > bp  &
                 LDblocks_tabb$CHR == chr)[1]
  idx2 = which(bimpooltabb$CHR == LDblocks_tabb$CHR[idx1] &
                 bimpooltabb$BP > LDblocks_tabb$START[idx1] &
                 bimpooltabb$BP < LDblocks_tabb$STOP[idx1])

  if(length(idx2) < 3){
    idx2 = which(bimpooltabb$BP < bp + 50000 & bimpooltabb$BP > bp - 50000 &
                   bimpooltabb$CHR == chr)
  }
  potential_negatives = intersect(bimpooltabb$SNP[idx2], pre_negative_snps)
  mafs_potential_negatives = maf_tabb$MAF[match(potential_negatives, maf_tabb$SNP)]
  snp_temp = potential_negatives[order(abs(mafs_potential_negatives - maf), decreasing = F)[1:2]]
  negative_snps = c(negative_snps, snp_temp)
  cat("Matching for SNP:", nn, "\n")
}

negative_snps = negative_snps[!is.na(negative_snps)]


positive_snps2 = c()
negative_snps2 = c()
pos_feature_tabb_list = c()
neg_feature_tabb_list = c()

for(numchr in 1:22){
  dff = data.frame(fread(paste0(opt$feature_file, numchr, ".txt")))
  pos_feature_tabb_list[[numchr]] = dff[match(intersect(positive_snps, dff$SNP), dff$SNP), -(1:4)]
  neg_feature_tabb_list[[numchr]] = dff[match(intersect(negative_snps, dff$SNP), dff$SNP), -(1:4)]
  positive_snps2 = c(positive_snps2, intersect(positive_snps, dff$SNP))
  negative_snps2 = c(negative_snps2, intersect(negative_snps, dff$SNP))
  cat("We are at chr:", numchr, "\n")
}

pos_feature_tabb = do.call(rbind, pos_feature_tabb_list)
neg_feature_tabb = do.call(rbind, neg_feature_tabb_list)
rownames(pos_feature_tabb) = positive_snps2
rownames(neg_feature_tabb) = negative_snps2




combined_feature_tabb = rbind(pos_feature_tabb, neg_feature_tabb)
rownames(combined_feature_tabb) = c(positive_snps2, negative_snps2)
labels = c(rep(1, nrow(pos_feature_tabb)), rep(0, nrow(neg_feature_tabb)))

##########################################################################################################

if(is.null(opt$numPCs)){
  annotation_feature_tabb = combined_feature_tabb
}else{
  prr = prcomp(combined_feature_tabb, tol=0.01)
  annotation_feature_tabb = prr$x[, 1:min(opt$numPCs, ncol(prr$x))]
}
rownames(annotation_feature_tabb) = c(positive_snps2, negative_snps2)

cat("Compute performance metrics to assess classofication accuracy")
cat("\n")

aucvec = c()
prcvec = c()
for(num_iter in 1:6){

  training_idx = sample(1:length(labels), floor(0.5*length(labels)), replace=F)
  test_idx = setdiff(1:length(labels), training_idx)

  training_labels = labels[training_idx]
  test_labels = labels[test_idx]

  training_data = annotation_feature_tabb[training_idx, ]
  test_data = annotation_feature_tabb[test_idx, ]

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

cat("Create cV2F")
cat("\n")

cv2f_list = c()
for(numchr in 1:22){
  bimtabb = data.frame(fread(paste0(opt$bimpath, numchr, ".bim")))
  temp_positive_snps = setdiff(positive_snps, bimtabb[,2])
  temp_negative_snps = setdiff(negative_snps, bimtabb[,2])

  temp_pos_feature_tabb = annotation_feature_tabb[match(temp_positive_snps, rownames(annotation_feature_tabb)), ]
  temp_neg_feature_tabb = annotation_feature_tabb[match(temp_negative_snps, rownames(annotation_feature_tabb)), ]
  temp_combined_feature_tabb = rbind(temp_pos_feature_tabb, temp_neg_feature_tabb)
  temp_labels = c(rep(1, nrow(temp_pos_feature_tabb)), rep(0, nrow(temp_neg_feature_tabb)))

  training_data = temp_combined_feature_tabb
  training_labels = temp_labels
  dff = data.frame(fread(paste0(opt$feature_file, numchr, ".txt")))
  test_data = dff[match(bimtabb[,2], dff$SNP), -(1:4)]
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
