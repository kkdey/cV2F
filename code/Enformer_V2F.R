library(data.table)
library(R.utils)
library(optparse)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

option_list <- list(
  make_option("--celltype_mode", type="character", default = "data/celltypes_mode.txt", help="Celltype and mode of mark")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

celltype_mode = opt$celltype_mode #all_all


bimfiles = list.files("/n/groups/price/kushal/extras/BIMS_hg38")
bimpooltabb = c()
for(numchr in 1:22){
  bimtabb = read.table(paste0("/n/groups/price/kushal/extras/BIMS_hg38", "/", "1000G.EUR.QC.", numchr, ".bim"))
  bimpooltabb = rbind(bimpooltabb, bimtabb[, c(1, 4, 2)])
  cat("We are at chr:", numchr, "\n")
}
colnames(bimpooltabb) = c("CHR", "BP", "SNP")


ll = list.files("/n/groups/price/kushal/Imperio/ALL_DEEP/Enformer/Z_SKEW_FEATURES")


maxz2_all = list()
meanz2_all = list()
rsids = c()
for(numchr in 1:22){
  tabb = data.frame(fread(paste0("/n/groups/price/kushal/Imperio/ALL_DEEP/Enformer/Z_SKEW_FEATURES", "/",
                           "chr", numchr, ".", celltype_mode, ".txt")))
  tabb = tabb[grep("rs", tabb$SNP), ]
  rsids = c(rsids, tabb$SNP)
  maxz2_all[[numchr]] = tabb$MAXZ
  meanz2_all[[numchr]] = tabb$MEANZ
  cat("We are at chr:", numchr, "\n")
}
maxz2_all_scores = unlist(maxz2_all)
meanz2_all_scores = unlist(meanz2_all)
common_snps = intersect(bimpooltabb$SNP, rsids)
maxz2_all_scores = maxz2_all_scores[match(common_snps, rsids)]
meanz2_all_scores = meanz2_all_scores[match(common_snps, rsids)]

cutoffs = c(floor(0.001*length(common_snps)),
            floor(0.01*length(common_snps)),
            floor(0.1*length(common_snps)))
rsids_max1 = common_snps[order(maxz2_all_scores, decreasing = T)[1:cutoffs[1]]]
rsids_max2 = common_snps[order(maxz2_all_scores, decreasing = T)[1:cutoffs[2]]]
rsids_max3 = common_snps[order(maxz2_all_scores, decreasing = T)[1:cutoffs[3]]]

rsids_mean1 = common_snps[order(meanz2_all_scores, decreasing = T)[1:cutoffs[1]]]
rsids_mean2 = common_snps[order(meanz2_all_scores, decreasing = T)[1:cutoffs[2]]]
rsids_mean3 = common_snps[order(meanz2_all_scores, decreasing = T)[1:cutoffs[3]]]

write.table(rsids_max1, file = paste0("/n/groups/price/kushal/Imperio/ALL_DEEP/Enformer/Celltypes_P_sigs/Enformer_Max_",
            celltype_mode, "_0.1percent.txt"), row.names = F, col.names = F, sep = "\t", quote=F)
write.table(rsids_max2, file = paste0("/n/groups/price/kushal/Imperio/ALL_DEEP/Enformer/Celltypes_P_sigs/Enformer_Max_",
                                      celltype_mode, "_1percent.txt"), row.names = F, col.names = F, sep = "\t", quote=F)
write.table(rsids_max3, file = paste0("/n/groups/price/kushal/Imperio/ALL_DEEP/Enformer/Celltypes_P_sigs/Enformer_Max_",
                                      celltype_mode, "_10percent.txt"), row.names = F, col.names = F, sep = "\t", quote=F)

write.table(rsids_mean1, file = paste0("/n/groups/price/kushal/Imperio/ALL_DEEP/Enformer/Celltypes_P_sigs/Enformer_Mean_",
                                      celltype_mode, "_0.1percent.txt"), row.names = F, col.names = F, sep = "\t", quote=F)
write.table(rsids_mean2, file = paste0("/n/groups/price/kushal/Imperio/ALL_DEEP/Enformer/Celltypes_P_sigs/Enformer_Mean_",
                                      celltype_mode, "_1percent.txt"), row.names = F, col.names = F, sep = "\t", quote=F)
write.table(rsids_mean3, file = paste0("/n/groups/price/kushal/Imperio/ALL_DEEP/Enformer/Celltypes_P_sigs/Enformer_Mean_",
                                      celltype_mode, "_10percent.txt"), row.names = F, col.names = F, sep = "\t", quote=F)


# c("all", )
