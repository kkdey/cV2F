

library(data.table)
annotpool_list = list()
for(numchr in 1:22){
  dff = data.frame(fread(paste0("/n/groups/price/kushal/ENCODE/data/ANNOTATIONS_hg38/Baselines/baselineLD_v2.2/",
                                "baselineLD.", numchr, ".annot.gz")))
  annotpool_list[[numchr]] = dff
  cat("We are at chr:", numchr, "\n")
}

annotpool = do.call(rbind, annotpool_list)

fwrite(annotpool, file = "/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/BaselineLDv2.2_features_Aug032022.txt",
       row.names=F, col.names=T, sep = "\t", quote=F)


library(data.table)
annotpool_list = list()
for(numchr in 1:22){
  dff = data.frame(fread(paste0("/n/groups/price/kushal/LDSC/1000G_EUR_Phase3_plink_hg38/",
                                "1000G.EUR.QC.", numchr, ".frq")))
  annotpool_list[[numchr]] = dff
  cat("We are at chr:", numchr, "\n")
}

annotpool = do.call(rbind, annotpool_list)

fwrite(annotpool, file = "/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/MAF_features_Aug032022.txt",
       row.names=F, col.names=T, sep = "\t", quote=F)
