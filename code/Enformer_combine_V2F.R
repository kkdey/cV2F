
ccre_features = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/cCRE/cCRE_v04_13Jun2022.txt"))

enformer_rsdir = "/n/groups/price/kushal/Imperio/ALL_DEEP/Enformer/Celltypes_P_sigs"
ll2 = as.character(sapply(list.files(enformer_rsdir), function(x) return(strsplit(x, ".txt")[[1]][1])))
annotpool_list = list()

for(numl in 1:length(ll2)){
  rsids = unique(read.table(paste0(enformer_rsdir, "/", ll2[numl],".txt"), header=F)[,1])
  annot = rep(0, nrow(ccre_features))
  annot[match(intersect(ccre_features$SNP, rsids), ccre_features$SNP)] = 1
  annotpool_list[[numl]] = annot
  cat("We are at biosample:", numl, "\n")
}

annotpool = do.call(cbind, annotpool_list)
colnames(annotpool) = ll2
annotpool2 = cbind.data.frame(ccre_features[,1:4], annotpool)
fwrite(annotpool2, file = "/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/Enformer_allelic_effect_features_Aug032022.txt",
       row.names=F, col.names=T, sep = "\t", quote=F)
