library(data.table)
library(R.utils)

option_list <- list(
  make_option("--vierstra_dir", type="character", default = "data/vierstra", help="Directory where you have downloaded the Vierstra allelic imbalance files"),
  make_option("--samp_v2f_file", type="character", default = "data/vierstra", help="Name of a sample V2F file for co-ordinates"),
  make_option("--out_v2f_file", type="character", default = "data/vierstra", help="Output of Vierstra V2F file")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)


vierstra_dir = opt$vierstra_dir

#vierstra_dir="/n/groups/price/kushal/ENCODE/data/Vierstra/Version2"
#samp_v2f_file="/n/groups/price/kushal/ENCODE/data/Deliverables/cCRE/cCRE_v04_13Jun2022.txt"
#out_v2f_file="/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/DNaseAI_Vierstra_allelic_effect_features_Aug032022.txt"

ccre_features = data.frame(fread(paste0(samp_v2f_file)))

dff = data.frame(fread(paste0(vierstra_dir, "/", "all.hets.combined.pvals.txt.gz")))
rsids_test = dff$dbsnp[grep("rs", dff$dbsnp)]

biosamples = read.table(paste0(vierstra_dir, "/", "sample_names2.txt"), header=F)[,1]
biosamples_out = c()
all_rsids = c()
annotpool_list = list()
kk=1
for(numl in 1:length(biosamples)){
  if(!file.exists(paste0(vierstra_dir, "/", "Celltypes_P_sigs", "/", "rsids_FDR10_", biosamples[numl], ".txt"))){
    next
  }
  if(file.info(paste0(vierstra_dir, "/", "Celltypes_P_sigs", "/", "rsids_FDR10_", biosamples[numl], ".txt"))$size > 0){
    rsids = read.delim(paste0(vierstra_dir, "/", "Celltypes_P_sigs", "/", "rsids_FDR10_", biosamples[numl], ".txt"))[,1]
    all_rsids = c(all_rsids, rsids)
    if(length(rsids) > 1000){
      annot = rep(0, nrow(ccre_features))
      annot[match(intersect(ccre_features$SNP, rsids_test), ccre_features$SNP)] = length(rsids_test)/nrow(ccre_features)
      annot[match(intersect(ccre_features$SNP, rsids), ccre_features$SNP)] = 1
      biosamples_out = c(biosamples_out, biosamples[numl])
      annotpool_list[[kk]] = annot
      kk=kk+1
    }
  }
  cat("We are at biosample:", numl, "\n")
}

annotpool = do.call(cbind, annotpool_list)
colnames(annotpool) = biosamples_out

all_rsids = unique(all_rsids)

annot = rep(0, nrow(ccre_features))
annot[match(intersect(ccre_features$SNP, rsids_test), ccre_features$SNP)] = length(rsids_test)/nrow(ccre_features)
annot[match(intersect(ccre_features$SNP, all_rsids), ccre_features$SNP)] = 1
annotpool$all = annot

annotpool2 = cbind.data.frame(ccre_features[,1:4], annotpool)
fwrite(annotpool2, file = paste0(out_v2f_file), row.names=F, col.names=T, sep = "\t", quote=F)



