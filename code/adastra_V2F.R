
library(data.table)
library(R.utils)

## Download link for ADASTRA Bill Cipher: https://adastra.autosome.org/bill-cipher/downloads
## Local directory:
## opt$adastra_dir = "/Users/kushaldey/Documents/ENCODE4AWG/data/ADASTRA/ADASTRA_BillCipher"
## Cluster directory:
## adastra_dir = "/n/groups/price/kushal/ENCODE/data/ADASTRA_BillCipher"
## opt$annot_dir = "/n/groups/price/kushal/ENCODE/data/ANNOTATIONS_hg38/ADASTRA_BillCipher"
## opt$base_dir = "/n/groups/price/kushal/ENCODE/data/ANNOTATIONS_hg38/Baselines/"

option_list <- list(
  make_option("--adastra_dir", type="character", default = "data/adastra", help="Directory where you have downloaded the ADASTRA files"),
  make_option("--annot_dir", type="character", default = "data/annotations/adastra", help="Directory where to save the ADASTRA V2F annotations"),
  make_option("--baseline_dir", type="character", default = "data/annotations/baselines", help="Directory with baseline annotations in S-LDSC")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)


#####################################  Extract all tested variants in ADASTRA   ######################################################################

adastra_dir = opt$adastra_dir
annot_dir = opt$annot_dir
base_dir = opt$base_dir


tf_files = list.files(paste0(adastra_dir, "/", "adastra.bill_cipher/release_dump", "/TF/"))
tf_names = as.character(sapply(tf_files, function(x) return(strsplit(x, "_HUMAN")[[1]][1])))
library(data.table)
merged_tabb_list = list()
for(numm in 1:length(tf_names)){
  tabb = data.frame(fread(paste0(adastra_dir, "/", "adastra.bill_cipher/release_dump", "/TF/", tf_files[numm])))
  idx = 1:nrow(tabb)
  if(length(idx) > 0){
    tabb2 = tabb[idx, ]
    tabb2$TF = tf_names[numm]
    merged_tabb_list[[numm]] =  tabb2
  }
  cat("We are at TF:", tf_files[numm], "\n")
}

merged_tabb = do.call(rbind, merged_tabb_list)
rsids = unique(merged_tabb$ID)
write.table(rsids, file = paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_tested_adastra.txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)

#####################################  Extract cell line level ASB variants in ADASTRA   ######################################################################

cell_lines = c("HEK293", "HeLa", "GM1", "K562", "A549", "HCT-116", "LNCaP", "MCF7", "SK-N-SH", "HepG2")
cl_files2 = list.files(paste0(adastra_dir, "/", "adastra.bill_cipher/release_dump", "/CL/"))
cl_names = as.character(sapply(cl_files2, function(x) return(strsplit(x, "_.tsv")[[1]][1])))

for(cc in cell_lines){
  cl_files = cl_files2[grep(paste0(cc), cl_names)]
  merged_tabb = c()
  for(numm in 1:length(cl_files)){
    tabb = data.frame(fread(paste0(adastra_dir, "/", "adastra.bill_cipher/release_dump", "/CL/", "/",
                                   cl_files[numm])))
    idx = unique(c(which(tabb$fdrp_bh_ref < 0.10), which(tabb$fdrp_bh_alt < 0.10)))
    if(length(idx) > 0){
      tabb2 = tabb[idx, ]
      tabb2$CL = cl_files[numm]
      merged_tabb = rbind(merged_tabb, tabb2)
    }
    cat("We are at CL:", cl_files[numm], "\n")
  }

  rsids = unique(merged_tabb$ID)
  write.table(rsids, file =  paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_", cc, ".txt"),
              row.names = F, col.names = F, sep = "\t", quote=F)
}


#####################################  Extract tissue level ASB variants in ADASTRA   ######################################################################

cl_files = list.files(paste0(adastra_dir, "/", "adastra.bill_cipher/release_dump", "/CL/"))
cl_names = as.character(sapply(cl_files, function(x) return(strsplit(x, "_.tsv")[[1]][1])))
library(data.table)
merged_tabb = c()
for(numm in 1:length(cl_names)){
  tabb = data.frame(fread(paste0(adastra_dir, "/", "adastra.bill_cipher/release_dump", "/CL/", "/",
                                 cl_files[numm])))
  idx = unique(c(which(tabb$fdrp_bh_ref < 0.10), which(tabb$fdrp_bh_alt < 0.10)))
  if(length(idx) > 0){
    tabb2 = tabb[idx, ]
    tabb2$CL = cl_names[numm]
    merged_tabb = rbind(merged_tabb, tabb2)
  }
  cat("We are at CL:", cl_files[numm], "\n")
}

write.table(merged_tabb, file = paste0(adastra_dir, "/", "AllCL_ADASTRA_sig_SNPs_FDR10_refalt.txt"),
            row.names = F, col.names = T, sep = "\t", quote=F)

snp_tabb1 = data.frame(fread(paste0(adastra_dir, "/", "AllCL_ADASTRA_sig_SNPs_FDR10_refalt.txt")))
rsids_all = unique(snp_tabb1$ID)
write.table(rsids_all, file = paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_", "all", ".txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)

snp_tabb1 = data.frame(fread(paste0(adastra_dir, "/", "AllCL_ADASTRA_sig_SNPs_FDR10_refalt.txt")))
ll = unique(snp_tabb1$CL)
fix_celltypes = c("CD", "GM12", "T_", "B_", "monocy", "neutro", "Mono", "Neutro", "Blood", "blood", "immune", "_B-", "_T-",
                  "-B-", "-T-", "lympho", "Lympho")
ll1 = unique(ll[unlist(sapply(fix_celltypes, function(x) return(grep(x, ll))))])
snp_tabb2 = snp_tabb1[ which(snp_tabb1$CL %in% ll1 == T), ]
rsids_all = unique(snp_tabb2$ID)
write.table(rsids_all, file = paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_", "blood", ".txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)

ll = unique(snp_tabb1$CL)
fix_celltypes = c("Liver", "liver", "HepG2", "Hepato", "hepato")
ll1 = unique(ll[unlist(sapply(fix_celltypes, function(x) return(grep(x, ll))))])
snp_tabb2 = snp_tabb1[ which(snp_tabb1$CL %in% ll1 == T), ]
rsids_all = unique(snp_tabb2$ID)
write.table(rsids_all, file = paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_", "liver", ".txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)

ll = unique(snp_tabb1$CL)
fix_celltypes = c("Lung", "lung", "pulmo", "Pulmo")
ll1 = unique(ll[unlist(sapply(fix_celltypes, function(x) return(grep(x, ll))))])
snp_tabb2 = snp_tabb1[ which(snp_tabb1$CL %in% ll1 == T), ]
rsids_all = unique(snp_tabb2$ID)
write.table(rsids_all, file = paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_", "lung", ".txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)

ll = unique(snp_tabb1$CL)
fix_celltypes = c("Skin", "skin", "kerat", "Kerat")
ll1 = unique(ll[unlist(sapply(fix_celltypes, function(x) return(grep(x, ll))))])
snp_tabb2 = snp_tabb1[ which(snp_tabb1$CL %in% ll1 == T), ]
rsids_all = unique(snp_tabb2$ID)
write.table(rsids_all, file = paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_", "skin", ".txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)

ll = unique(snp_tabb1$CL)
fix_celltypes =c("Intestine", "intestine", "bowel", "Bowel", "colon", "Colon", "Stomach", "stomach", "gut")
ll1 = unique(ll[unlist(sapply(fix_celltypes, function(x) return(grep(x, ll))))])
snp_tabb2 = snp_tabb1[ which(snp_tabb1$CL %in% ll1 == T), ]
rsids_all = unique(snp_tabb2$ID)
write.table(rsids_all, file = paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_", "gut", ".txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)


ll = unique(snp_tabb1$CL)
fix_celltypes =c("Kidney", "kidney", "renal", "Renal")
ll1 = unique(ll[unlist(sapply(fix_celltypes, function(x) return(grep(x, ll))))])
ll1 = ll1[-11]
snp_tabb2 = snp_tabb1[ which(snp_tabb1$CL %in% ll1 == T), ]
rsids_all = unique(snp_tabb2$ID)
write.table(rsids_all, file = paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_", "kidney", ".txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)


ll = unique(snp_tabb1$CL)
fix_celltypes = c("Brain", "brain", "neur", "Neur", "Spin", "spin", "_astro", "Astro")
ll1 = unique(ll[unlist(sapply(fix_celltypes, function(x) return(grep(x, ll))))])
snp_tabb2 = snp_tabb1[ which(snp_tabb1$CL %in% ll1 == T), ]
rsids_all = unique(snp_tabb2$ID)
write.table(rsids_all, file = paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_", "brain", ".txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)


ll = unique(snp_tabb1$CL)
fix_celltypes = c("heart", "Heart")
ll1 = unique(ll[unlist(sapply(fix_celltypes, function(x) return(grep(x, ll))))])
snp_tabb2 = snp_tabb1[ which(snp_tabb1$CL %in% ll1 == T), ]
rsids_all = unique(snp_tabb2$ID) ## very few picked
write.table(rsids_all, file = paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_", "heart", ".txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)


ll = unique(snp_tabb1$CL)
fix_celltypes = c("adipo")
ll1 = unique(ll[unlist(sapply(fix_celltypes, function(x) return(grep(x, ll))))])
snp_tabb2 = snp_tabb1[ which(snp_tabb1$CL %in% ll1 == T), ]
rsids_all = unique(snp_tabb2$ID) ## very few picked
write.table(rsids_all, file = paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_", "fat", ".txt"),
            row.names = F, col.names = F, sep = "\t", quote=F)



####################   Generate V2F annotations for ADASTRA corresponding to all biosamples  #######################################


bio_types = as.character(sapply(as.character(sapply(list.files(paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/"), pattern = "rsids_FDR"),
                                                    function(x) return(strsplit(x, "rsids_")[[1]][2]))), function(x) return(strsplit(x, ".txt")[[1]][1])))

rsids_test = read.table(paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_tested_adastra.txt"), header=F)[,1]

for(numchr in 1:22){
  base = data.frame(fread(paste0( base_dir, "/",
                                 "baseline_Epi_hg38", "/",
                                 "baselineLD.", numchr, ".annot.gz")))
  for(tt in 1:length(bio_types)){
    rsids_tt = read.table(paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/",
                                 "rsids_", bio_types[tt], ".txt"), header=F)[,1]
    length(rsids_tt)
    annot1 = rep(0, nrow(base))
    annot1[which(base$SNP %in% rsids_test == T)]=1

    annot2 = rep(0, nrow(base))
    annot2[which(base$SNP %in% rsids_tt == T)]=1

    newdf = cbind.data.frame(base[,1:4], annot2, annot1)
    colnames(newdf) = c(colnames(base)[1:4], "ANNOT", "Tested")

    if(!dir.exists(paste0(annot_dir, "/", "ADASTRA_", bio_types[tt]))){
      dir.create(paste0(annot_dir, "/", "ADASTRA_", bio_types[tt]))
    }
    write.table(newdf, file = gzfile(paste0(annot_dir, "/",
                                            "ADASTRA_", bio_types[tt], "/",
                                            "ADASTRA_", bio_types[tt], ".",
                                            numchr, ".annot.gz")),
                quote=FALSE, row.names=FALSE)
    cat("We are at tissue type:", bio_types[tt], "\n")
  }
  cat("We are at chr:", numchr, "\n")
}


#################################  TF all ASB V2F features    ############################################################

ccre_features = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/cCRE/cCRE_v04_13Jun2022.txt"))

rsids_test = unique(read.table(paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_tested_adastra.txt"), header=F)[,1])

snp_tabb1 = data.frame(fread(paste0(adastra_dir, "/", "AllCL_ADASTRA_sig_SNPs_FDR10_refalt.txt")))
ll = unique(snp_tabb1$CL)
annotpool = c()
biosamples = c()
for(numl in 1:length(ll)){
  rsids = snp_tabb1$ID[snp_tabb1$CL == ll[numl]]
  if(length(rsids) > 1000){
    annot = rep(length(rsids_test)/nrow(ccre_features), nrow(ccre_features))
    annot[match(intersect(ccre_features$SNP, rsids_test), ccre_features$SNP)] = 0
    annot[match(intersect(ccre_features$SNP, rsids), ccre_features$SNP)] = 1
    biosamples = c(biosamples, ll[numl])
    annotpool = cbind(annotpool, annot)
  }
  cat("We are at biosample:", numl, "\n")
}

bio_types = as.character(sapply(as.character(sapply(list.files(paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/"), pattern = "rsids_FDR"),
                                                    function(x) return(strsplit(x, "rsids_FDR10_")[[1]][2]))), function(x) return(strsplit(x, ".txt")[[1]][1])))

for(tt in 1:length(bio_types)){
  rsids = read.table(paste0(adastra_dir, "/", "Celltypes_P_sigs_1KG/",
                               "rsids_FDR10_", bio_types[tt], ".txt"), header=F)[,1]
  annot = rep(length(rsids_test)/nrow(ccre_features), nrow(ccre_features))
  annot[match(intersect(ccre_features$SNP, rsids_test), ccre_features$SNP)] = 0
  annot[match(intersect(ccre_features$SNP, rsids), ccre_features$SNP)] = 1
  annotpool = cbind(annotpool, annot)
  biosamples = c(biosamples, bio_types[tt])
}

colnames(annotpool) = biosamples

annotpool2 = cbind.data.frame(ccre_features[,1:4], annotpool)
fwrite(annotpool2, file = "/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/TF_ASB_allelic_effect_features_Aug032022.txt",
       row.names=F, col.names=T, sep = "\t", quote=F)


