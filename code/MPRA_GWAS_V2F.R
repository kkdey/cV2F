
library(data.table)
library(R.utils)
library(GenomicRanges)

option_list <- list(
  make_option("--mpra_dir", type="character", default = "data/mpra_gwas", help="Directory where you have downloaded the MPRA allelic imbalance files"),
  make_option("--bimdir", type="character", default = "data/mpra_gwas", help="bim file directory"),
  make_option("--out_file", type="character", default = "data/mpra_gwas", help="Output MPRA V2F file")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)
mpra_dir = opt$mpra_dir
bimdir = opt$bimdir

# mpra_dir = "/n/groups/price/kushal/ENCODE/data/MPRA"
# bimdir = "/n/groups/price/kushal/extras/BIMS"

celltypes = c("A549", "GM12878", "HEPG2", "K562", "SKNSH")
rsids_vec = c()
for(celltype in celltypes){
  ll = list.files(mpra_dir, pattern = celltype)
  snp_ids_all = c()
  for(numl in 1:length(ll)){
    tabb = read.table(paste0(mpra_dir, "/", ll[numl]), header=T)
    tabb2 = tabb[which(tabb$Skew.logFDR > -log(0.10)), ]
    snp_ids_all = c(snp_ids_all, tabb2$SNP)
  }

  chrs = paste0("chr", as.character(sapply(snp_ids_all, function(x) strsplit(x, ":")[[1]][1])))
  pos = as.numeric(sapply(snp_ids_all, function(x) strsplit(x, ":")[[1]][2]))


  rsids = rep(NA, length(snp_ids_all))
  for(numchr in 1:22){
    bimtabb = read.table(paste0(bimdir, "/", "1000G.EUR.QC.", numchr, ".bim"))

    gr1 = GRanges(seqnames = Rle(chrs),
                  ranges = IRanges(start=pos-1, end = pos+1))
    gr2 = GRanges(seqnames = Rle(paste0("chr", bimtabb[,1])),
                  ranges = IRanges(start=bimtabb[,4]-1, end = bimtabb[,4]+1))
    cc=findOverlaps(gr2, gr1)
    rsids[subjectHits(cc)] = bimtabb[queryHits(cc), 2]
    cat("We are at chr:", numchr, "\n")
  }

  rsids2 = rsids[!is.na(rsids)]
  if(!dir.exists(paste0(mpra_dir, "/", "Celltypes_P_sigs_1KG/"))){
    dir.create(paste0(mpra_dir, "/", "Celltypes_P_sigs_1KG/"))
  }

  rsids_vec = c(rsids_vec, rsids2)
  write.table(rsids2, file = paste0(mpra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_", celltype, ".txt"),
              quote=F, sep = "\t", col.names=F, row.names=F)
  cat("We are at celltype:", celltype, "\n")
}

rsids_vec = unique(rsids_vec)
write.table(rsids_vec, file = paste0(mpra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_all.txt"),
            quote=F, sep = "\t", col.names=F, row.names=F)



###################################  List all tested SNPs   ######################################################################

celltypes = c("A549", "GM12878", "HEPG2", "K562", "SKNSH")
rsids_vec = c()
for(celltype in celltypes){
  ll = list.files(mpra_dir, pattern = celltype)
  snp_ids_all = c()
  for(numl in 1:length(ll)){
    tabb = read.table(paste0(mpra_dir, "/", ll[numl]), header=T)
    snp_ids_all = c(snp_ids_all, tabb$SNP)
  }

  chrs = paste0("chr", as.character(sapply(snp_ids_all, function(x) strsplit(x, ":")[[1]][1])))
  pos = as.numeric(sapply(snp_ids_all, function(x) strsplit(x, ":")[[1]][2]))


  rsids = rep(NA, length(snp_ids_all))
  for(numchr in 1:22){
    bimtabb = read.table(paste0(bimdir, "/", "1000G.EUR.QC.", numchr, ".bim"))

    gr1 = GRanges(seqnames = Rle(chrs),
                  ranges = IRanges(start=pos-1, end = pos+1))
    gr2 = GRanges(seqnames = Rle(paste0("chr", bimtabb[,1])),
                  ranges = IRanges(start=bimtabb[,4]-1, end = bimtabb[,4]+1))
    cc=findOverlaps(gr2, gr1)
    rsids[subjectHits(cc)] = bimtabb[queryHits(cc), 2]
    cat("We are at chr:", numchr, "\n")
  }

  rsids2 = rsids[!is.na(rsids)]
  rsids_vec = c(rsids_vec, rsids2)
  cat("We are at celltype:", celltype, "\n")
}

rsids_vec = unique(rsids_vec)
write.table(rsids_vec, file = paste0(mpra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_tested.txt"),
            quote=F, sep = "\t", col.names=F, row.names=F)


################################  Generate MPRA GWAS V2F features   ##############################################################

ccre_features = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/cCRE/cCRE_v04_13Jun2022.txt"))
rsids_test = unique(read.table(paste0(mpra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_tested.txt"), header=F)[,1])
annotpool = c()
biosamples = c()
celltypes = as.character(sapply(as.character(sapply(list.files(paste0(mpra_dir, "/", "Celltypes_P_sigs_1KG/"), pattern = "rsids_FDR"),
                                                    function(x) return(strsplit(x, "rsids_FDR10_")[[1]][2]))), function(x) return(strsplit(x, ".txt")[[1]][1])))

for(numl in 1:length(celltypes)){
  rsids = unique(read.table(paste0(mpra_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_", celltypes[numl], ".txt"), header=F)[,1])
  annot = rep(length(rsids_test)/nrow(ccre_features), nrow(ccre_features))
  annot[match(intersect(ccre_features$SNP, rsids_test), ccre_features$SNP)] = 0
  annot[match(intersect(ccre_features$SNP, rsids), ccre_features$SNP)] = 1
  biosamples = c(biosamples, celltypes[numl])
  annotpool = cbind(annotpool, annot)
  cat("We are at biosample:", numl, "\n")
}

colnames(annotpool) = celltypes
annotpool2 = cbind.data.frame(ccre_features[,1:4], annotpool)
fwrite(annotpool2, file = "/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/MPRA_GWAS_allelic_effect_features_Aug032022.txt",
       row.names=F, col.names=T, sep = "\t", quote=F)

