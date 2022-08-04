library(data.table)
library(R.utils)
library(GenomicRanges)

option_list <- list(
  make_option("--mpra_gtex_dir", type="character", default = "data/mpra_gwas", help="Directory where you have downloaded the MPRA allelic imbalance files"),
  make_option("--bimdir", type="character", default = "data/mpra_gwas", help="bim file directory"),
  make_option("--out_file", type="character", default = "data/mpra_gwas", help="Output MPRA V2F file")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)
mpra_gtex_dir = opt$mpra_gtex_dir
bimdir = opt$bimdir

# mpra_gtex_dir = "/n/groups/price/kushal/ENCODE/data/MPRA_GTEx/Processed_Results/"
# bimdir = "/n/groups/price/kushal/extras/BIMS"


celltypes = c("HCT116", "HepG2", "SKNSH", "K562")
rsids_vec = c()
for(celltype in celltypes){
  ll = list.files(paste0(mpra_gtex_dir), pattern = celltype)
  snp_ids_all = c()
  for(numl in 1:length(ll)){
    tabb = read.table(paste0(mpra_gtex_dir, "/", ll[numl]), header=T)
    tabb2 = tabb[which(tabb$Skew.logFDR > -log(0.10)), ]
    snp_ids_all = c(snp_ids_all, tabb2$SNP)
  }

  chrs = paste0("chr", as.character(sapply(snp_ids_all, function(x) strsplit(x, ":")[[1]][1])))
  pos = as.numeric(sapply(snp_ids_all, function(x) strsplit(x, ":")[[1]][2]))

  library(GenomicRanges)
  rsids = rep(NA, length(snp_ids_all))
  for(numchr in 1:22){
    bimtabb = read.table(paste0(bimdir, "/", "1000G.EUR.QC.", numchr, ".bim"))
    #bimtabb = read.table(paste0("/n/groups/price/kushal/LDSC/1000G_EUR_Phase3_plink_hg38/1000G.EUR.QC.", numchr, ".bim"))
    gr1 = GRanges(seqnames = Rle(chrs),
                  ranges = IRanges(start=pos-1, end = pos+1))
    gr2 = GRanges(seqnames = Rle(paste0("chr", bimtabb[,1])),
                  ranges = IRanges(start=bimtabb[,4]-1, end = bimtabb[,4]+1))
    cc=findOverlaps(gr2, gr1)
    rsids[subjectHits(cc)] = bimtabb[queryHits(cc), 2]
    cat("We are at chr:", numchr, "\n")
  }

  rsids2 = rsids[!is.na(rsids)]
  if(!dir.exists(paste0(mpra_gtex_dir, "/", "Celltypes_P_sigs_1KG/"))){
    dir.create(paste0(mpra_gtex_dir, "/", "Celltypes_P_sigs_1KG/"))
  }

  rsids_vec = c(rsids_vec, rsids2)
  write.table(rsids2, file = paste0(mpra_gtex_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_", celltype, ".txt"),
              quote=F, sep = "\t", col.names=F, row.names=F)
  cat("We are at celltype:", celltype, "\n")
}

rsids_vec = unique(rsids_vec)
write.table(rsids_vec, file = paste0(mpra_gtex_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_FDR10_all.txt"),
            quote=F, sep = "\t", col.names=F, row.names=F)


############################  All tested variants in MPRA GTEx  #################################################

celltypes = c("HCT116", "HepG2", "SKNSH", "K562")
rsids_vec = c()
for(celltype in celltypes){
  ll = list.files(paste0(mpra_gtex_dir), pattern = celltype)
  snp_ids_all = c()
  for(numl in 1:length(ll)){
    tabb = read.table(paste0(mpra_gtex_dir, "/", ll[numl]), header=T)
    snp_ids_all = c(snp_ids_all, tabb$SNP)
  }

  chrs = paste0("chr", as.character(sapply(snp_ids_all, function(x) strsplit(x, ":")[[1]][1])))
  pos = as.numeric(sapply(snp_ids_all, function(x) strsplit(x, ":")[[1]][2]))

  library(GenomicRanges)
  rsids = rep(NA, length(snp_ids_all))
  for(numchr in 1:22){
    bimtabb = read.table(paste0(bimdir, "/", "1000G.EUR.QC.", numchr, ".bim"))
    #bimtabb = read.table(paste0("/n/groups/price/kushal/LDSC/1000G_EUR_Phase3_plink_hg38/1000G.EUR.QC.", numchr, ".bim"))
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
write.table(rsids_vec, file = paste0(mpra_gtex_dir, "/", "Celltypes_P_sigs_1KG/", "rsids_tested.txt"),
            quote=F, sep = "\t", col.names=F, row.names=F)

