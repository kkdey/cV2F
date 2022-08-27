library(data.table)
library(R.utils)
library(GenomicRanges)

option_list <- list(
  make_option("--gtex_dir", type="character", default = "data/gtex", help="Directory where you have downloaded the GTEx eQTL finemap data"),
  make_option("--bimdir", type="character", default = "data/", help="bim file directory"),
  make_option("--out_file", type="character", default = "data/gtex", help="Output MPRA V2F file")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)
gtex_dir = opt$gtex_dir
bimdir = opt$bimdir

# gtex_dir = "/n/groups/price/kushal/ENCODE/data/GTEx_finemapping/"
# bimdir = "/n/groups/price/kushal/extras/BIMS_hg38"

library(data.table)
tabb = data.frame(fread(paste0(gtex_dir, "/", "GTEx_49tissues_release1.tsv.gz")))
tissue_types = unique(tabb$V11)
for(tt in 1:length(tissue_types)){
  tabb2 = tabb[which(tabb$V11 == tissue_types[tt]), ]
  maxpip = tapply(tabb2$V17, tabb2$V5, max)
  df = data.frame("SNP" = names(maxpip),
                  "MaxCPP" = maxpip)
  write.table(df, file = paste0(gtex_dir, "/MaxPIP_BY_TISSUE", "/",
                                tissue_types[tt], ".txt"),
              row.names = F, col.names = T, sep = "\t", quote=F)
  cat("We are at tissue:", tissue_types[tt], "\n")
}


#######################  Restrict to SNPs with rsIDs  #######################################################

bimfile_all = c()
for(numchr in 1:22){
  bimfile = read.table(paste0(bimdir, "/", "1000G.EUR.QC.", numchr, ".bim"))
  bimfile_all = rbind(bimfile_all, bimfile)
  cat("We are at chr:", numchr, "\n")
}

gr_bimfiles <- GRanges(
  seqnames = Rle(paste0("chr", bimfile_all[,1])),
  ranges = IRanges(bimfile_all[,4]-1, end=bimfile_all[,4]+1),
  rsid= bimfile_all[,2])

ll = list.files(paste0(gtex_dir, "/", "MaxPIP_BY_TISSUE"), pattern=".txt")
ll2 = as.character(sapply(ll, function(x) return(strsplit(x, ".txt")[[1]][1])))
for(numl in 1:length(ll)){
  tabb = read.table(paste0(gtex_dir, "/", "MaxPIP_BY_TISSUE", "/", ll[numl]), header=T)
  chrvals = as.character(sapply(tabb[,1], function(x) return(strsplit(x, "_")[[1]][1])))
  posvals = as.numeric(sapply(tabb[,1], function(x) return(strsplit(x, "_")[[1]][2])))

  gr2 <- GRanges(
    seqnames = Rle(chrvals),
    ranges = IRanges(posvals-1, end=posvals+1),
    pip = tabb$MaxCPP)

  overlaps = findOverlaps(gr2, gr_bimfiles)
  df = cbind.data.frame(gr_bimfiles$rsid[subjectHits(overlaps)], gr2$pip[queryHits(overlaps)])
  maxcpp_rs = tapply(df[,2], df[,1], max)

  outdf = data.frame("SNP" = names(maxcpp_rs),
                     "MaxCPP" = maxcpp_rs)
  write.table(outdf, file = paste0(gtex_dir, "/", "MaxPIP_BY_TISSUE", "/",
                                   ll2[numl], ".rsid"),
              row.names = F, col.names = T, sep = "\t", quote=F)
  cat("We are at tissue:", ll[numl], "\n")
}


###################################  Generate GTEx finemap eQTL V2F scores   ##############################################################

ccre_features = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/cCRE/cCRE_v04_13Jun2022.txt"))

ll = list.files(paste0(gtex_dir, "/", "MaxPIP_BY_TISSUE"), pattern = ".rsid")
annotpool = c()
biosamples = c()
for(numl in 1:length(ll)){
  tabb = read.table(paste0(gtex_dir, "/", "MaxPIP_BY_TISSUE", "/", ll[numl]), header=T)
  pip_cut = c(0.1, 0.5, 0.9)
  for(cc in 1:length(pip_cut)){
    rsids = tabb$SNP[tabb$MaxCPP > pip_cut[cc]]
    annot = rep(0, nrow(ccre_features))
    annot[match(intersect(ccre_features$SNP, rsids), ccre_features$SNP)] = 1
    biosamples = c(biosamples, paste0(ll[numl], "_", pip_cut[cc]))
    annotpool = cbind(annotpool, annot)
  }
  cat("We are at tissue:", numl, "\n")
}

colnames(annotpool) = biosamples
annotpool2 = cbind.data.frame(ccre_features[,1:4], annotpool)
fwrite(annotpool2, file = "/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/GTEX_finemap_allelic_effect_features_Aug032022.txt",
       row.names=F, col.names=T, sep = "\t", quote=F)
