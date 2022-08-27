library(data.table)
library(R.utils)

option_list <- list(
  make_option("--vierstra_dir", type="character", default = "data/vierstra", help="Directory where you have downloaded the Vierstra files"),
  make_option("--biosample", type="character",  help="Name of the biosample for which to generate annotation")
)


opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

vierstra_dir = opt$vierstra_dir
biosample = opt$biosample

bimfiles = list.files("/n/groups/price/kushal/extras/BIMS_hg38")
bimpooltabb = c()
for(numchr in 1:22){
  bimtabb = read.table(paste0("/n/groups/price/kushal/extras/BIMS_hg38", "/", "1000G.EUR.QC.", numchr, ".bim"))
  bimpooltabb = rbind(bimpooltabb, bimtabb[, c(1, 4, 2)])
  cat("We are at chr:", numchr, "\n")
}
colnames(bimpooltabb) = c("CHR", "BP", "SNP")


#vierstra_dir = "/n/groups/price/kushal/ENCODE/data/Vierstra/Version1"
#biosample_names = as.character(sapply(list.files(paste0(vierstra_dir), pattern="refalt"), function(x) return(strsplit(x, ".refalt.txt")[[1]][1])))


tabb = data.frame(fread(paste0(vierstra_dir, "/", biosample, ".refalt.txt"), header=F))
idx = which(tabb[,10] > 30)
tabb = tabb[idx, ]
alpha_allelic_ratio = 75
beta_allelic_ratio = 75

p_imb = c()
for(numl in 1:nrow(tabb)){
  xx = tabb[numl, ]
  p_imb = c(p_imb, pmin(1 - pbb(as.numeric(xx[10] - xx[9]), as.numeric(xx[10]), alpha_allelic_ratio, beta_allelic_ratio),
                        1 - pbb(min(as.numeric(xx[9]), as.numeric(xx[10]) - 1),
                                as.numeric(xx[10]),
                                beta_allelic_ratio,
                                alpha_allelic_ratio)))
}
tabb$pval = p_imb

tempp = tabb[which(p.adjust(p_imb, method = "BH") < 0.20), ]

gr1 = GRanges(seqnames = Rle(tempp[,1]),
              ranges = IRanges(start=tempp[,2]-1, end = tempp[,3]+1))
gr2 = GRanges(seqnames = Rle(paste0("chr", bimpooltabb[,1])),
              ranges = IRanges(start=bimpooltabb[,2]-1, end = bimpooltabb[,2]+1))

cc=findOverlaps(gr2, gr1)

rsids = unique(bimpooltabb[queryHits(cc), 3])
write.table(rsids, file =  paste0(vierstra_dir, "/Celltypes_P_sigs/", "rsids_FDR20_", biosample, ".txt"),
            row.names = F, col.names=F, sep = "\t", quote=F)


tempp = tabb[which(p.adjust(p_imb, method = "BH") < 0.10), ]

gr1 = GRanges(seqnames = Rle(tempp[,1]),
              ranges = IRanges(start=tempp[,2]-1, end = tempp[,3]+1))
gr2 = GRanges(seqnames = Rle(paste0("chr", bimpooltabb[,1])),
              ranges = IRanges(start=bimpooltabb[,2]-1, end = bimpooltabb[,2]+1))

cc=findOverlaps(gr2, gr1)

rsids = unique(bimtabb[queryHits(cc), 2])
write.table(rsids, file =  paste0(vierstra_dir, "/Celltypes_P_sigs/", "rsids_FDR10_", biosample, ".txt"),
            row.names = F, col.names=F, sep = "\t", quote=F)


cat("We are at biosample:", biosample, "\n")
