
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

# vierstra_dir = "/n/groups/price/kushal/ENCODE/data/Vierstra/Version2"
#biosample_names = as.character(sapply(list.files(paste0(vierstra_dir), pattern="refalt"), function(x) return(strsplit(x, ".refalt.txt")[[1]][1])))
#write.table(biosample_names, file = paste0(vierstra_dir, "/", "sample_names2.txt"), col.names = F, row.names = F, sep = "\t", quote=F)

# biosample= "AG17017"

# vierstra_dir = "/n/groups/price/kushal/ENCODE/data/Vierstra/Version1"
# biosample_names = as.character(sapply(list.files(paste0(vierstra_dir), pattern="refalt"), function(x) return(strsplit(x, ".refalt.txt")[[1]][1])))
# write.table(biosample_names, file = paste0(vierstra_dir, "/", "sample_names2.txt"), col.names = F, row.names = F, sep = "\t", quote=F)

# biosample = "fSpinal_cord-DS20351"

tabb = data.frame(fread(paste0(vierstra_dir, "/", biosample, ".refalt.txt"), header=F))
idx = which(tabb[,10] > 35)
tabb = tabb[idx, ]
alpha_allelic_ratio = 12.4
beta_allelic_ratio = 105

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
rsids = tabb$V4[which(p.adjust(p_imb, method = "BH") < 0.10)]
rsids = rsids[grep("rs", rsids)]

if(!dir.exists(paste0(vierstra_dir, "/Celltypes_P_sigs/"))){
  dir.create(paste0(vierstra_dir, "/Celltypes_P_sigs/"))
}

write.table(rsids, file = paste0(vierstra_dir, "/Celltypes_P_sigs/", "rsids_FDR10_", biosample, ".txt"),
            row.names = F, col.names=F, sep = "\t", quote=F)

rsids = tabb$V4[which(p.adjust(p_imb, method = "BH") < 0.20)]
rsids = rsids[grep("rs", rsids)]

write.table(rsids, file = paste0(vierstra_dir, "/Celltypes_P_sigs/", "rsids_FDR20_", biosample, ".txt"),
            row.names = F, col.names=F, sep = "\t", quote=F)

