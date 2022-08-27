
library(data.table)
library(R.utils)

option_list <- list(
  make_option("--vierstra_dir", type="character", default = "data/vierstra", help="Directory where you have downloaded the Vierstra allelic imbalance files"),
  make_option("--annot_dir", type="character", default = "data/annotations/vierstra", help="Directory where to save the Vierstra allelic imbalance annotations")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

#####################################  Extract all tested variants in Vierstra DNase AI   ######################################################################

vierstra_dir = opt$vierstra_dir
annot_dir = opt$annot_dir

# vierstra_dir = "/n/groups/price/kushal/ENCODE/data/Vierstra/Version2"
# vierstra_dir = "/n/groups/price/kushal/ENCODE/data/Vierstra/Version1"

#########################  Generate a combo allelic imbalance data aross all biosamples  ####################
################# This part of the code several hours to run  #########################################

ll = list.files(paste0(vierstra_dir), pattern="refalt")
merged_altref = list()
for(numl in 1:100){
  tabb = data.frame(fread(paste0(vierstra_dir, "/", ll[numl]), header=F))
  idx = which(tabb[,10] > 35)
  tabb = tabb[idx, ]
  merged_altref[[numl]] = tabb[,c(9, 10)]
  cat("We are at biosample:", numl, "\n")
}

merged_altref_tabb = do.call(rbind, merged_altref)

mean_allelic_ratio = mean((merged_altref_tabb[,2] - merged_altref_tabb[,1])/merged_altref_tabb[,2])
sd_allelic_ratio = sd((merged_altref_tabb[,2] - merged_altref_tabb[,1])/merged_altref_tabb[,2])
alpha_allelic_ratio = mean_allelic_ratio*((1/sd_allelic_ratio^2) - 1)
beta_allelic_ratio = (1 - mean_allelic_ratio)*((1/sd_allelic_ratio^2) - 1)

## alpha_allelic ratio = 12.4 (version 2)
## beta_allelic_ratio = 105

## alpha_allelic ratio = 75 (version 1)
## beta_allelic_ratio = 75


library(Biobase)
library(TailRank)
