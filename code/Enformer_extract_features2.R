library(rhdf5)
library(data.table)
library(R.utils)
library(optparse)

options(echo=TRUE) # if you want see commands in output file
args <- commandArgs(trailingOnly = TRUE)
print(args)

option_list <- list(
  make_option("--enformer_dir", type="character", default = "data/enformer", help="Directory where Enformer SAD annotations are saved"),
  make_option("--output_dir", type="character", default = "data/enformer_features", help="Directory to save Enformer feature scores for variants"),
  make_option("--numchr", type="character", help="Chromosome number")
)

opt <- parse_args(OptionParser(option_list=option_list))
dput(opt)

enformer_dir = opt$enformer_dir #/n/groups/price/kushal/Imperio/ALL_DEEP/Enformer/1000-genomes/enformer
output_dir = opt$output_dir # "/n/groups/price/kushal/Imperio/ALL_DEEP/Enformer/Z_SKEW_FEATURES"
numchr = as.integer(opt$numchr)

enformer_lists = list.files(paste0(enformer_dir), pattern = ".h5")
h5ls(paste0("/n/groups/price/kushal/Imperio/ALL_DEEP/Enformer/1000-genomes/enformer", "/",
            "1000G.MAF_threshold=0.005.", numchr, ".h5"))
chr_vals = h5read(paste0(enformer_dir, "/",
                         "1000G.MAF_threshold=0.005.", numchr, ".h5"), name="chr")
pos_vals = h5read(paste0(enformer_dir, "/",
                         "1000G.MAF_threshold=0.005.", numchr, ".h5"), name="pos")
ref_vals = h5read(paste0(enformer_dir, "/",
                         "1000G.MAF_threshold=0.005.", numchr, ".h5"), name="ref")
alt_vals = h5read(paste0(enformer_dir, "/",
                         "1000G.MAF_threshold=0.005.", numchr, ".h5"), name="alt")
snp_vals = h5read(paste0(enformer_dir, "/",
                         "1000G.MAF_threshold=0.005.", numchr, ".h5"), name="snp")
target_labels = h5read(paste0(enformer_dir, "/",
                              "1000G.MAF_threshold=0.005.", numchr, ".h5"), name="target_labels")


K562_target_ids = as.numeric(unique(unlist(sapply(c("K562"), function(x) return(grep(x, target_labels))))))
GM12878_target_ids = as.numeric(unique(unlist(sapply(c("GM12878"), function(x) return(grep(x, target_labels))))))
