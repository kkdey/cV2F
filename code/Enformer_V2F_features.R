
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

target_ids = 1:5313 ## all tissues
blood_target_ids = unique(unlist(sapply(c("K562", "CD", "GM12878", "blood", "Blood", "Lymph", "lymph"),
                                   function(x) return(grep(x, target_labels)))))  ## blood tissues
brain_target_ids = unique(unlist(sapply(c("hippo", "astrocyt", "cereb", "corte", "brain", "Brain", "neur", "Neur"),
                                function(x) return(grep(x, target_labels)))))  ## brain tissues
heart_target_ids = unique(unlist(sapply(c("Heart", "heart", "cardi", "Cardi"),
                                        function(x) return(grep(x, target_labels)))))
lung_target_ids = unique(unlist(sapply(c("Lung", "lung", "pulmo"),
                                        function(x) return(grep(x, target_labels)))))
liver_target_ids = unique(unlist(sapply(c("Liver", "liver", "hepato", "HepG2"),
                                       function(x) return(grep(x, target_labels)))))
skin_target_ids = unique(unlist(sapply(c("Skin", "skin", "kerato", "Kerato"),
                                        function(x) return(grep(x, target_labels)))))
fat_target_ids = unique(unlist(sapply(c("adipo", "Adipo", "fat", "Fat"),
                                       function(x) return(grep(x, target_labels)))))
kidney_target_ids = setdiff(unique(unlist(sapply(c("Kidney", "kidney", "renal", "Renal"),
                                      function(x) return(grep(x, target_labels))))),
                            unique(unlist(sapply(c("Kidney", "kidney", "adrenal", "Adrenal"),
                                                 function(x) return(grep(x, target_labels))))))
gut_target_ids = unique(unlist(sapply(c("Intestine", "intestine", "bowel", "Bowel", "colon", "Colon", "Stomach", "stomach", "gut", "Gut"),
                                      function(x) return(grep(x, target_labels)))))

target_ids_list = list()
target_ids_list[["all"]] = target_ids
target_ids_list[["blood"]] = blood_target_ids
target_ids_list[["brain"]] = brain_target_ids
target_ids_list[["heart"]] = heart_target_ids
target_ids_list[["lung"]] = lung_target_ids
target_ids_list[["liver"]] = liver_target_ids
target_ids_list[["kidney"]] = kidney_target_ids
target_ids_list[["skin"]] = skin_target_ids
target_ids_list[["fat"]] = fat_target_ids
target_ids_list[["gut"]] = gut_target_ids

mode_ids_list = list()
mode_ids_list[["all"]] = 1:5313
mode_ids_list[["DNASE"]] = grep("DNASE", target_labels)
mode_ids_list[["CHIP"]] = grep("CHIP", target_labels)
mode_ids_list[["H3K4me3"]] = grep("H3K4me3", target_labels)
mode_ids_list[["H3K4me1"]] = grep("H3K4me1", target_labels)
mode_ids_list[["H3K27ac"]] = grep("H3K27ac", target_labels)

combo_ids_list = list()
kk=1
bios = c()
for(mm in 1:length(target_ids_list)){
  for(nn in 1:length(mode_ids_list)){
    combo_ids_list[[kk]] = intersect(target_ids_list[[mm]], mode_ids_list[[nn]])
    kk = kk + 1
    bios = c(bios, paste0(names(target_ids_list)[mm], "_", names(mode_ids_list)[nn]))
  }
}
names(combo_ids_list) = bios
combo_ids_list = combo_ids_list[which(sapply(combo_ids_list, length) > 10)]


NUMCHUNKS=100
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
all_chunks = chunk2(1:length(chr_vals), NUMCHUNKS)

meanz2_scores = vector(mode = "list", NUMCHUNKS)
maxz2_scores = vector(mode = "list", NUMCHUNKS)

for(numl in 1:NUMCHUNKS){
  tempp = h5read(paste0(enformer_dir, "/",
                        "1000G.MAF_threshold=0.005.", numchr, ".h5"), name="SAD", start = c(1, head(all_chunks[[numl]], 1)),
                 count = c(5313, length(all_chunks[[numl]])))
  mean_vecc = rowMeans(tempp)
  sd_vecc = apply(tempp, 1, sd)
  mean_matt = mean_vecc %*% t(rep(1, ncol(tempp)))
  sd_matt = sd_vecc %*% t(rep(1, ncol(tempp)))
  zmatt = (tempp - mean_matt)/sd_matt
  rm(sd_matt)
  rm(mean_matt)
  temp1 = c()
  temp2 = c()
  for(cc in 1:length(combo_ids_list)){
    zmatt_temp = zmatt[combo_ids_list[[cc]], ]
    temp1 = cbind(temp1, apply(abs(zmatt_temp), 2, mean))
    temp2 = cbind(temp2, apply(abs(zmatt_temp), 2, max))
  }
  meanz2_scores[[numl]] = temp1
  maxz2_scores[[numl]] = temp2
  cat("We are at chunk:", numl, "\n")
}

meanz2_scores_vec = do.call(rbind, meanz2_scores)
maxz2_scores_vec = do.call(rbind, maxz2_scores)
celltypes_mode = names(combo_ids_list)

for(numc in 1:length(celltypes_mode)){
  merged_outdf = cbind.data.frame(chr_vals, pos_vals, snp_vals, ref_vals, alt_vals, meanz2_scores_vec[,numc], maxz2_scores_vec[,numc])
  colnames(merged_outdf) = c("CHR", "BP", "SNP", "REF", "ALT", "MEANZ2", "MAXZ2")
  fwrite(merged_outdf, file = paste0(output_dir, "/",
                                     "chr", numchr,".", celltypes_mode[numc], ".txt"),
         row.names=F, col.names=T, sep = "\t", quote=F)
}


