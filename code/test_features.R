
celltypes = c("LCL", "K562", "SKNSH", "A549", "HepG2")
traitfiles = c( "RBC_PIP0_1.txt", "TBil_PIP0_1.txt", "Notfinemapped_PIP0_001.txt")

EEmat = matrix(0, length(celltypes), length(traitfiles))
for(numl in 1:length(celltypes)){
  feature_file = paste0("/n/groups/price/kushal/ENCODE/data/Deliverables/Feature_Tables/", celltypes[numl],
                 "_features_Aug032022.txt")
  data_features = data.frame(fread(feature_file))

  for(numt in 1:length(traitfiles)){
    positive_set = paste0("/n/groups/price/kushal/ENCODE/data/GWAS_traits/Finemap_SNPsets/", traitfiles[numt])
    positive_snps = read.table(positive_set, header=F)[,1]
    positive_snps = intersect(positive_snps, data_features$SNP)
    EEmat[numl, numt] = mean(colMeans(data_features[match(positive_snps, data_features$SNP), -(1:4)])/colMeans(data_features[, -(1:4)]))
  }
  cat("we are at cell type:", numl, "\n")
}

rownames(EEmat) = celltypes
colnames(EEmat) = traitfiles

