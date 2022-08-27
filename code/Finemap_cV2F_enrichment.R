
traits = c("RBC", "Lym", "TBil", "FEV1FVC")

ll=list.files("/n/groups/price/kushal/ENCODE/data/Deliverables/cV2F_hg38", pattern = ".cv2f.txt")

EEmatt = matrix(0, length(ll), length(traits))

for(numt in 1:length(traits)){
  finemap_snps = read.table(paste0("/n/groups/price/kushal/ENCODE/data/GWAS_traits/Finemap_SNPsets_PIP75/",
                     traits[numt], "_PIP0_75.txt"), header=F)[,1]
  for(numl in 1:length(ll)){
    cc = data.frame(fread(paste0("/n/groups/price/kushal/ENCODE/data/Deliverables/cV2F_hg38", "/",
                                 ll[numl])))
    EEmatt[numl, numt] = mean(cc[match(intersect(finemap_snps, cc$SNP), cc$SNP), 5])/mean(cc[,5])
  }
  cat("We are at trait:", numt, "\n")
}

rownames(EEmatt) = ll
colnames(EEmatt) = traits



