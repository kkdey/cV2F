library(data.table)

dff1 = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/GTEX_finemap_allelic_effect_features_Aug032022.txt"))
dff2 = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/LCL_AFR_Montgomery_allelic_effect_features_Aug032022.txt"))
dff3 = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/Enformer_allelic_effect_features_Aug032022.txt"))
dff4 = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/MPRA_GTEX_allelic_effect_features_Aug032022.txt"))
dff5 = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/TF_ASB_allelic_effect_features_Aug032022.txt"))
dff6 = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/Vierstra_v1_allelic_effect_features_Aug032022.txt"))
dff7 = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/cCRE/cCRE_celltypes_v04_13Jun2022.txt"))

for(numchr in 1:22){

  header = dff1[which(dff1$CHR == numchr), (1:4)]
  temp_dff1 = dff1[which(dff1$CHR == numchr), -(1:4)]
  temp_dff2 = dff2[which(dff2$CHR == numchr), -(1:4)]
  temp_dff3 = dff3[which(dff3$CHR == numchr), -(1:4)]
  temp_dff4 = dff4[which(dff4$CHR == numchr), -(1:4)]
  temp_dff5 = dff5[which(dff5$CHR == numchr), -(1:4)]
  temp_dff6 = dff6[which(dff6$CHR == numchr), -(1:4)]
  temp_dff7 = dff7[which(dff7$CHR == numchr), -(1:4)]

  merged_dff = cbind.data.frame(header,
                                temp_dff1,
                                temp_dff2,
                                temp_dff3,
                                temp_dff4,
                                temp_dff5,
                                temp_dff6,
                                temp_dff7)
  fwrite(merged_dff, file = paste0("/n/groups/price/kushal/ENCODE/data/Deliverables/Feature_Tables/Allfeatures_cCRE_V2F_perchr",
                                   ".", numchr, ".txt"), row.names=F, col.names=T, sep = "\t", quote=F)
  rm(temp_dff1)
  rm(temp_dff2)
  rm(temp_dff3)
  rm(temp_dff4)
  rm(temp_dff5)
  rm(temp_dff6)
  rm(temp_dff7)
  rm(merged_dff)
  cat("We are at chr:", numchr, "\n")

}







