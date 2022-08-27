

#################################  Curate fine-mapped SNP lists for taskfiles   #########################################################################

ll = list.files("/n/groups/price/kushal/ENCODE/data/GWAS_traits/Finemap_SNPsets_PIP90")
snpset = c()
for(numl in 1:length(ll)){
  tabb = read.table(paste0("/n/groups/price/kushal/ENCODE/data/GWAS_traits/Finemap_SNPsets_PIP90/", ll[numl]), header=F)
  if(nrow(tabb) > 0){
    snpset = c(snpset, tabb[,1])
  }
}
snpset = unique(snpset)
write.table(snpset, file = paste0("/n/groups/price/kushal/ENCODE/data/GWAS_traits/Finemap_SNPsets", "/",
                                  "ALL_combined_PIP0_90.txt"),
            row.names=F, col.names=F, sep = "\t", quote=F)



#########################  cV2F runs using cell-line specific features   ###############################################################################

cell_lines = c("LCL", "K562", "SKNSH", "A549", "HepG2")
traits = c("RBC", "Lym", "FEV1FVC", "TBil", "Brain", "Blood")
dff = c()
for(numc in 1:length(cell_lines)){
  for(numt in 1:length(traits)){
    dff = rbind(dff, c(paste0(traits[numt], "_PIP0_1.txt"),
            "Notfinemapped_PIP0_001.txt",
            paste0(cell_lines[numc], "_features_Aug032022.txt"),
            paste0("pos_", traits[numt], "_PIP0_1_neg_Notfinemapped_PIP0_001_", cell_lines[numc], ".metrics"),
            paste0("pos_", traits[numt], "_PIP0_1_neg_Notfinemapped_PIP0_001_", cell_lines[numc], ".cv2f.txt")))
  }
}


for(numc in 1:length(cell_lines)){
  for(numt in 1:length(traits)){
    dff = rbind(dff, c(paste0(traits[numt], "_PIP0_1.txt"),
                       "ALL_combined_PIP0_75.txt",
                       paste0(cell_lines[numc], "_features_Aug032022.txt"),
                       paste0("pos_", traits[numt], "_PIP0_1_neg_ALL_combined_PIP0_75.txt_", cell_lines[numc], ".metrics"),
                       paste0("pos_", traits[numt], "_PIP0_1_neg_ALL_combined_PIP0_75.txt_", cell_lines[numc], ".cv2f.txt")))
  }
}

write.table(dff, file = c("/n/groups/price/kushal/ENCODE/code/Figure8/cV2F/celltypecV2F_combinations.txt"),
            col.names = F, row.names = F, sep = "\t", quote=F)


#########################  cV2F runs using all features   ###############################################################################

cell_lines = c("LCL", "K562", "SKNSH", "A549", "HepG2")
traits = c("RBC", "Lym", "FEV1FVC", "TBil", "Brain", "Blood")
dff = c()
for(numc in 1:length(cell_lines)){
  for(numt in 1:length(traits)){
    dff = rbind(dff, c(paste0(traits[numt], "_PIP0_1.txt"),
                       "Notfinemapped_PIP0_001.txt",
                       "Allfeatures_cCRE_V2F_perchr.",
                       paste0("pos_", traits[numt], "_PIP0_1_neg_Notfinemapped_PIP0_001_", cell_lines[numc], ".metrics"),
                       paste0("pos_", traits[numt], "_PIP0_1_neg_Notfinemapped_PIP0_001_", cell_lines[numc], ".cv2f.txt")))
  }
}

dff = rbind(dff, c("ALL_combined_PIP0_90.txt",
             "Notfinemapped_PIP0_001.txt",
             "Allfeatures_cCRE_V2F_perchr.",
             "pos_ALL_combined_PIP0_90_neg_Notfinemapped_PIP0_001.metrics",
             "pos_ALL_combined_PIP0_90_neg_Notfinemapped_PIP0_001.cv2f.txt"),
           c("ALL_combined_PIP0_75.txt",
             "Notfinemapped_PIP0_001.txt",
             "Allfeatures_cCRE_V2F_perchr.",
             "pos_ALL_combined_PIP0_75_neg_Notfinemapped_PIP0_001.metrics",
             "pos_ALL_combined_PIP0_75_neg_Notfinemapped_PIP0_001.cv2f.txt"),
           c("ALL_combined_PIP0_1.txt",
             "Notfinemapped_PIP0_001.txt",
             "Allfeatures_cCRE_V2F_perchr.",
             "pos_ALL_combined_PIP0_1_neg_Notfinemapped_PIP0_001.metrics",
             "pos_ALL_combined_PIP0_1_neg_Notfinemapped_PIP0_001.cv2f.txt"))


for(numc in 1:length(cell_lines)){
  for(numt in 1:length(traits)){
    dff = rbind(dff, c(paste0(traits[numt], "_PIP0_1.txt"),
                       "ALL_combined_PIP0_75.txt",
                       "Allfeatures_cCRE_V2F_perchr.",
                       paste0("pos_", traits[numt], "_PIP0_1_neg_ALL_combined_PIP0_75_", cell_lines[numc], ".metrics"),
                       paste0("pos_", traits[numt], "_PIP0_1_neg_ALL_combined_PIP0_75_", cell_lines[numc], ".cv2f.txt")))
  }
}


write.table(dff, file = c("/n/groups/price/kushal/ENCODE/code/Figure8/cV2F/cV2F_combinations.txt"),
            col.names = F, row.names = F, sep = "\t", quote=F)


