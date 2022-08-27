

dff = data.frame(fread(paste0("/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/MPRA_GWAS_allelic_effect_features_Aug032022.txt")))
dff2 = dff[,-(1:4)]
ll=list.files("/n/groups/price/kushal/ENCODE/data/Deliverables/cV2F_hg38", pattern = ".cv2f.txt")

EEmatt = matrix(0, length(ll), ncol(dff2))

for(numl in 1:length(ll)){
  cc = data.frame(fread(paste0("/n/groups/price/kushal/ENCODE/data/Deliverables/cV2F_hg38", "/",
                               ll[numl])))
  for(numc in 1:ncol(dff2)){
    idx1 = c(which(dff2[,numc] == 0), which(dff2[,numc] == 1))
    EEmatt[numl, numc] =  (sum(cc[idx1,5]*dff2[idx1,numc])*length(idx1))/(sum(cc[idx1,5]) * sum(dff2[idx1,numc]))
  }
  cat("We are at file:", numl, "\n")
}

rownames(EEmatt) = ll
colnames(EEmatt) = colnames(dff2)
