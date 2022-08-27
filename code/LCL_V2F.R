
boost_features = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/Allele_effects_all_features_v04_6Nov2021.txt"))
dff = data.frame(fread("/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/MPRA_GTEX_allelic_effect_features_Aug032022.txt"))
outdf  = cbind.data.frame(dff[,1:4], boost_features[, grep("LCL", colnames(boost_features))])

fwrite(outdf, file = "/n/groups/price/kushal/ENCODE/data/Deliverables/AllV2F/LCL_AFR_Montgomery_allelic_effect_features_Aug032022.txt",
       row.names=F, col.names=T, sep = "\t", quote=F)




