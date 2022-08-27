

dff = rbind(


            c("RBC_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "K562_features_Aug032022.txt", "RBC_vs_Notfinemapped_K562"),
            c("WBC+Lym_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "K562_features_Aug032022.txt", "WBC+Lym_vs_Notfinemapped_K562"),
            c("RBC_PIP0_1.txt", "NonBlood_rsids_PIP0_1.txt", "K562_features_Aug032022.txt", "RBC_vs_NonBlood_K562"),
            c("WBC+Lym_PIP0_1.txt", "NonBlood_rsids_PIP0_1.txt", "K562_features_Aug032022.txt", "WBC+Lym_vs_NonBlood_K562"),
            c("RBC_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "K562_features_Aug032022.txt", "RBC_vs_ALL_K562"),
            c("WBC+Lym_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "K562_features_Aug032022.txt", "WBC+Lym_vs_ALL_K562"),

            c("RBC_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "LCL_features_Aug032022.txt", "RBC_vs_Notfinemapped_LCL"),
            c("WBC+Lym_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "LCL_features_Aug032022.txt", "WBC+Lym_vs_Notfinemapped_LCL"),
            c("RBC_PIP0_1.txt", "NonBlood_rsids_PIP0_1.txt", "LCL_features_Aug032022.txt", "RBC_vs_NonBlood_LCL"),
            c("WBC+Lym_PIP0_1.txt", "NonBlood_rsids_PIP0_1.txt", "LCL_features_Aug032022.txt", "WBC+Lym_vs_NonBlood_LCL"),
            c("RBC_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "LCL_features_Aug032022.txt", "RBC_vs_ALL_LCL"),
            c("WBC+Lym_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "LCL_features_Aug032022.txt", "WBC+Lym_vs_ALL_LCL"),

            c("RBC_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "HepG2_features_Aug032022.txt", "RBC_vs_Notfinemapped_HepG2"),
            c("WBC+Lym_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "HepG2_features_Aug032022.txt", "WBC+Lym_vs_Notfinemapped_HepG2"),
            c("RBC_PIP0_1.txt", "NonBlood_rsids_PIP0_1.txt", "HepG2_features_Aug032022.txt", "RBC_vs_NonBlood_HepG2"),
            c("WBC+Lym_PIP0_1.txt", "NonBlood_rsids_PIP0_1.txt", "HepG2_features_Aug032022.txt", "WBC+Lym_vs_NonBlood_HepG2"),
            c("RBC_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "HepG2_features_Aug032022.txt", "RBC_vs_ALL_HepG2"),
            c("WBC+Lym_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "HepG2_features_Aug032022.txt", "WBC+Lym_vs_ALL_HepG2"),


            c("RBC_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "A549_features_Aug032022.txt", "RBC_vs_Notfinemapped_A549"),
            c("WBC+Lym_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "A549_features_Aug032022.txt", "WBC+Lym_vs_Notfinemapped_A549"),
            c("RBC_PIP0_1.txt", "NonBlood_rsids_PIP0_1.txt", "A549_features_Aug032022.txt", "RBC_vs_NonBlood_A549"),
            c("WBC+Lym_PIP0_1.txt", "NonBlood_rsids_PIP0_1.txt", "A549_features_Aug032022.txt", "WBC+Lym_vs_NonBlood_A549"),
            c("RBC_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "A549_features_Aug032022.txt", "RBC_vs_ALL_A549"),
            c("WBC+Lym_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "A549_features_Aug032022.txt", "WBC+Lym_vs_ALL_A549"),


            c("RBC_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "SKNSH_features_Aug032022.txt", "RBC_vs_Notfinemapped_SKNSH"),
            c("WBC+Lym_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "SKNSH_features_Aug032022.txt", "WBC+Lym_vs_Notfinemapped_SKNSH"),
            c("RBC_PIP0_1.txt", "NonBlood_rsids_PIP0_1.txt", "SKNSH_features_Aug032022.txt", "RBC_vs_NonBlood_SKNSH"),
            c("WBC+Lym_PIP0_1.txt", "NonBlood_rsids_PIP0_1.txt", "SKNSH_features_Aug032022.txt", "WBC+Lym_vs_NonBlood_SKNSH"),
            c("RBC_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "SKNSH_features_Aug032022.txt", "RBC_vs_ALL_SKNSH"),
            c("WBC+Lym_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "SKNSH_features_Aug032022.txt", "WBC+Lym_vs_ALL_SKNSH"),




            c("Brain_rsids_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "SKNSH_features_Aug032022.txt", "Brain_vs_Notfinemapped_SKNSH"),
            c("Brain_rsids_PIP0_1.txt", "Blood_rsids_PIP0_1.txt", "SKNSH_features_Aug032022.txt", "Brain_vs_Blood_SKNSH"),
            c("Brain_rsids_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "SKNSH_features_Aug032022.txt", "Brain_vs_ALL_SKNSH"),
            c("Brain_rsids_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "A549_features_Aug032022.txt", "Brain_vs_Notfinemapped_A549"),
            c("Brain_rsids_PIP0_1.txt", "Blood_rsids_PIP0_1.txt", "A549_features_Aug032022.txt", "Brain_vs_Blood_A549"),
            c("Brain_rsids_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "A549_features_Aug032022.txt", "Brain_vs_ALL_A549"),
            c("Brain_rsids_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "HepG2_features_Aug032022.txt", "Brain_vs_Notfinemapped_HepG2"),
            c("Brain_rsids_PIP0_1.txt", "Blood_rsids_PIP0_1.txt", "HepG2_features_Aug032022.txt", "Brain_vs_Blood_HepG2"),
            c("Brain_rsids_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "HepG2_features_Aug032022.txt", "Brain_vs_ALL_HepG2"),
            c("Brain_rsids_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "K562_features_Aug032022.txt", "Brain_vs_Notfinemapped_K562"),
            c("Brain_rsids_PIP0_1.txt", "Blood_rsids_PIP0_1.txt", "K562_features_Aug032022.txt", "Brain_vs_Blood_K562"),
            c("Brain_rsids_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "K562_features_Aug032022.txt", "Brain_vs_ALL_K562"),
            c("Brain_rsids_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "LCL_features_Aug032022.txt", "Brain_vs_Notfinemapped_LCL"),
            c("Brain_rsids_PIP0_1.txt", "Blood_rsids_PIP0_1.txt", "LCL_features_Aug032022.txt", "Brain_vs_Blood_LCL"),
            c("Brain_rsids_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "LCL_features_Aug032022.txt", "Brain_vs_ALL_LCL"),

            c("TBil_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "HepG2_features_Aug032022.txt", "TBil_vs_Notfinemapped_HepG2"),
            c("TBil_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "HepG2_features_Aug032022.txt", "TBil_vs_ALL_HepG2"),
            c("TBil_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "K562_features_Aug032022.txt", "TBil_vs_Notfinemapped_K562"),
            c("TBil_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "K562_features_Aug032022.txt", "TBil_vs_ALL_K562"),
            c("TBil_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "LCL_features_Aug032022.txt", "TBil_vs_Notfinemapped_LCL"),
            c("TBil_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "LCL_features_Aug032022.txt", "TBil_vs_ALL_LCL"),
            c("TBil_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "SKNSH_features_Aug032022.txt", "TBil_vs_Notfinemapped_SKNSH"),
            c("TBil_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "SKNSH_features_Aug032022.txt", "TBil_vs_ALL_SKNSH"),
            c("TBil_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "A549_features_Aug032022.txt", "TBil_vs_Notfinemapped_A549"),
            c("TBil_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "A549_features_Aug032022.txt", "TBil_vs_ALL_A549"),


            c("FEV1FVC_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "HepG2_features_Aug032022.txt", "FEV1FVC_vs_Notfinemapped_HepG2"),
            c("FEV1FVC_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "HepG2_features_Aug032022.txt", "FEV1FVC_vs_ALL_HepG2"),
            c("FEV1FVC_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "K562_features_Aug032022.txt", "FEV1FVC_vs_Notfinemapped_K562"),
            c("FEV1FVC_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "K562_features_Aug032022.txt", "FEV1FVC_vs_ALL_K562"),
            c("FEV1FVC_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "LCL_features_Aug032022.txt", "FEV1FVC_vs_Notfinemapped_LCL"),
            c("FEV1FVC_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "LCL_features_Aug032022.txt", "FEV1FVC_vs_ALL_LCL"),
            c("FEV1FVC_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "SKNSH_features_Aug032022.txt", "FEV1FVC_vs_Notfinemapped_SKNSH"),
            c("FEV1FVC_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "SKNSH_features_Aug032022.txt", "FEV1FVC_vs_ALL_SKNSH"),
            c("FEV1FVC_PIP0_1.txt", "Notfinemapped_PIP0_001.txt", "A549_features_Aug032022.txt", "FEV1FVC_vs_Notfinemapped_A549"),
            c("FEV1FVC_PIP0_1.txt", "ALL_combined_rsids_PIP0_1.txt", "A549_features_Aug032022.txt", "FEV1FVC_vs_ALL_A549")

)

write.table(dff, file = c("/n/groups/price/kushal/ENCODE/data/celltypeV2F_combinations.txt"),
            col.names = F, row.names = F, sep = "\t", quote=F)

ll = list.files("/n/groups/price/kushal/ENCODE/data/AUC_cV2F")
oo_pool = c()
for(numl in 1:length(ll)){
  oo = get(load(paste0("/n/groups/price/kushal/ENCODE/data/AUC_cV2F", "/", ll[numl])))
  oo_pool = rbind(oo_pool, oo[1, ])
}

rownames(oo_pool) = ll



