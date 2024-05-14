source("downloadRseData.R")
source("findSignificantEvents.R")
source('plotToolsResults.R')

logfc_threshold=1.5
dpsi_threshold=0.1
abund_change_threshold=1
fdr_threshold=0.05

thresholds = list(logfc_threshold=logfc_threshold,
                   dpsi_threshold=dpsi_threshold,
                   abund_change_threshold=abund_change_threshold,
                   fdr_threshold=fdr_threshold)

# #==================== download and filter rse
# #prepareRse()
# 
# 
# # #=======================running tools
# # # # -- reading files
rse.gene.cytosk = readRDS('rse.gene.cytosk.rds', refhook = NULL)
rse.jxn.cytosk = readRDS('rse.jxn.cytosk.rds', refhook = NULL)
unique.tissues = unique(rse.jxn.cytosk@colData$tissue)
outputs_tissue = list()
for (tissue in unique.tissues){
  outputs_tissue[[tissue]] = runTools(rse.jxn.cytosk, tissue)
}
saveRDS(outputs_tissue,'outputs_tissue.rds')
outputs_tissue = readRDS('outputs_tissue.rds', refhook = NULL)
# 
# 
# #=================================dev
outputs_tissue = readRDS('outputs_tissue.rds', refhook = NULL)
unique.tissues = unique(rse.jxn.cytosk@colData$tissue)
outputs_dev_sign_info = list()
for (tissue in unique.tissues){
  outputs_dev_sign_info[[tissue]] = getJxnSignInfo(tools.outputs.list=outputs_tissue[[tissue]],
                                                   logfc_threshold=logfc_threshold,
                                                   dpsi_threshold=dpsi_threshold,
                                                   abund_change_threshold=abund_change_threshold,
                                                   fdr_threshold=fdr_threshold,
                                                   add_external_data=FALSE, file='')
}
saveRDS(outputs_dev_sign_info,'outputs_dev_sign_info.rds')
outputs_dev_sign_info = readRDS('outputs_dev_sign_info.rds', refhook = NULL)
plotResultsRepot(outputs_dev_sign_info, thresholds = thresholds)



#================================= tumor
#-- reading files
outputs_tissue = readRDS('outputs_tissue.rds', refhook = NULL)
unique.tissues = unique(rse.jxn.cytosk@colData$tissue)
outputs_gtex2tum = list()
outputs_norm2tum = list()
for (tissue in unique.tissues){
  outputs_gtex2tum[[tissue]] = getJxnSignInfo(outputs_tissue[[tissue]],
                                              logfc_threshold=logfc_threshold,
                                              dpsi_threshold=dpsi_threshold,
                                              abund_change_threshold=abund_change_threshold,
                                              fdr_threshold=fdr_threshold,
                                              add_external_data=TRUE, file='gtex2tum')
  outputs_norm2tum[[tissue]] = getJxnSignInfo(outputs_tissue[[tissue]],
                                              logfc_threshold=logfc_threshold,
                                              dpsi_threshold=dpsi_threshold,
                                              abund_change_threshold=abund_change_threshold,
                                              fdr_threshold=fdr_threshold,
                                              add_external_data=TRUE, file='norm2tum')
}
saveRDS(outputs_gtex2tum,'dev_vs_gtex2tum_tools.rds')
saveRDS(outputs_norm2tum,'dev_vs_norm2tum_tools.rds')

outputs_gtex2tum = readRDS('dev_vs_gtex2tum_tools.rds', refhook = NULL)
outputs_norm2tum = readRDS('dev_vs_norm2tum_tools.rds', refhook = NULL)


plotResultsRepot(outputs_gtex2tum, tumor=TRUE, file='gtex2tum',
                 thresholds = list(logfc_threshold=logfc_threshold,
                                   dpsi_threshold=dpsi_threshold,
                                   abund_change_threshold=abund_change_threshold,
                                   fdr_threshold=fdr_threshold))

plotResultsRepot(outputs_norm2tum, tumor=TRUE, file='norm2tum',
                 thresholds = list(logfc_threshold=logfc_threshold,
                                   dpsi_threshold=dpsi_threshold,
                                   abund_change_threshold=abund_change_threshold,
                                   fdr_threshold=fdr_threshold))



fisher_results_tissues_list = list()
for (tissue in unique.tissues){
  fisher.df = data.frame(tool_pair=character(), odds_ratio = numeric(), p_val = numeric())
  fisher.df = getFisher(fisher.df, outputs_dev_sign_info[[tissue]], one_to_all=FALSE,  ref_col='')
  fisher.df = getFisher(fisher.df, outputs_gtex2tum[[tissue]], one_to_all=TRUE, ref_col='sajr.norm.tumor')
  fisher.df$tool_pair = gsub("sajr.norm.tumor",  "sajr.gtex2tum",  fisher.df$tool_pair)
  fisher.df = getFisher(fisher.df, outputs_norm2tum[[tissue]], one_to_all=TRUE, ref_col='sajr.norm.tumor')
  fisher.df$tool_pair = gsub("sajr.norm.tumor",  "sajr.norm2tum",  fisher.df$tool_pair)
  fisher.df$q_val = p.adjust(fisher.df$p_val, method = "BH") 
  fisher_results_tissues_list[[tissue]] = fisher.df
}
fisher_results_tissues_list

plotFisherResults(fisher_results_tissues_list, thresholds= thresholds)






