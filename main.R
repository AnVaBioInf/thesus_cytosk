source("downloadRseData.R")
source("findSignificantEvents.R")
source('plotToolsResults.R')

logfc_threshold=1.5
dpsi_threshold=0.1
abund_change_threshold=0.5
fdr_threshold=0.1

#==================== download and filter rse
#prepareRse()


# #=======================running tools
# # # -- reading files
rse.gene.cytosk = readRDS('rse.gene.cytosk.rds', refhook = NULL)
rse.jxn.cytosk = readRDS('rse.jxn.cytosk.rds', refhook = NULL)
unique.tissues = unique(rse.jxn.cytosk@colData$tissue)
# outputs_tissue = list()
# for (tissue in unique.tissues){
#   outputs_tissue[[tissue]] = runTools(rse.jxn.cytosk, tissue)
# }
# saveRDS(outputs_tissue,'outputs_tissue.rds')
outputs_tissue = readRDS('outputs_tissue.rds', refhook = NULL)


#=================================dev
# outputs_tissue = readRDS('outputs_tissue.rds', refhook = NULL)
# unique.tissues = unique(rse.jxn.cytosk@colData$tissue)
# outputs_dev_sign_info = list()
# for (tissue in unique.tissues){
#   outputs_dev_sign_info[[tissue]] = getJxnSignInfo(tools.outputs.list=outputs_tissue[[tissue]],
#                                                    logfc_threshold=logfc_threshold,
#                                                    dpsi_threshold=dpsi_threshold,
#                                                    abund_change_threshold=abund_change_threshold,
#                                                    fdr_threshold=fdr_threshold,
#                                                    add_external_data=FALSE, file='')
# }
# saveRDS(outputs_dev_sign_info,'outputs_dev_sign_info.rds')
outputs_dev_sign_info = readRDS('outputs_dev_sign_info.rds', refhook = NULL)
# plotResultsRepot(outputs_dev_sign_info, thresholds = list(logfc_threshold=logfc_threshold,
#                                                            dpsi_threshold=dpsi_threshold,
#                                                            abund_change_threshold=abund_change_threshold,
#                                                            fdr_threshold=fdr_threshold))



#================================= tumor
# -- reading files
outputs_tissue = readRDS('outputs_tissue.rds', refhook = NULL)
# unique.tissues = unique(rse.jxn.cytosk@colData$tissue)
# outputs_gtex2tum = list()
# outputs_norm2tum = list()
# for (tissue in unique.tissues){
#   outputs_gtex2tum[[tissue]] = getJxnSignInfo(outputs_tissue[[tissue]],
#                                               logfc_threshold=logfc_threshold,
#                                               dpsi_threshold=dpsi_threshold,
#                                               abund_change_threshold=abund_change_threshold,
#                                               fdr_threshold=fdr_threshold,
#                                               add_external_data=TRUE, file='gtex2tum')
#   outputs_norm2tum[[tissue]] = getJxnSignInfo(outputs_tissue[[tissue]],
#                                               logfc_threshold=logfc_threshold,
#                                               dpsi_threshold=dpsi_threshold,
#                                               abund_change_threshold=abund_change_threshold,
#                                               fdr_threshold=fdr_threshold,
#                                               add_external_data=TRUE, file='norm2tum')
# }
# saveRDS(outputs_gtex2tum,'dev_vs_gtex2tum_tools.rds')
# saveRDS(outputs_norm2tum,'dev_vs_norm2tum_tools.rds')

outputs_gtex2tum = readRDS('dev_vs_gtex2tum_tools.rds', refhook = NULL)
outputs_norm2tum = readRDS('dev_vs_norm2tum_tools.rds', refhook = NULL)


fisher_results_tissues_list = list()
for (tissue in unique.tissues){
  fisher_results_tissues_list[[tissue]]$dev = getFisher(outputs_dev_sign_info[[tissue]],
                                                        one_to_all=FALSE, 
                                                        ref_col='')
  
  fisher_results_tissues_list[[tissue]]$gtex2tum = getFisher(outputs_gtex2tum[[tissue]],
                                                      one_to_all=TRUE, 
                                                      ref_col='sajr.norm.tumor')
  names(fisher_results_tissues_list[[tissue]]$gtex2tum) = gsub("sajr.norm.tumor", 
                                                               "sajr.gtex2tum", 
                                                               names(fisher_results_tissues_list[[tissue]]$gtex2tum))
  
  fisher_results_tissues_list[[tissue]]$norm2tum = getFisher(outputs_norm2tum[[tissue]],
                                                             one_to_all=TRUE, 
                                                             ref_col='sajr.norm.tumor')
  names(fisher_results_tissues_list[[tissue]]$norm2tum) = gsub("sajr.norm.tumor", 
                                                               "sajr.norm2tum", 
                                                               names(fisher_results_tissues_list[[tissue]]$norm2tum))
  
  fisher_results_tissues_list[[tissue]] = 
    Reduce(append, fisher_results_tissues_list[[tissue]])
  
}

plotFisherResults(fisher_results_tissues_list)

# p.adjust(p_values, method = "bonferroni")

# 
# plotResultsRepot(outputs_gtex2tum, tumor=TRUE, file='gtex2tum',
#                  thresholds = list(logfc_threshold=logfc_threshold,
#                                    dpsi_threshold=dpsi_threshold,
#                                    abund_change_threshold=abund_change_threshold,
#                                    fdr_threshold=fdr_threshold))
# 
# plotResultsRepot(outputs_norm2tum, tumor=TRUE, file='norm2tum',
#                  thresholds = list(logfc_threshold=logfc_threshold,
#                                    dpsi_threshold=dpsi_threshold,
#                                    abund_change_threshold=abund_change_threshold,
#                                    fdr_threshold=fdr_threshold))





