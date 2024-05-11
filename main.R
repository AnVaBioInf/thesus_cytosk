source("downloadRseData.R")
source("findSignificantEvents.R")
source('plotToolsResults.R')
# download and filter rse
#prepareRse()

# # -- reading files
# rse.gene.cytosk = readRDS('rse.gene.cytosk.rds', refhook = NULL)
# rse.jxn.cytosk = readRDS('rse.jxn.cytosk.rds', refhook = NULL)
# 
# unique.tissues = unique(rse.jxn.cytosk@colData$tissue)
# outputs_tissue = list()
# for (tissue in unique.tissues){
#   outputs_tissue[[tissue]] = getJxnSignInfo(rse.jxn.cytosk, tissue)
# }
# 
# saveRDS(outputs_tissue,'tool_outputs_all_tissues.rds')
# tool_outputs_all_tissues = readRDS('tool_outputs_all_tissues.rds', refhook = NULL)
# 
# plotResultsRepot(tool_outputs_all_tissues)

# 
# # tumor
# # # -- reading files
rse.gene.cytosk = readRDS('rse.gene.cytosk.rds', refhook = NULL)
rse.jxn.cytosk = readRDS('rse.jxn.cytosk.rds', refhook = NULL)

unique.tissues = unique(rse.jxn.cytosk@colData$tissue)[1]

outputs_gtex2tum = list()
outputs_norm2tum = list()
for (tissue in unique.tissues){
  outputs_gtex2tum[[tissue]] = getJxnSignInfo(rse.jxn.cytosk, tissue, add_external_data=TRUE)
  # outputs_norm2tum[[tissue]] = getJxnSignInfo(rse.jxn.cytosk, tissue, add_external_data=TRUE, file='norm2tum')
}

# outputs_gtex2tum
# 
# #
# 
# saveRDS(outputs_tissue,'tool_outputs_all_tissues_dev_tumor.rds')
# tool_outputs_all_tissues_dev_tumor = readRDS('tool_outputs_all_tissues_dev_tumor.rds', refhook = NULL)
# outputs.dev.cans = tool_outputs_all_tissues_dev_tumor



#plotResultsRepot(outputs_gtex2tum, tumor=TRUE, file='norm2tum')






