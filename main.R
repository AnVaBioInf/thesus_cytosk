source("downloadRseData.R")
source("findSignificantEvents.R")
source('plotToolsResults.R')
# download and filter rse
#prepareRse()

# # -- reading files
# rse.gene.cytosk = readRDS('rse.gene.cytosk.rds', refhook = NULL)
# rse.jxn.cytosk = readRDS('rse.jxn.cytosk.rds', refhook = NULL)
# 
unique.tissues = unique(rse.jxn.cytosk@colData$tissue)
# outputs_tissue = list()
# for (tissue in unique.tissues){
#   outputs_tissue[[tissue]] = getJxnSignInfo(rse.jxn.cytosk, tissue)
# }
# 
# saveRDS(outputs_tissue,'tool_outputs_all_tissues.rds')
tool_outputs_all_tissues = readRDS('tool_outputs_all_tissues.rds', refhook = NULL)

plotResultsRepot(tool_outputs_all_tissues)
