source("downloadRseData.R")
source("findSignificantEvents.R")

# prepareRse()

# -- reading files
rse.gene.cytosk = readRDS('rse.gene.cytosk.rds', refhook = NULL)
rse.jxn.cytosk = readRDS('rse.jxn.cytosk.rds', refhook = NULL)

unique.tissues = unique(rse.jxn.cytosk@colData$tissue)[1:2]

outputs_tissue = list()
for (tissue in unique.tissues){
  outputs_tissue[[tissue]] = getJxnSignInfo(rse.jxn.cytosk, tissue)
}