source('downloadRseData.R')
source('findSignificantEvents.R')

# prepareRse()

# -- reading files
rse.gene.cytosk = readRDS('rse.gene.cytosk.rds', refhook = NULL)
rse.jxn.cytosk = readRDS('rse.jxn.cytosk.rds', refhook = NULL)

#getJxnSignInfo(rse.jxn.cytosk, 'Brain')
