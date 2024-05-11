source("downloadRseData.R")
source("findSignificantEvents.R")

# prepareRse()

# -- reading files
rse.gene.cytosk = readRDS('rse.gene.cytosk.rds', refhook = NULL)
rse.jxn.cytosk = readRDS('rse.jxn.cytosk.rds', refhook = NULL)

# getting jxn significance in development
unique.tissues = unique(rse.jxn.cytosk@colData$tissue)
# # outputs_tissue = list()
# all.jxns.info.df.list = list()
# for (tissue in unique.tissues){
#   tools.outputs.list = runTools(rse.jxn.cytosk, tissue)
#   all.jxns.info.df.list[[tissue]] = mergeOutputs(tools.outputs.list)
# }

#saveRDS(all.jxns.info.df.list,'all_jxns_info_df_list.rds')
all.jxns.info.df.list = readRDS('all_jxns_info_df_list.rds', refhook = NULL)

#   #outputs_tissue[[tissue]] = getJxnSignInfo(all.jxns.info.df)
#outputs_dev_all_tissues = readRDS('tool_outputs_all_tissues.rds', refhook = NULL)
#plotResultsRepot(outputs_dev_all_tissues)


# cancer
gtex2tum <- read.csv("./gtex2tum.csv", header = TRUE, sep = ",")
norm2tum <- read.csv("./norm2tum.csv", header = TRUE, sep = ",")

gtex2tum = gtex2tum[,c('X', 'dpsi', 'qv')]
norm2tum = norm2tum[,c('X', 'dpsi', 'qv')]
names(gtex2tum) = c('junction_id_sajr', 'dPSI_gtex2tum', 'FDR_gtex2tum')
names(norm2tum) = c('junction_id_sajr', 'dPSI_norm2tum', 'FDR_norm2tum')

head(norm2tum)
head(all.jxns.info.df.list[['Brain']])

all.jxns.info.dev.cans.list = list()
for (tissue in unique.tissues){
  all.jxns.info.dev.cans.list[[tissue]] = merge(all.jxns.info.df.list[[tissue]], gtex2tum, 
                                                by = "junction_id_sajr", all = FALSE)
  all.jxns.info.dev.cans.list[[tissue]] = merge(all.jxns.info.dev.cans.list[[tissue]], norm2tum, 
                                                by = "junction_id_sajr", all = FALSE)
  all.jxns.info.dev.cans.list[[tissue]] = all.jxns.info.dev.cans.list[[tissue]][,c('junction_id_sajr', 'junction_id', 'logFC_dje', 'FDR_dje', 'dPSI_sajr', 'FDR_sajr',
                                           'abund_change_diego', 'FDR_diego', 'dPSI_gtex2tum', 'FDR_gtex2tum',
                                           'dPSI_norm2tum', 'FDR_norm2tum', 'gene_id', 'gene_name')]
  
}

# сделаем отдельно для norm2tum и gtex2tum
for (tissue in unique.tissues){
  sign.tf = abs(all.jxns.info.dev.cans.list[[tissue]]$dPSI_norm2tum) >= 0.2 & all.jxns.info.dev.cans.list[[tissue]]$FDR_norm2tum <= 0.05  
  all.jxns.info.dev.cans.list[[tissue]][sign.tf, 'junction_id_sajr']
  
  sign.jxns.info.list = findSignificantJxnsIds(all.jxns.info.dev.cans.list[[tissue]], 
                                               logfc_threshold=2, fdr_threshold=0.05, 
                                               dpsi_threshold=0.2, abund_change_threshold=1)
  sign.jxns.info.list
  
}
sajr.sign.tf = abs(all.jxns.info.dev.cans.list[[tissue]]$dPSI) >= 0.2 & all.jxns.info.dev.cans.list[[tissue]]$FDR_sajr <= 0.05  

all.jxns.info.dev.cans.list

# find significant


# находим общие в sajr... и в опухоли


grep("tum", names(all.jxns.info.dev.cans.list$Brain), value = TRUE)


#==================================================

