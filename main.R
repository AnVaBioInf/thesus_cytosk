source("downloadRseData.R")
source("findSignificantEvents.R")
source('plotToolsResults.R')
source('plotGenesExpression.R')
source('runTools.R')
library(SummarizedExperiment)

logfc_threshold=1.5
dpsi_threshold=0.1
abund_change_threshold=1
fdr_threshold=0.05


thresholds = list(logfc_threshold=logfc_threshold,
                  dpsi_threshold=dpsi_threshold,
                  abund_change_threshold=abund_change_threshold,
                  fdr_threshold=fdr_threshold)

#===============================================================================
# #==================== DOWNLOAD AND FILTER rse ================================
# ==============================================================================
# development
# prepareRse()
rse.gene.cytosk = readRDS('rds/rse.gene.cytosk.rds', refhook = NULL)
rse.jxn.cytosk = readRDS('rds/rse.jxn.cytosk.rds', refhook = NULL)

# # Breast normal tissue
# gtex.breast = prepareGeneRseAssay('BREAST', 'gene')
# gtex.breast$tissue = 'Breast_normal'
# saveRDS(gtex.breast,'rds/gtex.breast.rds')
gtex.breast = readRDS('rds/gtex.breast.rds')

# # BRCA
# prepareRse(project.id = 'BRCA', condition_col_name = "tissue",
#            file_name_jxn_rse = 'rse.jxn.brca.cytosk.rds', file_name_gene_rse = 'rse.gene.brca.cytosk.rds',
#            tumor=TRUE)
# rse.gene.brca.cytosk = readRDS('rds/rse.gene.brca.cytosk.rds')
# rse.jxn.brca.cytosk = readRDS('rds/rse.jxn.brca.cytosk.rds')

rse.gene.brca.cytosk = readRDS('rds/rse.gene.brca.cytosk.rds')
rse.jxn.brca.cytosk = readRDS('rds/rse.jxn.brca.cytosk.rds')

# # ==============================================================================
# # ============================= SAMPLES HEATMAP ================================
# # ==============================================================================
# plotHeatmapSamplesTissueAge(rse.gene.cytosk)
# table(rse.gene.brca.cytosk@colData$tissue)
# 
# # ==============================================================================
# # # ========================= GENE EXPRESSION PLOT =============================
# # ==============================================================================
# 
# # ========================== gene expression vs age ============================
# png(paste0('plots/cytosk_gene_expression_vs_age.png'), width = 45, height = 30, units = "cm", res = 700)
# plotScatterplotExpression(rse.gene.cytosk)
# dev.off()
# 
# # ============================ barplots ========================================
# # Extract assay data
# merged_rse = mergeRse(list(rse.gene.cytosk, rse.gene.brca.cytosk, gtex.breast))
# setParams(merged_rse)
# 
# png(paste0('plots/cytosk_gene_expression.png'), width = 45, height = 30, units = "cm", res = 700)
# plotBoxplotExpression(merged_rse)
# dev.off()

#==============================================================================
#=============================== RUNNING TOOLS ================================
#==============================================================================

#================================= development ================================
# # -- reading files
# unique.tissues = unique(rse.jxn.cytosk@colData$tissue)
# rse.jxn.cytosk@colData$age_group
# outputs_tissue = list()
# for (tissue in unique.tissues){
#   print(tissue)
#   if (tissue=='Testis'){
#     age_group= c('fetus', 'infant')
#     tissue_age = paste(tissue, paste(age_group, collapse='_'), sep = "_")
#     print(tissue_age)
#     outputs_tissue[[tissue_age]] = runTools(rse.jxn.cytosk, tissue, age_group = age_group,
#                                             reference_condition='fetus')
#     
#     age_group= c('infant', 'adult')
#     tissue_age = paste(tissue, paste(age_group, collapse='_'), sep = "_")
#     outputs_tissue[[tissue_age]] = runTools(rse.jxn.cytosk, tissue, age_group = age_group,
#                                             reference_condition='infant')
#   } else{
#     outputs_tissue[[tissue]] = runTools(rse.jxn.cytosk, tissue)
#   }
# }
# saveRDS(outputs_tissue,'rds/outputs_tissue.rds')
outputs_tissue = readRDS('rds/outputs_tissue.rds', refhook = NULL)
# 
# 

# #================================= tumor =====================================
outputs_gtex2tum = downloadExternalOutputs(file='gtex2tum')
outputs_norm2tum = downloadExternalOutputs(file='norm2tum')

#===============================================================================
#========================= FINDING SIGNIFICANT EVENTS ==========================
#===============================================================================


#================================= development =================================
outputs_tissue = readRDS('rds/outputs_tissue.rds', refhook = NULL)
outputs_dev_sign_info = list()
for (tissue in names(outputs_tissue)){
  print(tissue)
  outputs_dev_sign_info[[tissue]] = getJxnSignInfo(tools.outputs.list=outputs_tissue[[tissue]],
                                                   logfc_threshold=logfc_threshold,
                                                   dpsi_threshold=dpsi_threshold,
                                                   abund_change_threshold=abund_change_threshold,
                                                   fdr_threshold=fdr_threshold,
                                                   add_external_data=FALSE, file='')
}
saveRDS(outputs_dev_sign_info,'rds/outputs_dev_sign_info.rds')

outputs_dev_sign_info = readRDS('rds/outputs_dev_sign_info.rds', refhook = NULL)
plotResultsRepot(outputs_dev_sign_info, thresholds = thresholds)



# #================================= tumor =====================================
#-- reading files
outputs_tissue = readRDS('rds/outputs_tissue.rds', refhook = NULL)
outputs_gtex2tum = list()
outputs_norm2tum = list()

for (tissue in names(outputs_tissue)){
  print(tissue)
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
saveRDS(outputs_gtex2tum,'rds/dev_vs_gtex2tum_tools.rds')
saveRDS(outputs_norm2tum,'rds/dev_vs_norm2tum_tools.rds')

outputs_gtex2tum = readRDS('rds/dev_vs_gtex2tum_tools.rds', refhook = NULL)
outputs_norm2tum = readRDS('rds/dev_vs_norm2tum_tools.rds', refhook = NULL)


# ============================================================================
# =============================== METRICS PLOT ===============================
# ============================================================================

plotResultsRepot(outputs_gtex2tum, tumor=TRUE, file='gtex2tum', thresholds = thresholds)
plotResultsRepot(outputs_norm2tum, tumor=TRUE, file='norm2tum', thresholds = thresholds)


# ==============================================================================
# ============================= STATISTICAL TEST ===============================
# ==============================================================================
fisher_results_tissues_list = list()
for (tissue in names(outputs_tissue)){
  print(tissue)
  fisher.df = data.frame(tool_pair=character(), odds_ratio = numeric(), p_val = numeric())
  print(fisher.df)
  fisher.df = getFisher(fisher.df, outputs_dev_sign_info[[tissue]], one_to_all=FALSE,  ref_col='')
  print(fisher.df)
  fisher.df = getFisher(fisher.df, outputs_gtex2tum[[tissue]], one_to_all=TRUE, ref_col='sajr.norm.tumor')
  fisher.df$tool_pair = gsub("sajr.norm.tumor",  "sajr.gtex2tum",  fisher.df$tool_pair)
  fisher.df = getFisher(fisher.df, outputs_norm2tum[[tissue]], one_to_all=TRUE, ref_col='sajr.norm.tumor')
  fisher.df$tool_pair = gsub("sajr.norm.tumor",  "sajr.norm2tum",  fisher.df$tool_pair)
  fisher.df$q_val = p.adjust(fisher.df$p_val, method = "BH") 
  fisher_results_tissues_list[[tissue]] = fisher.df
}
plotFisherResults(fisher_results_tissues_list, thresholds=thresholds, log=TRUE)


# ==============================================================================
# ============================= COVARIDGES POLTS ===============================
# ==============================================================================
dev.off()
common_sign_jxns_outputs_gtex2tum = findCommonJxns(outputs_gtex2tum)
common_sign_jxns_outputs_norm2tum = findCommonJxns(outputs_norm2tum)

common_sign_jxns = merge(common_sign_jxns_outputs_gtex2tum,common_sign_jxns_outputs_norm2tum,
                         by = c('junction_id_sajr', 'tissue', 'junction_id', 'gene_name', 'gene_id',
                                'dPSI_sajr', 'FDR_sajr', 'logFC_dje', 'FDR_dje',
                                'abund_change_diego', 'FDR_diego'), all=TRUE)
common_sign_jxns = 
  common_sign_jxns[order(common_sign_jxns$gene_name, 
                         common_sign_jxns$junction_id,
                         common_sign_jxns$tissue), ]

write.table(common_sign_jxns, file = "jxns_common_dev_cancer.csv",
            sep = ",",  row.names = TRUE, col.names = TRUE)

# runJunxtionPlot(common_sign_jxns)

#=========================significant only in one tool
all.jxns.info = outputs_dev_sign_info[['Testis']]$all.jxns.info
unique.to.tool.ids = outputs_dev_sign_info[['Testis']]$sign.jxns.info.list$intersections$unique.to.tool

all.jxns.df = lapply(unique.to.tool.ids, function(tool) 
  all.jxns.info[all.jxns.info$junction_id_sajr %in% tool,])


all.jxns.df = lapply(all.jxns.df, function(x) {
  # Sort by absolute value of dPSI_sajr (descending)
  x <- x[order(abs(x$dPSI_sajr), decreasing = TRUE), ] 
  # Remove duplicates based on junction_id
  x <- x[!duplicated(x$junction_id), ] 
  return(x) # Return the modified data frame
})

all.jxns.df$dje = all.jxns.df$dje[abs(all.jxns.df$dje$logFC_dje)>=2,]
all.jxns.df$sajr = all.jxns.df$sajr[abs(all.jxns.df$sajr$dPSI_sajr)>=0.2,]

# for every gene
tissue = 'Testis'
par(mfrow = c(4,2), bty='n')
for (tool in all.jxns.df){
  for (jxn in unique(tool$junction_id)){
    print('gere')
    # setting number of plot rows to number of tissues where selected gene junctions were significant, but no more than 3
    # selecting significant junctions for the GENE in both, development and cancer
    sign.jxn.df = tool[tool$junction_id==jxn,]
    gene = unique(sign.jxn.df$gene_name)
    # for every tissue where any of selected junctions are significant
    covs.summed.gene = prepareCovs(gene,rse.jxn.cytosk,tissue)
    
    fetus.covs.summed.gene =  covs.summed.gene[['fetus.covs.summed.gene']]
    adult.covs.summed.gene = covs.summed.gene[['adult.covs.summed.gene']]
    
    gene.region.coords = strsplit(jxn,'[:-]')
    gene.region.coords = as.integer(c(gene.region.coords[[1]][2], gene.region.coords[[1]][3]))
    gene.region.coords = c(gene.region.coords[1], gene.region.coords[2])
    
    fetus.covs.summed.gene = set.colors(jxn, fetus.covs.summed.gene)
    adult.covs.summed.gene = set.colors(jxn, adult.covs.summed.gene)
    
    sign.jxn.tissue = sign.jxn.df[sign.jxn.df$tissue==tissue,]
    
    text = paste(paste(tissue,gene, jxn), " \n(",
                 'dPSI.sajr=',round(sign.jxn.df$dPSI_sajr, digits=2), 
                 ', FDR.sajr=', round(sign.jxn.df$FDR_sajr, digits=2),
                 ', logFC.dje=', round(sign.jxn.df$logFC_dje, digits=2),
                 ', FDR.dje=', round(sign.jxn.df$FDR_dje, digits=2),
                 ', sign.diego=', round(sign.jxn.df$abund_change_diego, digits=2),
                 ', FDR.diego=', round(sign.jxn.df$FDR_diego, digits=2), ")")
    
    plotReadCov(fetus.covs.summed.gene,
                junc.col = fetus.covs.summed.gene$cols,
                xlim=gene.region.coords,
                plot.junc.only.within = F,
                min.junc.cov.f = 0.05,
                sub='Before birth'
    )
    mtext(text, side=3, line=0.5, cex=0.7, adj=0) 
    
    plotReadCov(adult.covs.summed.gene,
                junc.col = fetus.covs.summed.gene$cols,
                xlim=gene.region.coords,
                plot.junc.only.within = F,
                min.junc.cov.f = 0.05,
                sub='After birth')
  }
}

