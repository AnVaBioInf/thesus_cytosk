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




#==============================================COVARIDGES========================
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
# for every gene
for (jxn in unique(common_sign_jxns$junction_id)){
  # setting number of plot rows to number of tissues where selected gene junctions were significant, but no more than 3
  # selecting significant junctions for the GENE in both, development and cancer
  sign.jxn.df = common_sign_jxns[common_sign_jxns$junction_id==jxn,]
  gene = unique(sign.jxn.df$gene_name)
  par(mfrow = c(min(length(unique(sign.jxn.df$tissue)),3),2), bty='n')
  # for every tissue where any of selected junctions are significant
  for (tissue in (unique(sign.jxn.df$tissue))){ 
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
      'dPSI.sajr=',round(sign.jxn.tissue$dPSI_sajr, digits=2), 
      'FDR.sajr=', round(sign.jxn.tissue$FDR_sajr, digits=2),
      'logFC.dje=', round(sign.jxn.tissue$logFC_dje, digits=2),
      'FDR.dje=', round(sign.jxn.tissue$FDR_dje, digits=2),
      'sign.diego=', round(sign.jxn.tissue$abund_change_diego, digits=2),
      'FDR.diego=', round(sign.jxn.tissue$FDR_diego, digits=2), "\n",
      'dPSI_gtex2tum',round(sign.jxn.tissue$dPSI_gtex2tum, digits=2),
      'FDR_gtex2tum',round(sign.jxn.tissue$FDR_gtex2tum, digits=2),
      'dPSI_norm2tum',round(sign.jxn.tissue$dPSI_norm2tum, digits=2),
      'FDR_norm2tum',round(sign.jxn.tissue$FDR_norm2tum, digits=2), ")")
    
    plotReadCov(fetus.covs.summed.gene,
                junc.col = fetus.covs.summed.gene$cols,
                xlim=gene.region.coords,
                plot.junc.only.within = F,
                min.junc.cov.f = 0.05,
                sub='Before birth'
    )
    mtext(text, side=3, line=1.5, cex=0.7, adj=0) 
    
    plotReadCov(adult.covs.summed.gene,
                junc.col = fetus.covs.summed.gene$cols,
                xlim=gene.region.coords,
                plot.junc.only.within = F,
                min.junc.cov.f = 0.05,
                sub='After birth')

  }
}



