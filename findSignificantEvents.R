downloadExternalOutputs = function(path = './',
                                   file.name,
                                   extenction='.csv'){
  tum <- read.csv(paste0(path, file.name, extenction), header = TRUE, sep = ",")
  tum = tum[,c('X', 'dpsi', 'qv')]
  names(tum) = c('junction_id_sajr', 
                 paste0('dPSI_sajr_', file.name), 
                 paste0('FDR_sajr_', file.name))
  print('Finished downloading external data')
  tum
}

mergeOutputs = function(output.list){
  diego.output = output.list$diego.output 
  dje.output = output.list$dje.output 
  sajr.output = output.list$sajr.output
  
  diego.output = diego.output[, c('junction', 'abundance_change', 'q_val', 'geneID', 'geneName')]
  dje.output = dje.output$dje.out[, c('junctionID', 'logFC', 'FDR')]
  sajr.output = sajr.output[,c('dPSI', 'FDR.sajr')]
  
  colnames(diego.output) = c("junction_id", "abund_change_diego", "FDR_diego", 'gene_id', 'gene_name')
  colnames(dje.output) = c("junction_id", "logFC_dje", "FDR_dje")
  colnames(sajr.output) = c("dPSI_sajr", "FDR_sajr")
  
  dje.output$junction_id = sub(":.$", "", dje.output$junction_id)
  dje.output$junction_id = sub("([0-9]{3,})(:)([0-9]{3,})", "\\1-\\3", dje.output$junction_id)
  sajr.output$junction_id = sub(':[-+*]:.*','', rownames(sajr.output))
  sajr.output$junction_id_sajr = rownames(sajr.output)
  
  output.merged.df = merge(dje.output, sajr.output, by = "junction_id", all = TRUE)
  output.merged.df = merge(output.merged.df, diego.output, by = "junction_id", all = TRUE)
  output.merged.df[,c('junction_id_sajr', 'junction_id', 'gene_id', 'gene_name',
                      'logFC_dje', 'FDR_dje', 'dPSI_sajr', 'FDR_sajr',
                      'abund_change_diego', 'FDR_diego')]
  
  print('Finished merging outputs')
  return(output.merged.df) 
}

#=================================================================================
#=================================================================================
#=================================================================================
findSignificantJxnsIds = function(jxns.significance.df, logfc_threshold, fdr_threshold, dpsi_threshold, abund_change_threshold,
                                  suffix=''){
  # 1. Define filtering conditions for each tool
  diego.sign.tf = abs(jxns.significance.df$abund_change_diego) >= abund_change_threshold & jxns.significance.df$FDR_diego <= fdr_threshold
  dje.sign.tf = abs(jxns.significance.df$logFC_dje) >= logfc_threshold & jxns.significance.df$FDR_dje <= fdr_threshold
  sajr.sign.tf = abs(jxns.significance.df$dPSI_sajr) >= dpsi_threshold & jxns.significance.df$FDR_sajr <= fdr_threshold  
  
  # 2. Extract significant junction IDs for each tool
  all.sign.jxns.diego = jxns.significance.df[diego.sign.tf, 'junction_id_sajr']
  all.sign.jxns.dje = jxns.significance.df[dje.sign.tf, 'junction_id_sajr']
  all.sign.jxns.sajr = jxns.significance.df[sajr.sign.tf, 'junction_id_sajr']
  
  all.sign.jxn.ids.tool.list = list(diego = all.sign.jxns.diego[!is.na(all.sign.jxns.diego)], 
                                    dje = all.sign.jxns.dje[!is.na(all.sign.jxns.dje)], 
                                    sajr = all.sign.jxns.sajr[!is.na(all.sign.jxns.sajr)])
  
  # Set the names of the list elements
  names(all.sign.jxn.ids.tool.list) <- c(paste0('diego', suffix), 
                                         paste0('dje', suffix), 
                                         paste0('sajr', suffix))
  
  print('Significant junction in tools found')
  return(all.sign.jxn.ids.tool.list)
}

compareOutputs =  function(list1, list2) {
  all_commons <- list()
  combined_list <- c(list1, list2)
  combined_names <- names(combined_list)
  num_vectors <- length(combined_list)
  
  # Level 0: Common to all vectors
  all_commons[["all"]] <- Reduce(intersect, combined_list)
  
  # Levels 1 to num_vectors - 1: Common to all except 1, 2, ... vectors
  for (exclusion_level in 1:(num_vectors - 1)) {
    excluded_combos <- combn(1:num_vectors, exclusion_level, simplify = FALSE)
    
    for (excluded_indices in excluded_combos) {
      included_indices <- setdiff(1:num_vectors, excluded_indices)
      selected_vectors <- combined_list[included_indices]
      selected_names <- combined_names[included_indices]
      intersection <- Reduce(intersect, selected_vectors)
      
      # Exclude elements found at higher levels
      intersection <- setdiff(intersection, unlist(all_commons))
      
      if (length(intersection) > 0) {
        # Name using included vectors
        key_name <- paste(selected_names, collapse = "_") 
        all_commons[[key_name]] <- intersection
      } else { # If intersection is empty
        key_name <- paste(selected_names, collapse = "_") 
        all_commons[[key_name]] <- character(0) # Assign character(0)
      }
    }
  }
  return(all_commons)
}


getJxnSignInfo = function(tools.outputs.list, 
                          logfc_threshold, dpsi_threshold, abund_change_threshold, fdr_threshold,
                          add_external_data=FALSE, file=''){
  all.jxns.info.df = mergeOutputs(tools.outputs.list, add_external_data, file)
  sign.jxns.info.list = findSignificantJxnsIds(all.jxns.info.df, logfc_threshold, 
                                               fdr_threshold, dpsi_threshold, abund_change_threshold, 
                                               add_external_data, file)
  print('Finished running')
  list(all.jxns.info = all.jxns.info.df, sign.jxns.info.list = sign.jxns.info.list)
}





# outputs_tissue = runTools(rse.jxn.cytosk, 'Brain')
# getJxnSignInfo(outputs_tissue)

makePairs = function(names, ref_col='', one_to_all){
  if (one_to_all){
    # ref.col.name <- grep(ref_col, colnames(df), value = TRUE)
    # Get the names of all other columns
    other.cols <- setdiff(names, ref_col)
    # Create a list of combinations
    all.pairs.comb = lapply(other.cols, function(col_name) c(ref_col, col_name))
  } else{
    all.pairs.comb = combn(names, 2, simplify = FALSE)
  }
  all.pairs.comb
}

getFisher = function(fisher.df, outputs_tissue, one_to_all=FALSE, ref_col=''){
  all.sign.jxns.tool = outputs_tissue$sign.jxns.info.list$all.single.tool
  all.jxns = outputs_tissue$all.jxns.info
  all.tool.pairs.comb = makePairs(names(all.sign.jxns.tool),
                                  one_to_all=one_to_all, 
                                  ref_col=ref_col)
  
  for (tool.pair in all.tool.pairs.comb){
    tool1.name = tool.pair[1]
    tool2.name = tool.pair[2]
    all.sign.jxns.tool.pair.ids = all.sign.jxns.tool[tool.pair]
    sign.jxns.intersect.ids = compareOutputs(all.sign.jxns.tool.pair.ids)
    
    
    sign.both.tools = length(sign.jxns.intersect.ids$intersections$all.tools[[1]])
    sign.tool.1 = length(sign.jxns.intersect.ids$intersections$unique.to.tool[[tool1.name]])
    sign.tool.2 = length(sign.jxns.intersect.ids$intersections$unique.to.tool[[tool2.name]])
    not.sign = nrow(all.jxns)-sign.both.tools-sign.tool.1-sign.tool.2
    
    contingency.table = matrix(c(sign.both.tools, sign.tool.1, sign.tool.2, not.sign), ncol = 2, byrow = TRUE)
    rownames(contingency.table) <- c(paste0(toupper(tool1.name)," Significant"), 
                                     paste0(toupper(tool1.name)," Not Significant"))
    colnames(contingency.table) <- c(paste0(toupper(tool2.name)," Significant"), 
                                     paste0(toupper(tool2.name)," Not Significant"))
    
    result = fisher.test(contingency.table)
    new_row = data.frame(tool_pair = paste(tool.pair, collapse = " & "), 
                         odds_ratio = result$estimate, 
                         p_val = result$p.value)
    fisher.df = rbind(fisher.df, new_row)
  }
  rownames(fisher.df) = c(1:nrow(fisher.df))
  fisher.df
}

findCommonJxns = function(outputs_tum){
  unique.tissue = names(outputs_tum)
  common_sign_jxns_outputs_tum = list()
  for (tissue in unique.tissue){
    intersections = outputs_tum[[tissue]]$sign.jxns.info.list$intersections
    intersections = Reduce(append, intersections)
    intersections = intersections[names(intersections) != 'sajr_norm_tumor']
    all.jxns.info = outputs_tum[[tissue]]$all.jxns.info
    intersections = lapply(intersections, function(tool) 
      all.jxns.info[all.jxns.info$junction_id_sajr %in% tool,])
    intersections = do.call(rbind, intersections)
    intersections = intersections[order(abs(intersections$dPSI_sajr), decreasing = TRUE), ]
    intersections = intersections[!duplicated(intersections$junction_id), ]
    if (nrow(intersections) > 0) {
      intersections$tissue = tissue
    }
    common_sign_jxns_outputs_tum[[tissue]]=intersections
  }
  common_sign_jxns_outputs_tum = do.call(rbind, common_sign_jxns_outputs_tum)
  common_sign_jxns_outputs_tum
}

#getFisher(outputs_dev_sign_info[['Brain']])



#=================================
#######################-----coverage
#=================================

bigWig2Cov = function(bw){
  bw = as.data.frame(bw) # columns: seqnames, start, end, width, strand, score
  bw = bw[order(bw$start),] # ordering bw df my start coordinate column
  start = bw$start[1] # starting coordinate of the gene
  stop = bw$end[nrow(bw)] # end coordinate of the gene
  
  cov = rep(0,stop-start+1) # vector with 0s of length the gene
  
  for(i in 1:nrow(bw)){ # for every record in bw
    cov[(bw$start[i]:bw$end[i])-start+1] = bw$score[i] # assigning score to every position of gene
    # every position of a range is asigned with the same score
  }
  list(cov=cov,x=start:stop) # coverage, start and end gene coordinates on a chromosome
}

# sid - sample id
# jxn - rse object
# gene.grange - granges of genes of interest
import = function(sample.id, gene.grange){
  folder.path = '/home/an/Manananggal/Input/bigWig/'
  sample.path = paste0(folder.path, sample.id, '.bw')
  # download and subset bw file
  bw = rtracklayer::import.bw(sample.path, which=gene.grange) 
  bw
}

get.counts = function(sample.id, rse.gene.tissue){
  # rse for one sample
  rse.tissue.samples = rse.gene.tissue[, sample.id]
  all.genejxn.info = cbind(as.data.frame(rse.gene.tissue@rowRanges)[,c('start','end','strand')],
                           counts = rse.gene.tissue@assays@data$counts[,sample.id]) # ?? зачем еще раз выбирать sample.id?
  all.genejxn.info
}

getRecountCov = function(sample.id, rse.gene.tissue, gene.grange){
  bw = import(sample.id, gene.grange)
  # coverage on a gene, start:stop 
  cov.list.sample = bigWig2Cov(bw) 
  cov.list.sample$juncs =  get.counts(sample.id, rse.gene.tissue)
  # coverages and counts for the entire gene (all jxns)
  cov.list.sample
}

sumCovs = function(gene.cov.samples.list){
  # launched for each of the conditions (2ce)
  # launched for every tissue where significant junction was found
  # creating a template of merged object
  cov.merged = gene.cov.samples.list[[1]] # read coverages for the first sample
  
  # Extract "val" elements from each sublist
  cov.list <- lapply(gene.cov.samples.list, `[[`, "cov")
  # Sum corresponding elements using Reduce
  cov.merged$cov <- Reduce("+", cov.list)
  
  juncs.list <- lapply(gene.cov.samples.list, `[[`, "juncs")
  counts.list <- lapply(juncs.list, `[[`, "counts")
  cov.merged$juncs$counts <- Reduce("+", counts.list)
  cov.merged
}


filter.data = function(condition.cov.list, xlim, min.junc.cov,plot.junc.only.within, min.junc.cov.f){
  #print(condition.cov.list)
  # choosing only x inside the range of interest
  x.in.range.tf = condition.cov.list$x >= xlim[1]-500 &
    condition.cov.list$x <= xlim[2]+500
  #??? x is from 1 to what?
  
  condition.cov.list$x = condition.cov.list$x[x.in.range.tf]
  condition.cov.list$cov = condition.cov.list$cov[x.in.range.tf]
  
  condition.cov.list$juncs =
    condition.cov.list$juncs[(condition.cov.list$juncs$start == xlim[1] |
                                condition.cov.list$juncs$end == xlim[2]) &
                               condition.cov.list$juncs$counts > min.junc.cov.f,]
  
  condition.cov.list$cov[c(1,length(condition.cov.list$cov))] = 0 # assigning cov 0 to first and last elements of cov
  condition.cov.list
}

prepareCovs = function(gene, rse.gene.cytosk, tissue){
  # gene.grange is needed by rtracklayer to filter bw files
  gene.grange = rse.gene.cytosk@rowRanges[rse.gene.cytosk@rowRanges$gene_name==gene,]
  rse.gene = rse.gene.cytosk[rse.gene.cytosk@rowRanges$gene_name==gene,]
  rse.gene = rse.gene[,rse.gene@colData$tissue==tissue]
  
  rse.gene = rse.gene[,rse.gene@colData$age_group %in% c('fetus','adult')]
  adult.samples.ids = rownames(rse.gene@colData[rse.gene@colData$age_group == 'adult',])
  fetus.samples.ids = rownames(rse.gene@colData[rse.gene@colData$age_group == 'fetus',])
  
  all.samples.ids = c(adult.samples.ids, fetus.samples.ids)
  # -- coverages
  # covearge for each sample, output is a list of lists with read coverages, start:end positions, juncs df
  # gene.cov.all.samples.list - named list for each tissue sample
  gene.cov.all.samples.list = 
    lapply(all.samples.ids, function(samples.id){getRecountCov(samples.id, rse.gene, gene.grange)}) 
  # assigning elements of the list sample ids names
  names(gene.cov.all.samples.list) = all.samples.ids
  
  # --- merging
  # sum coverage in each condition
  fetus.covs.summed.gene = sumCovs(gene.cov.all.samples.list[fetus.samples.ids])
  adult.covs.summed.gene = sumCovs(gene.cov.all.samples.list[adult.samples.ids])
  
  list(fetus.covs.summed.gene=fetus.covs.summed.gene, adult.covs.summed.gene=adult.covs.summed.gene)
} 
