downloadExternalOutputs = function(path = 'sign_norm_vs_all_tumor/',
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

mergeOutputs = function(output.list, add_external_data, file){
  diego.output = output.list$diego.output 
  dje.output = output.list$dje.output 
  sajr.output = output.list$sajr.output

  diego.output = diego.output[, c('junction', 'abundance_change', 'q_val', 'geneID', 'geneName')]
  dje.output = dje.output$dje.out[, c('junctionID', 'logFC', 'FDR')]
  sajr.output = sajr.output[,c('dPSI', 'FDR.sajr')]
  
  colnames(diego.output) = c("junction_id", "abund_change_diego", "FDR_diego", 'gene_id', 'gene_name')
  colnames(dje.output) = c("junction_id", "logFC_dje", "FDR_dje")
  colnames(sajr.output) = c("dPSI_sajr", "FDR_sajr")
  
  # consideting that diego substracts from base-test
  #diego.output$abund_change_diego = (diego.output$abund_change_diego)*(-1)
  
  dje.output$junction_id = sub(":.$", "", dje.output$junction_id)
  dje.output$junction_id = sub("([0-9]{3,})(:)([0-9]{3,})", "\\1-\\3", dje.output$junction_id)
  sajr.output$junction_id = sub(':[-+*]:.*','', rownames(sajr.output))
  sajr.output$junction_id_sajr = rownames(sajr.output)
  
  output.merged.df = merge(dje.output, sajr.output, by = "junction_id", all = TRUE)
  output.merged.df = merge(output.merged.df, diego.output, by = "junction_id", all = TRUE)
  output.merged.df[,c('junction_id_sajr', 'junction_id', 'gene_id', 'gene_name',
                      'logFC_dje', 'FDR_dje', 'dPSI_sajr', 'FDR_sajr',
                      'abund_change_diego', 'FDR_diego')]
  
  if (add_external_data){
    tum = downloadExternalOutputs(file.name=file)
    output.merged.df = merge(output.merged.df, tum, by = "junction_id_sajr", all = TRUE)
  }
  
  print('Finished merging outputs')
  return(output.merged.df) 
}

compareOutputs = function(jxn.ids.list) {
  jxn.ids.list = lapply(jxn.ids.list, function(x) x[!is.na(x)])
  # 1. Find elements present in all vectors
  all.tools = Reduce(intersect, jxn.ids.list)
  
  # 2. Find elements unique to each vector
  unique.to.tool = lapply(seq_along(jxn.ids.list), function(i) {
    unique = setdiff(jxn.ids.list[[i]], unname(unlist(jxn.ids.list[-i])))
    return(unique)
  })
  names(unique.to.tool) = names(jxn.ids.list)
  
  # 3. Find elements unique to subsets of jxn.ids.list (excluding those in all)
  only.pair.tools = list()
  tool.pairs = combn(seq_along(jxn.ids.list), 2, simplify = FALSE)  # Get all pairs of vector indices
  for (pair in tool.pairs) {
    i = pair[[1]]
    j = pair[[2]]
    common.pair = intersect(jxn.ids.list[[i]], jxn.ids.list[[j]])
    common.pair = setdiff(common.pair, all.tools)
    pair.name = paste(names(jxn.ids.list)[c(i, j)], collapse = "&")
    only.pair.tools[[pair.name]] = common.pair
  }
  all.tools = list(all.tools)
  names(all.tools) = paste(names(jxn.ids.list), collapse = "&")
  
  # ids2df = function(ids) jxn.sign.df[jxn.sign.df$junction_id_sajr %in% ids, ]
  
  # sign.jxns.info = 
  #   list(list(all.tools = all.tools),
  #        only.pair.tools = only.pair.tools,
  #        unique.to.tool = unique.to.tool)
  
  sign.jxns.info =list(
    unique.to.tool = unique.to.tool,
    only.pair.tools = only.pair.tools,
    all.tools = all.tools)
  
  # jxn.ids.list = lapply(jxn.ids.list, ids2df)
  
  return(list(
    all.single.tool =  jxn.ids.list,
    intersections = sign.jxns.info ))
}

# сделать ее универсальной и все сохранить в список
# то есть на вход один параметр за раз
findSignificantJxnsIds = function(jxns.significance.df, logfc_threshold, fdr_threshold, dpsi_threshold, abund_change_threshold,
                                  add_external_data, file.name){
  # 1. Define filtering conditions for each tool
  diego.sign.tf = abs(jxns.significance.df$abund_change_diego) >= abund_change_threshold & jxns.significance.df$FDR_diego <= fdr_threshold
  dje.sign.tf = abs(jxns.significance.df$logFC_dje) >= logfc_threshold & jxns.significance.df$FDR_dje <= fdr_threshold
  sajr.sign.tf = abs(jxns.significance.df$dPSI_sajr) >= dpsi_threshold & jxns.significance.df$FDR_sajr <= fdr_threshold  
  
  # 2. Extract significant junction IDs for each tool
  all.sign.jxns.diego = jxns.significance.df[diego.sign.tf, 'junction_id_sajr']
  all.sign.jxns.dje = jxns.significance.df[dje.sign.tf, 'junction_id_sajr']
  all.sign.jxns.sajr = jxns.significance.df[sajr.sign.tf, 'junction_id_sajr']
  
  all.sign.jxn.ids.tool.list = list(diego = all.sign.jxns.diego, 
                                    dje = all.sign.jxns.dje, 
                                    sajr = all.sign.jxns.sajr)
  
  sign.jxns.info = compareOutputs(all.sign.jxn.ids.tool.list)
  if (add_external_data){
    dpsi_file = paste0('dPSI_sajr_', file.name)
    fdr_file = paste0('FDR_sajr_', file.name)
    
    tum.sign.tf = abs(jxns.significance.df[, dpsi_file]) >= dpsi_threshold & jxns.significance.df[, fdr_file] <= fdr_threshold  
    all.sign.jxns.tum = jxns.significance.df[tum.sign.tf, 'junction_id_sajr']
    all.sign.jxns.tum = all.sign.jxns.tum[!is.na(all.sign.jxns.tum)]
    sign.jxns.info.d = lapply(sign.jxns.info$intersections, function(sub) 
      lapply(sub, function(x) intersect(x, all.sign.jxns.tum)))
    sign.jxns.info$intersections = sign.jxns.info.d
    tumor_name = paste0('sajr_', file.name)
    sign.jxns.info$intersections$unique.to.tool[[tumor_name]] = 
      setdiff(all.sign.jxns.tum, unname(unlist(sign.jxns.info.d)))
    sign.jxns.info$all.single.tool[[tumor_name]] = all.sign.jxns.tum
  } 
  print('Finished intersecting junction ids between tools')
  return(sign.jxns.info)
  
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


#===============================================================================
#============================ FISHER EXACT TEST ================================
#===============================================================================
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
  all.sign.jxns.ids.list.tool = outputs_tissue$sign.jxns.info.list$all.single.tool
  
  # from df overall number of input junctions
  all.jxns.df = outputs_tissue$all.jxns.info
  
  all.tool.pairs.comb = makePairs(names(all.sign.jxns.ids.list.tool),
                                  one_to_all=one_to_all, 
                                  ref_col=ref_col)
  for (tool.pair in all.tool.pairs.comb){
    tool1.name = tool.pair[1]
    tool2.name = tool.pair[2]
    all.sign.jxns.tool.pair.ids = all.sign.jxns.ids.list.tool[tool.pair]
    sign.jxns.intersect.ids = compareOutputs(all.sign.jxns.tool.pair.ids)
    
    
    sign.both.tools = length(sign.jxns.intersect.ids$intersections$all.tools[[1]])
    sign.tool.1 = length(sign.jxns.intersect.ids$intersections$unique.to.tool[[tool1.name]])
    sign.tool.2 = length(sign.jxns.intersect.ids$intersections$unique.to.tool[[tool2.name]])
    not.sign = nrow(all.jxns.df)-sign.both.tools-sign.tool.1-sign.tool.2
    
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

#getFisher(outputs_dev_sign_info[['Brain']])

#===================================
#===================================
#================================================================================

addLeftmostRightmostCoord = function(common_sign_jxns_dev_tum, rse.jxn){
  sites = makeSites(rse.jxn@rowRanges)
  lefmost_rightmost_info_df = 
    sites$jxns.same.start.or.end.info[,c('gene_id', 'gene_name', 'leftmost', 'rightmost')]
  lefmost_rightmost_info_df$junction_id_sajr = rownames(lefmost_rightmost_info_df)
  
  common_sign_jxns_dev_tum = 
    merge(common_sign_jxns_dev_tum, lefmost_rightmost_info_df, 
          by = c("junction_id_sajr", 'gene_id', 'gene_name'), all.x = TRUE)
  
  common_sign_jxns_dev_tum$lr_most = paste(common_sign_jxns_dev_tum$leftmost, 
                                           common_sign_jxns_dev_tum$rightmost, 
                                               sep = "-")
  common_sign_jxns_dev_tum
}



findCommonJxns = function(outputs_tum, rse.jxn){
  common_sign_jxns_dev_tum = list()
  for (tissue in names(outputs_tum)){
    intersections = outputs_tum[[tissue]]$sign.jxns.info.list$intersections
    all.jxns.info = outputs_tum[[tissue]]$all.jxns.info
    
    intersections = Reduce(append, intersections)
    intersections = intersections[ !grepl("tum", names(intersections)) ]
    
    intersections = lapply(intersections, function(jxns.found.by.tool) 
      all.jxns.info[all.jxns.info$junction_id_sajr %in% jxns.found.by.tool,])
    intersections = do.call(rbind, intersections)
    
    intersections = intersections[order(abs(intersections$dPSI_sajr), decreasing = TRUE), ]
    intersections = intersections[!duplicated(intersections$junction_id), ]
    if (nrow(intersections) > 0) {
      intersections$tissue = tissue
    }
    common_sign_jxns_dev_tum[[tissue]] = intersections
    
  }
  common_sign_jxns_dev_tum = do.call(rbind, common_sign_jxns_dev_tum)
  
  common_sign_jxns_dev_tum = addLeftmostRightmostCoord(common_sign_jxns_dev_tum, rse.jxn)
  
  common_sign_jxns_dev_tum = 
    common_sign_jxns_dev_tum[order(
      common_sign_jxns_dev_tum$gene_name,
      common_sign_jxns_dev_tum$tissue,
      common_sign_jxns_dev_tum$lr_most, 
      decreasing = TRUE), ]

  common_sign_jxns_dev_tum
}












