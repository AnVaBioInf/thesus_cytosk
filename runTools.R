# the project goal is to search for AS events in development
# DIEGO, DJexpress and SAJR tools are used for that (work with junction counts)
# in this file functions to make input for those tools are gathered
library(reticulate)
library(SAJR)
library(DJExpress)
library(recount3)
library(rtracklayer)

rse2countDf = function(rse){
  counts = as.matrix(assay(rse, "counts"))
  counts = as.data.frame(counts)
  counts
}

# проверить!
filterRse = function(rse, tissue, age_group=c('fetus', 'adult'), mim_row_sum=10, min_variance=0){
  rse = rse[, (rse@colData$tissue %in% tissue) & (rse@colData$age_group %in% age_group)]
  rse = rse[apply(rse@assays@data$counts, 1, sum) >= mim_row_sum &
              apply(rse@assays@data$counts, 1, var) >= min_variance, ]
  rse
}

findConditionIds = function(rse_filtered, condition_col_name, condition_name){
  rownames(rse_filtered@colData[rse_filtered@colData[,condition_col_name] == condition_name,])
}

#=================================DIEGO=======================================
# http://legacy.bioinf.uni-leipzig.de/Software/DIEGO/
# -- a file (table of splice junction supports per sample)
makeAFile = function(rse.filtered, tissue, age_group, path_input){
  junction_table = rse2countDf(rse.filtered)
  sample_ids = colnames(junction_table)
  # location of the splice junction or exon, which must have the string:number-number format
  junction_table$junction = sub('(:[-+*]$)', '' , rownames(junction_table))
  #  the type of splice junction
  junction_table$type = 'N_w'
  junction_table$geneID = rse.filtered@rowRanges$gene_id
  junction_table$geneName = rse.filtered@rowRanges$gene_name
  junction_table = junction_table[,c('junction', 'type', sample_ids, 'geneID', 'geneName')]
  junction_table = junction_table[order(junction_table$geneID),]
  write.table(junction_table,
              paste0(path_input, "/junction_table_", paste(tissue, collapse = "_"), "_",
                     paste0(age_group, collapse = "_"), ".txt"),
              sep = "\t",
              row.names = FALSE, col.names = TRUE, quote = FALSE)
}

# --b file (condition to sample relation in the format: condition tab-delimiter sampleName)
makeBFile = function(rse.filtered, tissue, age_group, path_input,
                     condition_col_name){
  group_table = rse.filtered@colData[,condition_col_name,drop=F]
  group_table$sample_id = rownames(group_table)
  group_table = group_table[order(group_table[,condition_col_name]), ]
  write.table(group_table,
              paste0(path_input, "/group_table_", paste(tissue, collapse = "_"), "_",
                     paste(age_group, collapse = "_"), ".txt"),
              sep = "\t",
              row.names = FALSE, col.names = FALSE, quote = FALSE)
}

# makeDiegoBashFile = function(tissue, path_input, path_output, reference_condition,
#                              min_support, min_samples, FDR_threashold, FC_threashold){
#   # Define the Bash script content
#   diego_bash_file_content = paste0(
#   '#!/bin/bash
#   python /home/an/DIEGO/diego.py \\\
#   -a ', path_input, '/junction_table_', tissue, '.txt \\\
#   -b ', path_input, '/group_table_', tissue, '.txt \\\
#   -x ', reference_condition, ' \\\
#   --minsupp ', min_support,' \\\
#   --minsamples ', min_samples,' \\\
#   --significanceThreshold ', FDR_threashold,' \\\
#   --foldchangeThreshold ', FC_threashold,' \\\
#   > ', path_output, '/DIEGO_output_', tissue, '.txt')
#   filename = paste0(path_input,'/run_diego_', tissue,'.sh')
#   writeLines(diego_bash_file_content, filename)
#   system(paste("chmod +x", filename))
# }

makeDiegoInputFiles = function(rse.filtered, tissue, age_group, path_input,
                               condition_col_name){
  makeAFile(rse.filtered, tissue, age_group, path_input)
  makeBFile(rse.filtered, tissue, age_group, path_input, condition_col_name)
  # makeDiegoBashFile(tissue, path_input, path_output, reference_condition,
  #                   min_support, min_samples, FDR_threashold, FC_threashold)
}

readDiegoOutput = function(path_output, tissue, age_group){
  file_path = paste0(path_output, '/DIEGO_output_', paste(tissue, collapse = "_"), "_",
                     paste(age_group, collapse = "_"), '.txt')
  diego_output = read.delim(file_path, sep = "\t")
  diego_output
}

#' conda environment should be created before running this function and
#' DIEGO should be installed
#' conda create -n DIEGO_1 numpy=1.9 scipy matplotlib
#' wget http://legacy.bioinf.uni-leipzig.de/Software/DIEGO/DIEGO.tar.gz
#' tar -xzf DIEGO.tar.gz
runDiego = function(rse.filtered, tissue, age_group, reference_condition, path_input, path_output,
                    min_support, min_samples, FC_threshold, FDR_threshold, condition_col_name){
  makeDiegoInputFiles(rse.filtered, tissue, age_group, path_input, condition_col_name)
  use_condaenv("DIEGO_1")
  # considering that DIEGO counts abundance change as control-test
  conditions = unique(rse.filtered@colData[,condition_col_name])
  test_condition = setdiff(conditions, reference_condition)
  print(test_condition)
  system2(py_exe(), c("/home/an/DIEGO/diego.py",
                      paste0('-a ', path_input, '/junction_table_', paste(tissue, collapse = "_"),
                             "_", paste(age_group, collapse = "_"), '.txt'),
                      paste0('-b ', path_input, '/group_table_', paste(tissue, collapse = "_"),
                             "_", paste(age_group, collapse = "_"), '.txt'),
                      paste0('-x ', test_condition),
                      paste0('--minsupp ', min_support),
                      paste0('--minsamples ', min_samples),
                      paste0('--foldchangeThreshold ', FC_threshold),
                      paste0('--significanceThreshold ', FDR_threshold),
                      paste0('> ', path_output, '/DIEGO_output_', paste(tissue, collapse = "_"),
                             "_", paste(age_group, collapse = "_"), '.txt')
  )
  )
  diego.output = readDiegoOutput(path_output, tissue, age_group)
  diego.output
}

#===================================DJexpress==================================
# replace +/- strand with 1/2 (0: undefined, 1: +, 2: -)
# STAR manual, p.12 on output splice junctions file
makeDjeCoordinates <- function(coordinates_vector) {
  strand_to_numb_dict <- c("+" = 1, "-" = 2, "*" = 0)
  strand_signs <- sub(".*(?=.$)", "", coordinates_vector, perl = TRUE)
  coordinates_vector =  sub("([0-9]{3,})(-)([0-9]{3,})", "\\1:\\3", coordinates_vector)
  coordinates_vector = sapply(seq_along(coordinates_vector), function(i) {
    sub(".$", strand_to_numb_dict[strand_signs[i]], coordinates_vector[i])
  })
  coordinates_vector
}

# instead of DJEimport() etc, because input data differ (recount3 jxns instead of STAR raw out file)
makePrepOutObj = function(rse.filtered, tissue, reference_sample_ids){
  JunctExprfilt = rse2countDf(rse.filtered)
  rownames(JunctExprfilt) = makeDjeCoordinates(rownames(JunctExprfilt))
  featureID = rownames(JunctExprfilt)
  groupID = rse.filtered@rowRanges$gene_name
  
  isFromGroup1 = colnames(JunctExprfilt) %in% reference_sample_ids
  design       = cbind(1, ifelse(isFromGroup1, 0, 1))
  
  list(JunctExprfilt=JunctExprfilt, featureID=featureID, groupID=groupID, design=design)
}


runDJExpress = function(rse.filtered, tissue, condition_col_name, reference_condition, FDR_threshold, logFC_threshold){
  reference_sample_ids = findConditionIds(rse.filtered, condition_col_name, reference_condition)
  prep_out = makePrepOutObj(rse.filtered, tissue, reference_sample_ids)
  anlz_out = DJEanalyze(prepare.out = prep_out,
                         Group1 = reference_sample_ids,
                         FDR = FDR_threshold,
                         logFC = logFC_threshold)
  anlz_out
}


# ================================SAJR=========================================
makeSites = function(gene.jxns.info){ # junxtions of a gene
  gene.jxns.info = unique(as.data.frame(gene.jxns.info)) # choosing only unique rows nothing changes
  if(nrow(gene.jxns.info)==0)
    return(NULL) # if there are no junctions, return null
  
  # adding extra columns
  gene.jxns.info$rightmost = gene.jxns.info$leftmost =  NA  # adding 2 columns for rightmost and leftmost coordinates of a junction
  gene.jxns.info$id = rownames(gene.jxns.info) # adding a column with coordinates of each junction for gene i ("chrX:71910743-71913070:+" "chrX:71910743-71958862:+")
  
  # dublicating df
  jxns.same.start.info = jxns.same.end.info = gene.jxns.info  # duplicating jxns dataframe
  jxns.same.start.info$side = 'l' # adding a column 'side'
  jxns.same.end.info$side = 'r'
  
  # in duplicated dfs formatting rownames
  rownames(jxns.same.start.info) = paste0(rownames(jxns.same.start.info),':',jxns.same.start.info$gene_name,':l')
  rownames(jxns.same.end.info) = paste0(rownames(jxns.same.end.info),':',jxns.same.end.info$gene_name,':r')
  
  # forming a list, where key is a junction, and corresponding element is all junctions with same start/end coordinate
  jxns.same.coordinate.list = list()
  # filling in rigtmost and leftmost columns and the list
  for(junction in 1:nrow(gene.jxns.info)){ # number of rows in gene.jxns.info
    # same start (including the junction!)
    # finding junctions with the start coordinate same to the junction
    jxns.same.start.tf = gene.jxns.info$start == gene.jxns.info$start[junction]  # true/false vector for the junction
    # adding to the list [the junction coordinate] as a key, and all junctions that have the same start as the junction
    jxns.same.coordinate.list[[rownames(jxns.same.start.info)[junction]]] =
      rownames(gene.jxns.info)[jxns.same.start.tf] # assigning junction.info rownames! They fill further be used to select rows from counts df.
    
    # filling in leftmost column = start coordinate
    jxns.same.start.info$leftmost[junction]  = unique(gene.jxns.info$start[jxns.same.start.tf])
    # searching for the rightmost coordinate amongst junctions
    jxns.same.start.info$rightmost[junction] =  max(gene.jxns.info$end[jxns.same.start.tf])
    
    # same end (including the junction!)
    jxns.same.end.tf = gene.jxns.info$end == gene.jxns.info$end[junction]
    jxns.same.coordinate.list[[rownames(jxns.same.end.info)[junction]]] =
      rownames(gene.jxns.info)[jxns.same.end.tf]
    
    jxns.same.end.info$leftmost[junction]  = min(gene.jxns.info$start[jxns.same.end.tf])
    jxns.same.end.info$rightmost[junction] = unique(gene.jxns.info$end[jxns.same.end.tf])
  }
  
  # concatenating dataframes
  jxns.same.start.or.end.info = rbind(jxns.same.start.info, jxns.same.end.info)
  
  list(jxns.same.start.or.end.info = jxns.same.start.or.end.info,
       jxns.same.coordinate.list = jxns.same.coordinate.list[rownames(jxns.same.start.or.end.info)]) # elements to jxns.same.coordinate.list were asigned l after r. However, in df junctions.info.start.end.fixated first there is information for all l, than for all r. To make it corresponding, we will reorder the list
}


makeSAJRgene = function(rse.gene.filtered, min.cov=10){ # min.cov was set based on binomial dispersion
  jxns.sites.list = makeSites(rse.gene.filtered@rowRanges) # list, contaning jxns_s and js2j
  gene.jxns.counts = as.matrix(rse.gene.filtered@assays@data$counts)
  
  if(is.null(jxns.sites.list))
    return(NULL)
  jxns.same.start.or.end.info = jxns.sites.list$jxns.same.start.or.end.info
  jxns.same.coordinate.list = jxns.sites.list$jxns.same.coordinate.list
  
  # inclusion ratio calculation
  # -- creating count matrixes
  
  # inclusion segments - reads that map to the segment itself (cassette exon)
  # exclusion segments - reads that map to junction between upstream and downstream segments
  # all segments - reads that map to junction between upstream and downstream segments + reads that map to the junction
  jxns.incl.count.mtrx = jxns.incl.and.excl.count.mtrx =
    matrix(0, nrow=nrow(jxns.same.start.or.end.info), ncol=ncol(gene.jxns.counts))
  colnames(jxns.incl.count.mtrx) = colnames(jxns.incl.and.excl.count.mtrx) = colnames(gene.jxns.counts) # samples
  rownames(jxns.incl.count.mtrx) = rownames(jxns.incl.and.excl.count.mtrx) = rownames(jxns.same.start.or.end.info)
  
  # -- filling in count matrixes
  for(junction in 1:nrow(jxns.same.start.or.end.info)){
    # counts for each sample are stored in gene_jxns_counts. New matrix's columns are asigned names same to colnames(gene_jxns_counts).
    # when filling in matrix, counts are taken from gene_jxns_counts df, so counts for samples in matrix and df gene_jxns_counts are filled in correctly.
    
    # selecting all raw of counts for junctions (number of inclusion)
    jxns.incl.count.mtrx[junction,] = gene.jxns.counts[jxns.same.start.or.end.info$id[junction], ] # id column contains same junction ids as in count matrix
    jxns.incl.and.excl.count.mtrx[junction,] =
      apply(gene.jxns.counts[jxns.same.coordinate.list[[junction]], , drop=F], 2, sum)
    # drop=F prevents R from converting single column to vector. But here is not a single column?
    # rows from count.gene df are selected, and we sum counts for excluded (alternative) junctions for each column (sample) - 2 (apply to column), sum (function to apply)
  }
  
  # inclusion ratio - the proportion of transcripts that contains a segment
  jxn.count.incl.ratio.mtrx = jxns.incl.count.mtrx/jxns.incl.and.excl.count.mtrx # matrix
  jxn.count.incl.ratio.mtrx[jxns.incl.and.excl.count.mtrx < min.cov] = NA # frequency of inclusion is not defined (min.cov = 10, set based on binomial distribution dispersion)
  
  # exclusion counts. Substract included junctions from all
  jxn.count.excl.matr = jxns.incl.and.excl.count.mtrx - jxns.incl.count.mtrx
  
  sajr.gene = list(seg=jxns.incl.and.excl.count.mtrx, i=jxns.incl.count.mtrx,
                   e=jxn.count.excl.matr, ir=jxn.count.incl.ratio.mtrx)
  class(sajr.gene)='sajr'
  sajr.gene
}


makeSAJR = function(rse.filtered){
  sajr = NULL
  gene.ids = unique(rse.filtered@rowRanges$gene_id)
  # for each cytoskeleton gene
  for(gene.id in gene.ids){
    # selecting information only for the gene
    rse.gene.filtered = rse.filtered[rse.filtered@rowRanges$gene_id==gene.id,]
    # making sajr object for the geme
    sajr.gene =  makeSAJRgene(rse.gene.filtered)
    
    if(is.null(sajr.gene))
      next
    
    # making combined sajr object for all genes
    if(is.null(sajr))
      sajr = sajr.gene
    else{
      for(element in names(sajr)){ # combining seg with seg, e with e, etc for each element, for every gene
        sajr[[element]] = rbind(sajr[[element]],sajr.gene[[element]])}
    }
  }
  sajr
}

calculateMetrics = function(sajr, reference.sample.ids){
  reference.indices = match(reference.sample.ids, colnames(sajr$ir))
  dPSI = apply(sajr$ir, 1, function(x) mean(x[-reference.indices],na.rm=T)-mean(x[reference.indices],na.rm=T))
  #logFC = apply(sajr$ir, 1, function(x) log2(mean(x[adult.samples.ids],na.rm=T)/mean(x[fetus.samples.ids],na.rm=T)))
  dPSI #, logFC=logFC)
}

runSAJR = function(rse.filtered, tissue, condition_col_name, reference_condition){
  sajr.tissue = makeSAJR(rse.filtered)
  mod = rse.filtered@colData[,condition_col_name]
  mod = list(group=factor(mod)) # ~ model data
  reference.sample.ids = findConditionIds(rse.filtered, condition_col_name, reference_condition)
  
  alt.glm = as.data.frame( fitSAGLM(sajr.tissue, terms(x ~ group, keep.order=T),mod,return.pv=T) )
  alt.glm$dPSI = calculateMetrics(sajr.tissue, reference.sample.ids)
  alt.glm$FDR.sajr = p.adjust(alt.glm$group,m='BH')
  names(alt.glm)[names(alt.glm) == "group"] = "p.value.sajr"
  
  alt.glm
}
#==========================================================================================================
#============================Running tools and processing outputs==========================================
#==========================================================================================================
runTools = function(rse, tissue, age_group=c('fetus','adult'), 
                    condition_col_name='age_group', reference_condition='fetus',
                    path_input = '/home/an/Documents/GitHub/thesus_cytosk/DIEGO_input',
                    path_output = '/home/an/Documents/GitHub/thesus_cytosk/DIEGO_output',
                    min_support = 1, # minimum jxn count for a splice site to be considered.
                    min_samples = 1, #  minimum number of samples that must show the minimum support
                    FC_threshold = 1, # ? ratio of read counts of a splice junction in one condition compared to another
                    FDR_threshold = 0.05, # adjusted p-value threshold
                    logFC_threshold = 0
){
  rse.filtered = filterRse(rse, tissue, age_group)
  diego.output = runDiego(rse.filtered, tissue, age_group, reference_condition, path_input, path_output,
                          min_support, min_samples, FC_threshold, FDR_threshold, condition_col_name)
  print('diego running complete')
  dje.output = runDJExpress(rse.filtered, tissue, condition_col_name, reference_condition, FDR_threshold, logFC_threshold)
  
  print('dje running complete')
  sajr.output = runSAJR(rse.filtered, tissue, condition_col_name, reference_condition)
  print('sajr running complete')
  list(diego.output=diego.output, dje.output=dje.output, sajr.output=sajr.output)
}

