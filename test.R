# the project goal is to search for AS events in development
# DIEGO, DJexpress and SAJR tools are used for that (work with junction counts)
# in this file functions to make input for those tools are gathered
library(recount3)
library(reticulate)


rse2countDf = function(rse){
  counts = as.matrix(assay(rse, "counts"))
  counts = as.data.frame(counts)
  counts
}

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
  system2(py_exe(), c("/home/an/DIEGO/diego.py",
                      paste0('-a ', path_input, '/junction_table_', paste(tissue, collapse = "_"),
                             "_", paste(age_group, collapse = "_"), '.txt'),
                      paste0('-b ', path_input, '/group_table_', paste(tissue, collapse = "_"),
                             "_", paste(age_group, collapse = "_"), '.txt'),
                      paste0('-x ', reference_condition),
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
  list(diego.output=diego.output)
}

#==============================================================================
#=============================== running tools ================================
#==============================================================================
# # -- reading files

rse.jxn.cytosk = readRDS('rds/rse.jxn.cytosk.rds') # CHANGE TO YOUR DIRECTORY!

unique.tissues = unique(rse.jxn.cytosk@colData$tissue)
outputs_tissue = list()
for (tissue in unique.tissues[[1]]){
  outputs_tissue[[tissue]] = runTools(rse.jxn.cytosk, tissue)
}

assay(rse.jxn.cytosk, 'counts')
