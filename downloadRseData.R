## In this file project ERP109002 rse data are downloaded, filtered by genes (only cytoskeleton genes are left)
## and gene names are added to rse@rowRanges doi: 10.1038/s41586-019-1338-5
library(recount3)

# merging age groups
fetus = c("4wpc", "5wpc", "6wpc", "7wpc", "8wpc", "9wpc",  "10wpc", "11wpc", "12wpc", "13wpc",
          "14wpc", "16wpc", "18wpc", "19wpc", "20wpc")
newborn = c('0dpb', '4dpb', '6dpb', '15dpb', '18dpb', '19dpb', '34dpb') # typically refers to a baby from birth up to about 2 months of age
infant = c('94dpb', '127dpb', '221dpb', '226dpb', '270dpb', '271dpb', '6mpb', '1ypb') # infant - child before the age of 12 months
toddler = c('2ypb') # child between 2 and 3 years
school = c('4ypb', '7ypb','8ypb')
teen = c("13ypb","14ypb","16ypb","17ypb" ) # 13-19 years
yo_25_35 = c("25ypb","28ypb","29ypb","32ypb" )
yo_36_45 = c("39ypb" )
yo_46_55 = c( "46ypb", "50ypb", "53ypb" , "54ypb", "55ypb")
yo_56_65 = c("58ypb")

toddler.to.adult = c(toddler, school, teen, yo_25_35, yo_36_45, yo_46_55, yo_56_65)
adult = c(yo_25_35, yo_36_45, yo_46_55, yo_56_65)

cytoskeleton.genes = rbind(data.frame(gene_name=c('CYFIP1','CYFIP2','NCKAP1','NCKAP1L','ABI1','ABI2','ABI3','WASF1','WASF2','WASF3','BRK1'),group='WAVE'),
                           data.frame(gene_name=c('NHS','NHSL1','NHSL2','KIAA1522'),group='NHS'),
                           data.frame(gene_name=c('ARPC1A','ARPC1B','ARPC2','ARPC3','ARPC4','ARPC5','ACTR2','ACTR3','ACTR3B'),group='Arp2/3'))

downloadRse = function(project.id, type){
  human_projects = recount3::available_projects()
  proj_info = subset(human_projects,
                     project == project.id & project_type == "data_sources")
  rse = create_rse(proj_info, type=type) 
  rse
}

filterRseGenes = function(gene_names_df=cytoskeleton.genes, rse.gene){
  # -- gene rse only with cytoskeleton genes --
  rse.gene = rse.gene[rse.gene@rowRanges$gene_name %in% gene_names_df$gene_name,]
  # adding group names
  rse.gene@rowRanges$group = gene_names_df$group[
    match(rse.gene@rowRanges$gene_name,gene_names_df$gene_name)]
  # leaving only genes on main chr
  rse.gene =
    rse.gene[startsWith(as.character(rse.gene@rowRanges@seqnames),'chr'),]
  # sorting
  rse.gene =
    rse.gene[order(rse.gene@rowRanges@seqnames, rse.gene@rowRanges@ranges), ]
  
  # subsetting seqinfo
  new_seqinfo = Seqinfo(seqnames = as.character(unique(rse.gene@rowRanges@seqnames)),
                        seqlengths = NA,
                        isCircular = NA)
  seqlevels(rse.gene@rowRanges) = seqlevels(new_seqinfo)
  seqinfo(rse.gene@rowRanges) = new_seqinfo
  rse.gene
}

normaliseCoutsCPM = function(rse.gene){
  # normalisation (CPM)
  # raw reads mapped to the transcript / sum counts per sample *10**6
  rse.gene@assays@data$cpm = sweep(x = rse.gene@assays@data$raw_counts, MARGIN = 2,
                                   STATS = colSums(rse.gene@assays@data$raw_counts), FUN = `/`)*10**6
  # MARGIN = 1 means row; MARGIN = 2 means column.
  # STATS - the value(s) that should be added or subtracted
  # FUN The operation that has to be done
  rse.gene
}

annotateJxns = function(rse.gene.filtered, rse.jxn){
  # --making junxtion "annotation"
  # intersecting ranges
  # selecting only cytoskeleton genes in rse.jxn and adding gene column to it (making jxn "annotation")
  # The function compares each interval in the query object with the intervals in the subject object and identifies any overlaps.
  overlaps = findOverlaps(query = rse.jxn@rowRanges,
                          subject = rse.gene.filtered@rowRanges,
                          type='within') # the start and end positions of the query range must fall within the start and end positions of the subject range.
  
  # selecting only jxns of cytoskeleton genes and assing gene info
  rse.jxn = rse.jxn[overlaps@from,]
  rse.jxn@rowRanges$gene_id = rse.gene.filtered@rowRanges$gene_id[overlaps@to]
  rse.jxn@rowRanges$gene_name = rse.gene.filtered@rowRanges$gene_name[overlaps@to]
  list(rse.gene = rse.gene.filtered, rse.jxn = rse.jxn)
}

removeJxnDublicates = function(rse.jxn){
  # -- removing duplicates
  # checking
  table(duplicated(rse.jxn@rowRanges@ranges@NAMES))
  nms = unique(rse.jxn@rowRanges@ranges@NAMES)
  rse.jxn = rse.jxn[nms,]
  rse.jxn
}

splitSraLibraryCol = function(rse){
  # metadata
  # from sample metadata extracting sra.library_name, containing information on tissue, age etc
  split.sra.library_name = strsplit(rse@colData$sra.library_name, "\\.")
  # create a matrix to store the extracted values
  ncols = lengths(split.sra.library_name)[1]
  matrix.sra.library_name = matrix(nrow = nrow(rse@colData), ncol = ncols)
  colnames(matrix.sra.library_name) = c("sample_id", "organism", "tissue", "age", "gender")
  
  for (i in 1:nrow(matrix.sra.library_name)) {
    matrix.sra.library_name[i,] = split.sra.library_name[[i]]
  }
  # add the extracted columns to the colData
  rse@colData = cbind(rse@colData, matrix.sra.library_name)
  rse
}

transformCS2weeks = function(rse){
  # -- CS - Carnegie Stage to weeks
  cs.to.days = c('CS13' = 32, 'CS14' = 33, 'CS15' = 36, 'CS16' = 39, 'CS17' = 41, 'CS18' = 44,
                 'CS19' = 46, 'CS20' = 49, 'CS21' = 51, 'CS22' = 53, 'CS23' = 56) # dpc
  # CS to weeks
  cs.to.w = round(cs.to.days/4)
  cs.to.w = paste(cs.to.w, "w", sep="") # adding w, that's how all values are stored in age column age+measure
  names(cs.to.w) = names(cs.to.days) # asigning CS stages in correspondance to each week
  
  # replacing CS with w
  cs.idx = which(rse@colData$age %in% names(cs.to.w)) # index of cs values to replace
  cs.values = rse@colData$age[cs.idx] # selecting all cs values
  rse@colData$age = replace(rse@colData$age, cs.idx, cs.to.w[cs.values])
  
  # replace w (weeks) with wpc (weeks post conseption) (w is wpc according to figure 1 in doi: 10.1038/s41586-019-1338-5)
  rse@colData$age = gsub('w$', 'wpc', rse@colData$age)
  rse
}

combineTissues = function(rse, tissue.pairs.to.replace.list){
  # combining some tissues
  for (tissue.pair in tissue.pairs.to.replace.list){
    rse@colData$tissue = sub(tissue.pair[1], tissue.pair[2], rse@colData$tissue)
  }
  rse
}

getAgeOrder = function(rse){
  # getting ordered list of ages
  measures_order <- c("wpc", "dpb", "mpb", "ypb")
  unique_ages = unique(rse@colData$age)
  ages_order = c()
  for (i in c("wpc", "dpb", "mpb", "ypb")){
    time_group_values = grep(paste0('^', i, '|', i, '$'), unique_ages, value=TRUE)
    period_numbers = as.numeric(gsub("\\D", "", time_group_values))
    ages_order = c(ages_order,
                   time_group_values[order(period_numbers)])
  }
  ages_order
}

addSpecificAgeGroupColomn = function(rse){
  
  rse@colData$age_group_specific[rse@colData$age %in% fetus]    = rse@colData$age[rse@colData$age %in% fetus]
  rse@colData$age_group_specific[rse@colData$age %in% newborn]  = "newborn"
  rse@colData$age_group_specific[rse@colData$age %in% infant]   = "infant"
  rse@colData$age_group_specific[rse@colData$age %in% toddler]  = "toddler"
  rse@colData$age_group_specific[rse@colData$age %in% school]   = "school"
  rse@colData$age_group_specific[rse@colData$age %in% teen]     = "teen"
  rse@colData$age_group_specific[rse@colData$age %in% yo_25_35] = "25-35 y.o."
  rse@colData$age_group_specific[rse@colData$age %in% yo_36_45] = "36-45 y.o."
  rse@colData$age_group_specific[rse@colData$age %in% yo_46_55] = "46-45 y.o."
  rse@colData$age_group_specific[rse@colData$age %in% yo_56_65] = "56-55 y.o."
  
  rse
}

addConditionColumn = function(rse){
  # -- Adding column age group
  # fetus
  # heart develops earlier than other organs, so period before birth was limited to first 10 weeks
  rse@colData$age_group[!(rse@colData$tissue %in% c('Heart', 'Ovary')) & rse@colData$age %in% fetus] = "fetus"
  rse@colData$age_group[rse@colData$tissue == 'Heart' & rse@colData$age %in% 
                          c("4wpc", "5wpc",  "6wpc" , "7wpc",  "8wpc",  "9wpc",  "10wpc")] = "fetus"
  
  # adult
  rse@colData$age_group[!(rse@colData$tissue %in% c('Testies', 'Kidney', 'Ovary')) &
                                     rse@colData$age %in% toddler.to.adult] = "adult"
  # including infant in adult group for kidney, since otherwise there are not enough of samples
  rse@colData$age_group[rse@colData$tissue == 'Kidney' & rse@colData$age %in% c(infant, toddler.to.adult)] = "adult"
  # excluding teen from Testis, considering developmental peculiarities of testis development
  rse@colData$age_group[rse@colData$tissue == 'Testis' & rse@colData$age %in% adult] = "adult"
  rse
}

formatAnnotation = function(rse, tissue.pairs.to.replace.list){
  rse = splitSraLibraryCol(rse)
  rse = combineTissues(rse, tissue.pairs.to.replace.list)
  rse = transformCS2weeks(rse)
  rse = addSpecificAgeGroupColomn(rse)
  rse = addConditionColumn(rse)
  # sorting by tissue (for Brain and Cerebelum to be together)
  rse = rse[,order(rse@colData$tissue)]
  rse
}

save2RDS = function(rse.gene, rse.jxn, path){
  # -- saving files
  saveRDS(rse.gene,'rse.gene.cytosk.rds')
  saveRDS(rse.jxn,'rse.jxn.cytosk.rds')
}

prepareGeneRseAssay = function(project.id, type='gene'){
  rse.gene = downloadRse(project.id, type='gene')
  rse.gene = normaliseCoutsCPM(rse.gene)
  rse.gene = filterRseGenes(cytoskeleton.genes, rse.gene)
  rse.gene
}

prepareRse = function(project.id = 'ERP109002', 
                      # 24 cytoskeleton genes
                      cytoskeleton.genes = rbind(data.frame(gene_name=c('CYFIP1','CYFIP2','NCKAP1','NCKAP1L','ABI1','ABI2','ABI3','WASF1','WASF2','WASF3','BRK1'),group='WAVE'),
                                                 data.frame(gene_name=c('NHS','NHSL1','NHSL2','KIAA1522'),group='NHS'),
                                                 data.frame(gene_name=c('ARPC1A','ARPC1B','ARPC2','ARPC3','ARPC4','ARPC5','ACTR2','ACTR3','ACTR3B'),group='Arp2/3')),
                      tissue.pairs.to.replace.list = list(c("Forebrain", "Brain"),
                                                          c("Hindbrain", "Cerebellum"),
                                                          c("KidneyTestis", "Kidney"))
                      ){
  rse.jxn = downloadRse(project.id, type='jxn')
  rse.gene = prepareGeneRseAssay(project.id, type='gene')

  rse = annotateJxns(rse.gene, rse.jxn)
  rse.jxn = rse$rse.jxn
  rse.gene = rse$rse.gene 
  
  rse.jxn = removeJxnDublicates(rse.jxn)
  rse.jxn = formatAnnotation(rse.jxn, tissue.pairs.to.replace.list)
  rse.gene = formatAnnotation(rse.gene, tissue.pairs.to.replace.list)
  rse.jxn = rse.jxn[,rse.jxn@colData$age_group %in% c('adult', 'fetus')]
  save2RDS(rse.gene, rse.jxn, path='./')
}

mergeRse = function(gene.rse.list){
  assay.data.list = lapply(gene.rse.list, function(x) assay(x, 'cpm'))
  merged.assay.data = do.call(cbind, assay.data.list)
  col.data.list = lapply(gene.rse.list, function(x) colData(x, 'cpm')[,'tissue',drop=F])
  merged.col.data = do.call(rbind, col.data.list)
  # Create the new RangedSummarizedExperiment
  merged_rse <- SummarizedExperiment(
    assays = list(cpm = merged.assay.data), 
    rowRanges = rowRanges(gene.rse.list[[1]]), 
    colData = merged.col.data
  )
  merged_rse
}
