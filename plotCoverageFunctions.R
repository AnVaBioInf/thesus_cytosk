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


#============================================
#=================covs
#====================================

set.colors = function(jxn, covs){
  # colors
  cols = ifelse(sub(':.$', '', rownames(covs$juncs)) %in% jxn,
                'orange2', '#377EB8')
  covs$juncs$cols = cols
  covs
}


#' Plots read coverage
#'
#' @param r read coverage; output of \code{\link{getReadCoverage}}
#' @param min.junc.cov numeric, plots only junctions (introns) with coverage not less than \code{min.junc.cov}
#' @param min.junc.cov.f numeric, plots only junctions (introns) with coverage not less than \code{min.junc.cov.f} of maximal coverage in the region
#' @param plot.junc.only.within logical, plot only juction with both ends within the region, FALSE plots all junctions with at least one end within region. NA plot all junctions overlapping the region.
#' @param reverse reverse x coordinates
#' @param junc.col colour for junction line. Individual color could be specified for each junction
#' @param junc.lwd line width for jucntion line
#' @param ... other parameters for plot function
#'
#' @export

plotReadCov = function(condition.cov.list,
                       min.junc.cov=0,
                       min.junc.cov.f=0,
                       plot.junc.only.within=FALSE,
                       xlim,
                       reverse=FALSE,
                       #    junc.col='blue',
                       junc.lwd=3,
                       bottom.mar=0,...){
  
  condition.cov.list = filter.data(condition.cov.list, xlim, min.junc.cov,plot.junc.only.within, min.junc.cov.f)
  
  # create a graph
  plot(condition.cov.list$x,
       condition.cov.list$cov,
       t='n',
       xlab = 'Chromosome coordinate',
       ylab = 'Coverage',
       cex.axis = 0.7,
       cex.lab = 0.9,
       xlim=c(xlim[1]-500, xlim[2]+500),
       ...)
  
  # plot verticale lines
  polygon(condition.cov.list$x,
          condition.cov.list$cov,
          col = 'gray',
          border=NA)
  
  if(nrow(condition.cov.list$juncs)>0) {
    for(i in 1:nrow(condition.cov.list$juncs)){
      plotArc(condition.cov.list$juncs$start[i],
              condition.cov.list$juncs$end[i],
              condition.cov.list$juncs$counts[i],
              col=condition.cov.list$juncs$cols[i],
              lwd=junc.lwd)
    }
  }
}

#' Plots parabolic arc
#'
#' @param from,to x coordinates of arc
#' @param top highest point of arc
#' @param n number of points
#' @param y.base bottom coordinate of arc
#' @param ... other parameters of lines functoin
plotArc = function(from,to,top,n=100,y.base=0,...){
  len = to - from
  x = seq(from=0,to=len,length.out = n)
  y = x*4*top/len - x^2*(4*top/len^2) 
  # This equation represents a downward-facing parabola that starts at 0, reaches its peak at top, and ends at 0 again.
  lines(x+from,y+y.base,...)
}

# for every gene
runJunxtionPlot = function(common_sign_jxns){
  for (jxn in unique(common_sign_jxns$junction_id)){
    # setting number of plot rows to number of tissues where selected gene junctions were significant, but no more than 3
    # selecting significant junctions for the GENE in both, development and cancer
    sign.jxn.df = common_sign_jxns[common_sign_jxns$junction_id==jxn,]
    gene = unique(sign.jxn.df$gene_name)
    par(mfrow = c(min(length(unique(sign.jxn.df$tissue)),4),2), bty='n')
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
                   ', FDR.sajr=', round(sign.jxn.tissue$FDR_sajr, digits=2),
                   ', logFC.dje=', round(sign.jxn.tissue$logFC_dje, digits=2),
                   ', FDR.dje=', round(sign.jxn.tissue$FDR_dje, digits=2),
                   ', sign.diego=', round(sign.jxn.tissue$abund_change_diego, digits=2),
                   ', FDR.diego=', round(sign.jxn.tissue$FDR_diego, digits=2), "\n",
                   ' dPSI_gtex2tum =',round(sign.jxn.tissue$dPSI_gtex2tum, digits=2),
                   ', FDR_gtex2tum =',round(sign.jxn.tissue$FDR_gtex2tum, digits=2),
                   ', dPSI_norm2tum =',round(sign.jxn.tissue$dPSI_norm2tum, digits=2),
                   ', FDR_norm2tum =',round(sign.jxn.tissue$FDR_norm2tum, digits=2), ")")
      
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
}