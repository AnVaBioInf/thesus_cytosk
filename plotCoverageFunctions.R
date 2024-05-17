
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