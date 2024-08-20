#====================================
#--------------coverage plots--------
#====================================

# converting bigwig file to junction coverage
bigWig2Cov = function(bw){
  bw = as.data.frame(bw) # columns: seqnames, start, end, width, strand, score
  bw = bw[order(bw$start),] # ordering bw df my start coordinate column
  start = bw$start[1] # starting coordinate of the gene
  stop = bw$end[nrow(bw)] # end coordinate of the gene
  r = rep(0,stop-start+1) # null vector with the length of a gene
  for(i in 1:nrow(bw)){ # for every record in bw
    r[(bw$start[i]:bw$end[i])-start+1] = bw$score[i] # assigning score (coverage) to every position of gene
  }
  list(cov=r,x=start:stop) # coverage, start and end gene coordinates on a chromosome
}


getRecountCov = function(sample.id, rse.jxn.filtered, path='/home/an/Manananggal/Input/bigWig/'){
  sample.path = paste0(path, sample.id, '.bw') # bw file name for a sample
  bw = rtracklayer::import.bw(sample.path, which=rse.jxn.filtered@rowRanges) # download and subset bw file for a sample
  r = bigWig2Cov(bw)  # get juntion coverage, start:stop gene coordinates on a chromosome
  
  r$juncs =  
    cbind(as.data.frame(rse.jxn.filtered@rowRanges)[,c('start','end','strand')],
                     score=rse.jxn.filtered@assays@data$counts[,sample.id])
  r
}

sumCovs = function(l){
  r = l[[1]] # read coverages for the first sample
  juncs = 
    unique(do.call(rbind,unname(lapply(l,function(c)c$juncs[,1:3]))))
  juncs$score = 0
  juncs[rownames(r$juncs),'score'] = r$juncs$score
  
  for(i in 2:length(l)){
    if(!all(r$x==l[[i]]$x))
      stop("objects should cover identicall intervals")
    r$cov = r$cov + l[[i]]$cov
    juncs[rownames(l[[i]]$juncs),'score'] = 
      juncs[rownames(l[[i]]$juncs),'score'] + l[[i]]$juncs$score
  }
  r$juncs = juncs
  r
}

loadEnsGTF = function(f,features=NULL){
  r = read.table(f,sep='\t')
  if(!is.null(features ))
    r = r[r$V3 %in% features,]
  a = lapply(strsplit(r$V9,';\\s?',perl=T),function(x){x=strsplit(x,'[ =]');setNames(sapply(x,'[',2),sapply(x,'[',1))})
  names = unique(unlist(lapply(a,names)))
  a = do.call(rbind,lapply(a,'[',names))
  colnames(a) = names
  r = r[,c(1:5,7)]
  colnames(r) = c('chr_id','type','feature','start','stop','strand')
  cbind(r,a)
}


gtf = loadEnsGTF('/home/an/Manananggal/Input/ref_annotation/filtered_ann_v26.gtf')
gtf[gtf$gene_name=='ABI1',]

plotTranscripts = function(a,
                           ylim=c(0,length(unique(a$transcript_id))),
                           xlim=c(ifelse(a$strand[1]=='+',min(a$start),max(a$stop)),ifelse(a$strand[1]=='+',max(a$stop),min(a$start))),
                           xlab=a$chr_id[1],
                           new=TRUE,yspace=0.8,exon.col='black',cds.col='black',
                           text.cex = 0.7,
                           ...){
  
  if(!is.na(exon.col))
    a$exon.col = exon.col
  if(!is.na(cds.col))
    a$cds.col = cds.col
  
  transc = split(a,a$transcript_id)
  transc = transc[order(sapply(transc,function(x){max(x$stop)-min(x$start)}))]
  if(new)
    plot(1,t='n',xlim=xlim,ylim=ylim,yaxt='n',ylab='',xlab=xlab,...)
  ystep = (ylim[2]-ylim[1])/length(transc)
  
  for(i in 1:length(transc)){
    y = ylim[1] + ystep*i - ystep/2
    t = transc[[i]]
    lines(c(min(t$start),max(t$stop)),c(y,y))
    f = t$feature == 'exon'
    if(sum(f)>0)
      rect(t$start[f],y-ystep/2*yspace,t$stop[f],y+ystep/2*yspace,col = 'white',border = t$exon.col[f])
    f = t$feature == 'CDS'
    if(sum(f)>0)
      rect(t$start[f],y-ystep/2*yspace,t$stop[f],y+ystep/2*yspace,col = t$cds.col[f],border = t$cds.col[f])
  }
  text(par('usr')[2], # Use the right x-coordinate of the plotting region
       seq(ylim[1] + ystep/2, by = ystep, length.out = length(transc)),
       sapply(transc, function(x) x$transcript_name[1]),
       adj = c(0, 0.5), # Adjust text to the left (0) and center vertically (0.5)
       xpd = TRUE, 
       cex = text.cex)
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

plotReadCov = function(r,min.junc.cov=0,min.junc.cov.f=0,plot.junc.only.within=FALSE,ylim=NULL,
                       xlim=range(r$x),reverse=FALSE,junc.col='blue',junc.lwd=3,bottom.mar=0,...){
  f = r$x >= xlim[1] & r$x <=xlim[2]
  r$x = r$x[f]
  r$cov = r$cov[f]
  if(nrow(r$juncs)>0)
    r$juncs$col = junc.col
  r$juncs = r$juncs[r$juncs$start <= xlim[2] & r$juncs$end >=xlim[1] & r$juncs$score >= min.junc.cov,]
  if(!is.na(plot.junc.only.within)){
    if(plot.junc.only.within){
      r$juncs = r$juncs[r$juncs$start > xlim[1] & r$juncs$end < xlim[2],]
    }else{
      r$juncs = r$juncs[(r$juncs$start > xlim[1] & r$juncs$start < xlim[2]) | (r$juncs$end > xlim[1] & r$juncs$end < xlim[2]),]
    }
  }
  
  start = r$x[1]
  end = r$x[length(r$x)]
  r$cov[c(1,length(r$cov))] = 0
  if(is.null(ylim)){
    ylim = c(0,max(r$cov,ifelse(nrow(r$juncs)>0,max(r$juncs$score),1)))
  }
  r$juncs = r$juncs[r$juncs$score >= min.junc.cov.f * ylim[2],]
  if(reverse)
    xlim=rev(xlim)
  plot(r$x,r$cov,t='n',ylim=ylim,xlim=xlim,yaxt='n', ...)
  axis(2,at=c(0,ylim[2]),labels = c('',ylim[2]))

  polygon(r$x,r$cov,col = 'gray',border=NA)
  if(nrow(r$juncs)>0)
    for(i in 1:nrow(r$juncs)){
      start = r$juncs$start[i]
      stop = r$juncs$end[i]
      # to make junction max height always within region
      if(start<xlim[1] & stop >= xlim[2]){
        start = xlim[1] - mean(xlim)
        stop = xlim[2] + mean(xlim)
      }else if (start < xlim[1]){
        start = max(start,xlim[1] - (stop - xlim[1]))
      }else if (stop > xlim[2]){
        stop = min(stop,xlim[2] + (xlim[2]-start))
      }
      plotArc(start,stop,r$juncs$score[i],col=r$juncs$col[i],lwd=junc.lwd)
    }
  invisible(ylim)
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

