# install.packages("ggVennDiagram")
library(ggplot2)
library(ggVennDiagram)
library(grid)
library(gridExtra)
library(eulerr)
library(RColorBrewer)


getNrowsNcols = function(outputs_tissue, tumor){
  nrow = length(outputs_tissue)
  ncol = colnames(outputs_tissue[[1]]$all.jxns.info)
  ncol = length(ncol)-4  # Remove first two columns
  if(tumor){
    ncol=ncol-2
  }
  # Create the plot matrix
  total_numb_of_plots = nrow*ncol
  layout_matrix = matrix(1:total_numb_of_plots, nrow = nrow, ncol = ncol, byrow = TRUE)
  layout_matrix = cbind(layout_matrix, rep(total_numb_of_plots+1, nrow))
  print(layout_matrix)
  layout_matrix
}



setPlotParameters <- function(bottom_page_margin=2, left_page_margin=1, top_page_margin=2, right_page_margin=0,
                              bottom_plot_margin=2.4, left_plot_margin=0, top_plot_margin=0, right_plot_margin=0,
                              title_axis_distance=2, axis_label_distance=0.7, axis_line_distance=0,
                              font_size = 0.7, tick_length=-0.4,
                              layout_matrix, pty="s") {
  plot.new()
  par(tcl=tick_length, # -0.2 specifies the length of the ticks as a fraction of the height of a line of text. A negative value means the ticks will point inwards towards the plot, while a positive value would make them point outwards.
      mar = (c(bottom_plot_margin, left_plot_margin, top_plot_margin, right_plot_margin) + 0.1),  # mar = c(bottom, left, top, right)
      oma = (c(bottom_page_margin, left_page_margin, top_page_margin, right_page_margin) + 0.1),  # Bottom, left, top, right
      mgp = c(title_axis_distance, axis_label_distance, axis_line_distance),  # Title, label, and line distances
      xpd=TRUE,
      pty = pty, # square plots
      cex.axis = font_size)
  layout(mat = layout_matrix)
}

col = c(sajr = "#984EA3",
        dje = "orange3",
        diego = "#5DADE2",
        'dje&sajr' = "#FF9900",
        'diego&sajr' = '#E8FF00',
        'diego&dje' = '#F781BF',
        'diego&dje&sajr' = "#4DAF4A")


col_tum = c('sajr.norm.tumor' = '#979A9A',
            sajr = "#8F00FF",
            dje = "darkorange4",
            diego = "deepskyblue",
            'dje&sajr' = "#FF5733",
            'diego&sajr' = '#FFFF00',
            'diego&dje' = '#FF00FF',
            'diego&dje&sajr' = "#7FFF00")

makeDotplots = function(tf, all.jxns, intersections, tissue, log = '', tumor){
  x_lim_dict = list('dPSI_sajr' = c(-1,1), 
                    'dPSI_gtex2tum' = c(-1,1),
                    'dPSI_norm2tum' = c(-1,1),
                    'logFC_dje' = c(-4,4), 
                    'abund_change_diego' = c(-3,3))
  
  ticks_dict = list('dPSI_sajr' = seq(-1,1,by=0.5), 
                    'dPSI_gtex2tum' = seq(-1,1,by=0.5), 
                    'dPSI_norm2tum' = seq(-1,1,by=0.5), 
                    'logFC_dje'=seq(-4,4,by=2), 
                    'abund_change_diego' = seq(-3,3,by=1.5))
  if (tumor){
    col1_name <- grep("tum", colnames(all.jxns[,tf]), value = TRUE)
    # Get the names of all other columns
    other_cols <- setdiff(colnames(all.jxns[,tf]), col1_name)
    # Create a list of combinations
    comb <- lapply(other_cols, function(col_name) {
      all.jxns[, c(col1_name, col_name)]
    })
  }
  
  else{
    comb = combn(all.jxns[,tf], 2, simplify = FALSE)
  }
  
  lapply(comb, function(x) {
    par.1 = colnames(x)[1]
    par.2 = colnames(x)[2]
    
    plot(x,
         xaxt = "n", yaxt='n',    # Suppress x-axis ticks and labels
         xlab = '', 
         log=log,
         xlim=x_lim_dict[[par.1]], ylim=x_lim_dict[[par.2]],
         type = 'p', col = 'lightgrey',
         pch = 16 , cex = 0.7, lwd = 1,
         bty = "L")

    
    # tick axis + titles
    if(all(colnames(x) %in% names(ticks_dict))){
      axis(1, at=ticks_dict[[par.1]], labels = FALSE)
      axis(2, at=ticks_dict[[par.2]], labels = TRUE)
      
      # Fit linear regression
      a = x[complete.cases(x),]
      lm_model = lm(a[, par.2] ~ a[, par.1])
      abline(lm_model, col = "#a72127", lwd = 1, xpd=FALSE)
      
      if (par("mfg")[1]==par("mfg")[3]){
        axis(1, at=ticks_dict[[par.1]], labels = TRUE)
        axis(2, at=ticks_dict[[par.2]], labels = TRUE)
      }

    } else {
      axis(1, labels = TRUE)
      axis(2, labels = TRUE)
    }

    if (par("mfg")[1]==par("mfg")[3]) {
      mtext(side=1, text = par.1, line = 2, cex= 0.7)
    }
    if (par("mfg")[2]==1) {
      mtext(side=2, text = tissue, line = 3.5, cex= 1)
    }
    # points
    if (tumor) col=col_tum
    for (ids in names(col)){
      i = which(all.jxns$junction_id_sajr %in% intersections[[ids]])
      points(x[i,], col=col[ids], pch = 16, cex=1.1)
      # gene labeles
      if (ids=='diego&dje&sajr'){
        jxns = x[i,]
        labels = all.jxns[i,'gene_name']
        # Add text labels to the points
        n_labels = min(5, length(labels))
        if (n_labels > 0) {
          # Calculate label positions based on modulo 3
          label_positions <- (seq_len(n_labels) - 1) %% 3 + 1
          label_positions <- c(4, 3, 1)[label_positions]  # Map 1 to right, 2 to below, 3 to above

          text(jxns[1:n_labels, par.1], jxns[1:n_labels, par.2],
               labels = labels[1:n_labels], pos = label_positions, cex = 0.8)
        }
      }
    }
    # Calculate Spearman correlation
    a = x[complete.cases(x),]

    result = cor.test(a[,par.1], a[,par.2], method = "spearman")
    corr.coef = round(result$estimate, digits=2)
    mtext(paste0('rho = ', corr.coef), side=3,
          col = "black", cex=0.5)
  })
}

plotGraphs = function(all.jxns, intersections, tissue, tumor){
  col.fdr.if = grepl("FDR", colnames(all.jxns))
  col.metrics.if = !grepl("FDR|gene|id", colnames(all.jxns))  
  makeDotplots(col.metrics.if, all.jxns, intersections, tissue=tissue, tumor=tumor)
  makeDotplots(col.fdr.if, all.jxns, intersections, tissue=tissue, log='xy',tumor=tumor)

}

plotEulerDiagram = function(outputs_tissue, title, thresholds_text){
  p=list()
  for (tissue in names(outputs_tissue)){
    intersections = Reduce(append, outputs_tissue[[tissue]]$sign.jxns.info.list$intersections)
    intersections = intersections[c('sajr',
                            'dje',
                            'diego',
                            'dje&sajr',
                            'diego&sajr',
                            'diego&dje',
                            'diego&dje&sajr')]
    n_jxns = nrow(outputs_tissue[[tissue]]$all.jxns.info)
    data = sapply(intersections, 
                  function(x) length(x))
    order = names(intersections)
    fit = euler(data,
                shape = "circle",   # Force circles
                control = list(area.prop = TRUE))
    p[[tissue]] = plot(fit,
                       quantities = TRUE,
                       fills = list(fill = col), 
                       newpage = FALSE,
                       legend = TRUE,       # Remove legend
                       main = list(label=paste0("#jxns = ",n_jxns),
                                   fontsize=7)
                      )
  }
  # rc = ceiling(sqrt(length(outputs_tissue)))
  do.call(grid.arrange, c(p, ncol = 1, top = title, bottom = thresholds_text))  
}

plotVennDiagram = function(outputs_tissue, title, thresholds_text){
  p=list()
  for (tissue in names(outputs_tissue)){
    n_jxns = nrow(outputs_tissue[[tissue]]$all.jxns.info)
    p[[tissue]] = ggVennDiagram(outputs_tissue[[tissue]]$sign.jxns.info.list$all.single.tool, 
                                label_alpha = 0, 
                                label_col = "white") +
      labs(title = paste0(tissue, '( #jxns = ', n_jxns, ' )' )) +
      theme(legend.position = "none")
  }
  rc = ceiling(sqrt(length(outputs_tissue)))
  do.call(grid.arrange, c(p, ncol = 2, top = title, bottom = thresholds_text))
}

setTitles = function(thresholds){
  thresholds_text = paste0('Output filtration thresholds. logFC >= ', thresholds$logfc_threshold,
                           ', dPSI >= ', thresholds$dpsi_threshold,
                           ', abundance change >= ', abund_change_threshold,
                           ', FDR <= ', fdr_threshold)
  thresholds_text
}

addLegend = function(labels, pch, col, pt.cex=1){
  plot(1, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend('center', inset = c(0, 0),  # Adjust inset as needed
         bty='n', xpd=TRUE,
         legend = labels,
         col = col, pch = pch, # Use filled circle and star 
         pt.cex=pt.cex, cex=1,
         horiz = FALSE, ncol=1)
}

plotResultsRepot = function(outputs_tissue, tumor=FALSE, file='', thresholds){
  layout_matrix = getNrowsNcols(outputs_tissue,tumor)
  setPlotParameters(layout_matrix = layout_matrix)
  for (tissue in names(outputs_tissue)){
    intersect = Reduce(append, outputs_tissue[[tissue]]$sign.jxns.info.list$intersections)
    plotGraphs(outputs_tissue[[tissue]]$all.jxns.info,intersect, tissue, tumor)
  }
  if (tumor==TRUE){
    tool_names = names(col_tum)[grep("tum", names(col_tum), invert = TRUE)]
    addLegend(labels=c(file, paste(tool_names, "&", file), "not significant"), 
              col=c(col_tum, 'lightgrey'), pch=16, pt.cex=2)
    title=paste0(file, " and development. Tool comparison. (Base conditions: before birth and norm accordingly)")
    mtext(side=3, text = title, outer=TRUE, cex= 0.7, line=1)
    }
  else{
    addLegend(labels=c(names(col), "not significant"), 
              col=c(col, 'lightgrey'), pch=16, pt.cex=2)
    title = "Development (before*-after birth). Tool comparison (Base condition: before birth)"
    mtext(side=3, text = title, outer=TRUE, cex= 0.9, line=1)  }
  thresholds_text = setTitles(thresholds)
  mtext(side=1,  text = thresholds_text, outer=TRUE, cex= 0.7, line=1)
  
  plotVennDiagram(outputs_tissue, title, thresholds_text)
  if (tumor==FALSE){
    plotEulerDiagram(outputs_tissue, title = title, thresholds_text = thresholds_text)
  }
}


plotFisherResults = function(fisher_results_tissues_list, thresholds){
  # Create the plot matrix
  nrow = length(fisher_results_tissues_list)
  ncol = 2
  plots = c(seq(1:nrow), rep(nrow+1,nrow))
  layout_matrix = matrix(plots, nrow = nrow, ncol = ncol, byrow = FALSE)
  print(layout_matrix)
  setPlotParameters(layout_matrix=layout_matrix,
                    pty = "m", left_page_margin=20, bottom_page_margin = 4,
                    right_page_margin = 15, bottom_plot_margin = 2)
  
  names = fisher_results_tissues_list[[1]]$tool_pair
  nbars = length(names)
  col=c(sajr = "#984EA3",
        dje = "orange3",
        diego = "#5DADE2",
        sajr = "#8F00FF",
        dje = "darkorange4",
        diego = "deepskyblue",
        sajr = "#8F00FF",
        dje = "darkorange4",
        diego = "deepskyblue")

  lapply(names(fisher_results_tissues_list), function(tissue) {
    fisher_results_tissue = fisher_results_tissues_list[[tissue]]
    odds_ratio = fisher_results_tissue$odds_ratio
    q_val = fisher_results_tissue$q_val
    
    bp = barplot(odds_ratio,
                 col= col,
                 ylim=c(-3,150),
                 ylab = 'odds ratio',
                 xaxt = "n",
                 cex.lab = 0.7,
                 las=2)
    
    # Add stars based on significance
    for (i in 1:length(q_val)) {
      if (q_val[i] <= 0.05) {
        text(bp[i], odds_ratio[i] + 15, "*", cex = 1.5) 
      }
    }
    axis(1, at=bp, labels = FALSE)
    if (par("mfg")[1]==nrow){
      text(bp, par("usr")[3] - 5, labels = names, srt = 15, 
           adj = c(1,1), xpd = TRUE, cex = 0.7)
    }
    mtext(side = 2, text = "odds ratio", line = 1.7, cex = 0.6) 
    mtext(side=2, text = tissue, cex= 0.7, line=2.7)
  })
  
  addLegend(labels=c(names, 'q-value <=0.05'),
            pch = c(rep(15, length(names)), 8),
            col=c(col, 'black'), pt.cex=1)
  
  thresholds_text = setTitles(thresholds)
  mtext(side=1,  text = thresholds_text, outer=TRUE, cex= 0.7, line=1)
  mtext(side=3, text = 'Fisher Exact Test Results for each pair of methods & condition comaparisons', outer=TRUE, cex= 0.7)
  
  
}




#============================================
#=================covs
#====================================

set.colors = function(jxn, covs){
  # colors
  cols = ifelse(sub(':.$', '', rownames(covs$juncs)) %in% jxn,
                'orange', 'skyblue')
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

