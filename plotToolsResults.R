# install.packages("ggVennDiagram")
library(ggplot2)
library(ggVennDiagram)
library(grid)
library(gridExtra)
library(eulerr)
library(RColorBrewer)

order = c("4wpc", "5wpc", "6wpc", "7wpc", "8wpc", "9wpc",  "10wpc", "11wpc", "12wpc", "13wpc",
          "14wpc", "16wpc", "18wpc", "19wpc", "20wpc", "newborn", "infant", "toddler", "school", 
          "teen", "25-35 y.o.", "36-45 y.o.", "46-45 y.o.", "56-55 y.o.") # setting new order
order = setNames(1:length(order), order)

# --- plotting
tissue.col=c('Brain'="#3399CC",
             'Cerebellum'="#33CCFF",
             'Heart'="#CC0000",
             'Kidney'="#CC9900",
             'Liver'="#339900",
             'Ovary'="#CC3399",
             'Testis'="#FF6600",
             'BRCA' = 'darkgrey',
             'Breast normal'='white')


# ----------------------------------------------------------------------------------
# ------------------------ gene expression vs time graphs -------------------------
# ----------------------------------------------------------------------------------

findYlim = function(rse.gene.cytosk){
  cpm = as.data.frame(t(rse.gene.cytosk@assays@data$cpm))
  cpm$tissue = rse.gene.cytosk@colData$tissue
  
  a = by(cpm[, sapply(cpm, is.numeric)], cpm$tissue, 
         function(x) as.data.frame(
           apply(x, 2, quantile, probs = c(0.25, 0.75)), simplify = FALSE))
  a = do.call(cbind, a)
  max.up.quantile = max.col(a)[2]
  cpm.max.quantiles = a[,max.up.quantile]
  ylim.max = cpm.max.quantiles[2]+1.5*(cpm.max.quantiles[2]-cpm.max.quantiles[1])
  c(0,ylim.max)
}

fitCurve = function(gene.id, age.group.specific, gene.tissue.counts, col){
  # Fit a curve (example using loess)
  model = loess(gene.tissue.counts[[gene.id]] ~ as.numeric(age.group.specific))
  # Generate x-values for prediction
  x_pred = seq(min(age.group.specific), max(age.group.specific), length.out = 100)
  # Predict y-values using the model
  y_pred = predict(model, newdata = data.frame(age.group.specific = x_pred))
  # Add the curve to the plot
  lines(x_pred, y_pred, col = col, lwd = 2)  # Adjust lwd for line thickness
}

createGrid = function(gene.rse){
  numb.graphs = round(sqrt(length(gene.rse@rowRanges)))
  numb.empty = numb.graphs**2 - length(gene.rse@rowRanges)
  list(numb.graphs=numb.graphs, numb.empty=numb.empty)
}

setParams = function(gene.rse){
  numb = createGrid(gene.rse)
  par(mfrow = c(numb$numb.graphs, numb$numb.graphs),
      oma = c(5, 3, 1, 1),  # bottom, left, top, right.
      mar = c(1, 2, 1, 1),  # bottom, left, top, right.
      cex.axis = 0.7,  # Adjust the value as needed
      bty="l")
}

setAxis = function(x.value, gene.name, gene.rse){
  numb = createGrid(gene.rse)
  title(main = gene.name, line = 0.2, cex=0.8)  # Place title 3 lines above the plot
  axis(1, at = 1:length(x.value), labels = FALSE)  # x-axis ticks
  axis(2, labels = FALSE)  # y-axis ticks
  
  boxplot.coord = par("mfg")
  # Add y-axis only for rightmost plots
  if (boxplot.coord[2] == 1) {
    axis(2, las = 1)  
    mtext("CPM", side = 2, line = 2.4, las = 0, cex = par("cex.axis"))
  }
  # Add x-axis only for bottom plots
  if ((boxplot.coord[1] == 5) |
      (boxplot.coord[2] == numb$numb.graphs & boxplot.coord[1] == (numb$numb.graphs-numb$numb.empty)) ){
    axis(1, at = 1:length(x.value), labels = x.value, las = 2) 
  }
  grid(nx = NULL, ny = NULL)
}

plotScatterplotExpression = function(gene.rse){
  setParams(gene.rse)
  gene.ids.ordered = gene.rse@rowRanges[order(gene.rse@rowRanges$gene_name),
                                        c('gene_id', 'gene_name')]
  gene.ids.ordered = setNames(gene.ids.ordered$gene_id, gene.ids.ordered$gene_name)
  tissues = unique(gene.rse@colData$tissue)
  
  for (gene.name in names(gene.ids.ordered)){
    # Create a new plot for each gene
    gene.id = gene.ids.ordered[gene.name]
    plot(
      0, 0,  # Placeholder values, will be replaced by actual data
      xlim = c(1, length(order)),  # Set x-axis limits based on number of tissues
      ylim = findYlim(rse.gene.cytosk),  # Set y-axis limits based on gene expression
      type = "n",  # Start with an empty plot
      xaxt = "n",  # Suppress default x-axis
      yaxt = "n",  # Suppress default x-axis
      # main = gene.rse@rowRanges[gene.rse@rowRanges$gene_id==gene.id,]$gene_name  # Set title to gene name
    )
    setAxis(names(order), gene.name, gene.rse)
    
    for (tissue in tissues){
      samples = rownames(gene.rse@colData[gene.rse@colData$tissue==tissue,])
      col = tissue.col[tissue]
      age.group.specific = 
        gene.rse@colData[gene.rse@colData$tissue==tissue,]$age_group_specific
      age.group.specific = order[age.group.specific]
      
      gene.tissue.counts = as.data.frame(t(gene.rse@assays@data$cpm[gene.id,samples,drop=F]))
      gene.tissue.counts$age.group.specific = age.group.specific
      
      points(age.group.specific, gene.tissue.counts[[gene.id]], col=col)
      fitCurve(gene.id, age.group.specific, gene.tissue.counts, col=col)
    }
  }
  # Plot the legend in the last cell
  plot(x=0, y=0, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend('bottom', legend = tissues, col = tissue.col, pch = 16,
         bty = "n", y.intersp = 0.6, xpd = TRUE,
         inset = c(0, -0.2))
}

# boxplots
plotBoxplotExpression = function(gene.rse, xlab = "Tissue", ...){
  # Box: The box represents the interquartile range (IQR), which contains the middle 50% of the data. The bottom and top edges of the box correspond to the first quartile (Q1) and third quartile (Q3), respectively.
  # Median Line: A horizontal line inside the box that marks the median (Q2) of the data.
  # Whiskers: Lines extending from the box that represent the range of the data, excluding outliers.
  setParams(gene.rse)
  cpm = as.data.frame(t(gene.rse@assays@data$cpm))
  cpm$tissue <- gene.rse@colData[rownames(cpm),'tissue']
  tissues = unique(gene.rse@colData$tissue)
  exclude_indices <- grepl("brca|tumor|metastatic", tissues, ignore.case = TRUE)
  tissues = c(tissues[!exclude_indices], tissues[exclude_indices])
  cpm$tissue = factor(cpm$tissue, levels = tissues)
  gene.ids.ordered = gene.rse@rowRanges[order(gene.rse@rowRanges$gene_name),
                                        c('gene_id', 'gene_name')]
  # make a named list
  gene.ids.ordered = setNames(gene.ids.ordered$gene_id, gene.ids.ordered$gene_name)
  # Create a vector of colors for each box in the plot
  # Create boxplot
  for (gene.name in names(gene.ids.ordered)){
    gene.id = gene.ids.ordered[gene.name]
    boxplot(cpm[[gene.id]] ~ cpm[["tissue"]],
            xlab = xlab, ylab = "CPM",
            ylim = findYlim(gene.rse),
            las=2,
            xaxt = "n", yaxt = "n", 
            col=tissue.col[tissues])
    setAxis(tissues, gene.name, gene.rse)
  }
}













#==================================================================================
#=============================== results comparison================================
#==================================================================================
x_lim_dict = list('dPSI_sajr' = c(-1,1), 
                  'dPSI_gtex2tum' = c(-1,1),
                  'dPSI_norm2tum' = c(-1,1),
                  'logFC_dje' = c(-4,4), 
                  'abund_change_diego' = c(3,-3))

ticks_dict = list('dPSI_sajr' = seq(-1,1,by=0.5), 
                  'dPSI_gtex2tum' = seq(-1,1,by=0.5), 
                  'dPSI_norm2tum' = seq(-1,1,by=0.5), 
                  'logFC_dje'=seq(-4,4,by=2), 
                  'abund_change_diego' = seq(-3,3,by=1.5))

axis_names_dict = c('dPSI_sajr' = 'dPSI sajr', 
                    'dPSI_gtex2tum' = 'dPSI gtex2tum',
                    'dPSI_norm2tum' = 'dPSI norm2tum',
                    'logFC_dje' = 'logFC dje', 
                    'abund_change_diego' = 'AC diego')

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


getNrowsNcols = function(ncol,nrow){
  # Create the plot matrix
  total_numb_of_plots = nrow*ncol
  layout_matrix = matrix(1:total_numb_of_plots, nrow = nrow, ncol = ncol, byrow = TRUE)
  layout_matrix = cbind(layout_matrix, rep(total_numb_of_plots+1, nrow))
  print(layout_matrix)
  layout_matrix
}

setPlotParameters <- function(bottom_page_margin=3, left_page_margin=2, top_page_margin=2, right_page_margin=3,
                              bottom_plot_margin=1.5, left_plot_margin=4, top_plot_margin=0, right_plot_margin=0,
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

addSpearmanCorr = function(x,y){
  result = cor.test(x, y, method = "spearman")
  corr.coef = round(result$estimate, digits=2)
  p.value = ifelse(result$p.value <= 0.05, "<= 0.05", "> 0.05")
  mtext(paste0('rho = ', corr.coef,  ', p.val ', p.value), side=3,
        col = "black", cex=0.5)
}

addRegressionCurve = function(x,y){
  lm_model = lm(y ~ x)
  ci=predict.lm(lm_model,interval='confidence')
  c=col2rgb("red")
  polygon(c(x,rev(x)),
          c(ci[,2],rev(ci[,3])),
          border=NA,
          col=rgb(c[1],c[2],c[3],0.2*255,maxColorValue = 255))
  lines(x, ci[,1], col = "red", lwd = 1)
}

getColumnCombinations = function(df, tumor){
  if (tumor){
    col1_name = grep("tum", colnames(df), value = TRUE)
    # Get the names of all other columns
    other_cols = setdiff(colnames(df), col1_name)
    # Create a list of combinations
    comb = lapply(other_cols, function(col_name) {
      df[, c(col1_name, col_name)]
    })
  }
  else{
    comb = combn(df, 2, simplify = FALSE)
  }
  comb
}

makeDotplots = function(tf, all.jxns, intersections, tissue, tumor, log = ''){
  comb = getColumnCombinations(all.jxns[,tf], tumor)
  lapply(comb, function(x) {
    par.1 = colnames(x)[1]
    par.2 = colnames(x)[2]
    plot(x,
         xaxt = "n", yaxt='n',    # Suppress x-axis ticks and labels
         xlab = '', ylab = axis_names_dict[par.2],
         log=log,
         xlim=x_lim_dict[[par.1]], ylim=x_lim_dict[[par.2]],
         type = 'p', col = 'lightgrey',
         pch = 16 , cex = 0.7, lwd = 1,
         bty = "L")
    # points
    if (tumor) col=col_tum
    for (tool in names(col)){
      sign.jxns.tool = which(all.jxns$junction_id_sajr %in% intersections[[tool]])
      points(x[sign.jxns.tool,], col=col[tool], pch = 16, cex=1.1)
      # gene labeles
      if (tool=='diego&dje&sajr'){
        jxns = x[sign.jxns.tool,]
        if (nrow(jxns)==0) next
        print(jxns)
        labels = all.jxns[sign.jxns.tool,'gene_name']
        # Add text labels to the points
        text(jxns[, par.1], jxns[, par.2],
             labels = labels, pos = c(2,4), cex = 0.7, font = 2,
             xpd=TRUE)
        }
    }
    
    # tick axis + titles
    if(all(colnames(x) %in% names(ticks_dict))){
      axis(1, at=ticks_dict[[par.1]], labels = FALSE)
      axis(2, at=ticks_dict[[par.2]], labels = TRUE)
      if (par("mfg")[1]==par("mfg")[3]){
        axis(1, at=ticks_dict[[par.1]], labels = TRUE)
        axis(2, at=ticks_dict[[par.2]], labels = TRUE)
      }
    } else {
      axis(1, labels = TRUE)
      axis(2, labels = TRUE)
    }
    # adding x axis names (tool metrix)
    if (par("mfg")[1]==par("mfg")[3]) {
      mtext(side=1, text = axis_names_dict[par.1], line = 2, cex= 0.7)
    }
    # adding tissue names to rows
    if (par("mfg")[2]==1) {
      mtext(side=2, text = tissue, line = 3, cex= 1)
    }
    df = x[complete.cases(x),]
    df = df[order(df[,1]),]
    addSpearmanCorr(df[,1],df[,2])
    addRegressionCurve(df[,1],df[,2])
  })
}

plotGraphs = function(outputs.prepr.list, tumor, cols, log=''){
  nrow = length(outputs.prepr.list)
  ncol = sum(cols)
  
  layout_matrix = getNrowsNcols(ncol,nrow)
  setPlotParameters(layout_matrix = layout_matrix)
  
  for (tissue in names(outputs.prepr.list)){
    intersect = Reduce(append, outputs.prepr.list[[tissue]]$sign.jxns.info.list$intersections)
    makeDotplots(cols, outputs.prepr.list[[tissue]]$all.jxns.info, 
                 intersect, tissue=tissue, tumor=tumor, log)
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
  mtext(side=1,  text = thresholds_text, outer=TRUE, cex= 0.7, line=2)
}


plotResultsRepot = function(outputs.prepr.list, tumor=FALSE, file='', thresholds){
  
  col.fdr.if = grepl("FDR", colnames(outputs.prepr.list[[1]]$all.jxns.info))
  col.metrics.if = !grepl("FDR|gene|id", colnames(outputs.prepr.list[[1]]$all.jxns.info))
  
  png("metrics_plot.png", width = 30, height = 30, units = "cm", res = 700)
  plotGraphs(outputs.prepr.list, tumor, col.fdr.if)
  dev.off()
  
  png("fdr_plot.png", width = 30, height = 30, units = "cm", res = 700)
  plotGraphs(outputs.prepr.list, tumor, col.metrics.if, log='')
  dev.off()
  
  
  # plotVennDiagram(outputs.prepr.list, title, thresholds_text)
  # if (tumor==FALSE){
  #   plotEulerDiagram(outputs.prepr.list, title = title, thresholds_text = thresholds_text)
  # }
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
  legend('center',  # Adjust inset as needed
         bty='n', xpd=TRUE,
         legend = labels,
         col = col, pch = pch, # Use filled circle and star 
         pt.cex=pt.cex, cex=1,
         horiz = FALSE, ncol=1)
}





































#===========fisher

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


##################### gene expression
# -----------------------------------------------------------------------------------
# ------------------------- samples occurrence heatmap ------------------------------
# -----------------------------------------------------------------------------------
plotHeatmapSamples = function(){
  par(oma = c(2, 2, 0, 0))  # bottom, left, top, right.
  # occurrence of samples
  sample_occurance = as.data.frame.matrix(
    table(rse.gene.cytosk@colData$tissue, rse.gene.cytosk@colData$age_group_specific) )
  # The table() function takes these two vectors as input and creates a contingency table. This table shows the frequency distribution of cells across different combinations of tissue types and age groups.
  sample_occurance = sample_occurance[,names(order)]
  sample_occurance = t(sample_occurance)
  # # making a column for tissues
  # sample_occurance$tissue = rownames(sample_occurance)
  # sample_occurance
  
  # Set color palette (e.g., blue to red)
  colors <- colorRampPalette(c("white", "red"))(256)
  
  # Create the heatmap
  image(1:ncol(sample_occurance), 1:nrow(sample_occurance), 
        t(sample_occurance), col = colors, axes = FALSE, xlab = "", ylab = "")
  
  # Add text labels with data values
  text(x = col(sample_occurance), 
       y = row(sample_occurance), 
       labels = sample_occurance, 
       col = "black")
  
  # Add axes and labels
  axis(1, at = 1:ncol(sample_occurance), labels = colnames(sample_occurance), las = 2)
  axis(2, at = 1:nrow(sample_occurance), labels = rownames(sample_occurance), las = 2)
}



# plotHeatmapSamples()
# 
# plotScatterplotExpression(rse.gene.cytosk, tissue.col)
# plotBarplotExpression(rse.gene.cytosk, tissue.col)

# gene.ids.ordered = gene.info[order(gene.info$gene_name),]$gene_id
# tissues = unique(gene.rse@colData$tissue)
# ylim=findYlim(rse.gene.cytosk)
# 
# cpm = as.data.frame(t(gene.rse@assays@data$cpm))
# cpm$tissue = gene.rse@colData$tissue
# cpm$tissue = factor(cpm$tissue)

