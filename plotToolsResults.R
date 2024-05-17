# install.packages("ggVennDiagram")
library(ggplot2)
library(ggVennDiagram)
library(grid)
library(gridExtra)
library(eulerr)
library(RColorBrewer)


#==================================================================================
#=============================== results comparison================================
#==================================================================================
getNrowsNcols = function(ncol,nrow){
  # Create the plot matrix
  total_numb_of_plots = nrow*ncol
  layout_matrix = matrix(1:total_numb_of_plots, nrow = nrow, ncol = ncol, byrow = TRUE)
  layout_matrix = cbind(layout_matrix, rep(total_numb_of_plots+1, nrow))
  print(layout_matrix)
  layout_matrix
}

setPlotParameters <- function(bottom_page_margin=3, left_page_margin=2, top_page_margin=2, right_page_margin=3,
                              bottom_plot_margin=1.5, left_plot_margin=4, top_plot_margin=0, right_plot_margin=1,
                              title_axis_distance=2, axis_label_distance=0.7, axis_line_distance=0,
                              font_size = 0.9, tick_length=-0.4,
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

addSpearmanCorr = function(df){
  df = df[complete.cases(df),]
  df = df[order(df[,1]),]
  x = df[,1]
  y = df[,2]
  result = cor.test(x, y, method = "spearman")
  corr.coef = round(result$estimate, digits=2)
  p.value = ifelse(result$p.value <= 0.05, "<= 0.05", "> 0.05")
  mtext(paste0('rho = ', corr.coef,  ', p.val ', p.value), side=3,
        col = "black", cex=0.7)
}

addRegressionCurve = function(df){
  df = df[complete.cases(df),]
  df = df[order(df[,1]),]
  x = df[,1]
  y = df[,2]
  lm_model = lm(y ~ x)
  ci=predict.lm(lm_model,interval='confidence')
  c=col2rgb("#8b0000") 
  polygon(c(x,rev(x)),
          c(ci[,2],rev(ci[,3])),
          border=NA,
          col=rgb(c[1],c[2],c[3],0.2*255,maxColorValue = 255))
  lines(x, ci[,1], col = '#8b0000', lwd = 2)
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

# calculate corr only between significant values
makeDotplots = function(all.jxns, intersections, cols.tf, tissue, tumor, log, 
                        show_all_xtick_labels, 
                        add_regression_curve, add_spearman_corr,
                        axis_cex, title_cex, point_label_cex, col){
  
  lim_dict = list('dPSI' = list(lim = c(-1,1), ticks =  seq(-1,1,by=0.5)),
                  'logFC' = list(lim=c(-4,4), ticks = seq(-4,4,by=2)),
                  'abund_change' = list(lim=c(2,-2), ticks = seq(-2,2,by=1)),
                  'FDR'=list(lim=NULL))
  
  if (show_all_xtick_labels) par(mar = (c(4, 4, 0, 0) + 0.1))
  comb = getColumnCombinations(all.jxns[,cols.tf], tumor)
  lapply(comb, function(x) {
    par.1 = colnames(x)[1]
    par.2 = colnames(x)[2]
    xlims = lim_dict[sapply(names(lim_dict), function(tool) grepl(tool, par.1))][[1]]
    ylims = lim_dict[sapply(names(lim_dict), function(tool) grepl(tool, par.2))][[1]]
    
    plot(x,
         xaxt = "n", yaxt='n',    # Suppress x-axis ticks and labels
         xlab = '', ylab = gsub("_", " ", par.2),
         log=log,
         xlim=xlims$lim, ylim=ylims$lim,
         type = 'p', col = 'lightgrey',
         pch = 16 , cex = axis_cex, lwd = 1,
         bty = "L")
    # points
    for (tool in names(intersections)){
      sign.jxns.tool = which(all.jxns$junction_id_sajr %in% intersections[[tool]])
      points(x[sign.jxns.tool,], col=col[tool], pch = 16, cex=1.1)
      # gene labeles
      if (tool=='diego&dje&sajr'){
        # removing junctions with l/r end
        all = all.jxns[sign.jxns.tool,]
        all = all[order(all$dPSI_sajr),]
        all = all[!duplicated(all$junction_id), ]
        sign.jxns = x[rownames(all),]
        
        if (nrow(sign.jxns)==0) next
        labels = all$gene_name
        # Add text labels to the points
        text(sign.jxns[, par.1], sign.jxns[, par.2],
             labels = labels, pos = c(2,4), cex = point_label_cex, font = 2,
             xpd=TRUE)
        }
    }
    if (show_all_xtick_labels){
      axis(1, labels = TRUE)
      axis(2, labels = TRUE)
    } else {
      axis(1, at=xlims$ticks, labels = FALSE)
      axis(2, at=ylims$ticks, labels = TRUE)
      if (par("mfg")[1]==par("mfg")[3]){
        axis(1, at=xlims$ticks, labels = TRUE)
        axis(2, at=ylims$ticks, labels = TRUE)
      }
    }
    # adding x axis names (tool metrix)
    if (par("mfg")[1]==par("mfg")[3]) {
      mtext(side=1, text = par.1, line = 2, cex= axis_cex)
    }
    # adding tissue names to rows
    if (par("mfg")[2]==1) mtext(side=2, text = tissue, line = 3, cex= title_cex)  

    if (add_regression_curve) addRegressionCurve(x)
    if (add_spearman_corr) addSpearmanCorr(x)
      
  })
}

plotGraphs = function(outputs.prepr.list, cols.tf, tumor,  col, file, log='', 
                      show_all_xtick_labels = FALSE, 
                      add_regression_curve=TRUE, add_spearman_corr=TRUE,
                      axis_cex = 0.7, title_cex=1, point_label_cex=0.7, thresholds=''){
  nrow = length(outputs.prepr.list)
  ncol = ifelse(tumor, sum(cols.tf)-1, sum(cols.tf))
  layout_matrix = getNrowsNcols(ncol,nrow)
  setPlotParameters(layout_matrix = layout_matrix)
  
  for (tissue in names(outputs.prepr.list)){
    intersections = Reduce(append, outputs.prepr.list[[tissue]]$sign.jxns.info.list$intersections)
    makeDotplots(outputs.prepr.list[[tissue]]$all.jxns.info, 
                 intersections, cols.tf, tissue, tumor, log, show_all_xtick_labels, 
                 add_regression_curve, add_spearman_corr,
                 axis_cex, title_cex, point_label_cex, col=col)
  }
  if (tumor==TRUE){
    tool_names = names(col)[grep("tum", names(col), invert = TRUE)]
    addLegend(labels=c(file, paste(tool_names, "&", file), "not significant"),
              col=c(col, 'lightgrey'), pch=16, pt.cex=2)
    title=paste0(file, " and development. Tool comparison. (Base conditions: before birth and norm accordingly)")
    mtext(side=3, text = title, outer=TRUE, cex= title_cex, line=1)
  }
  else{
    tool_names = names(col)[grep("tum", names(col), invert = TRUE)]
    addLegend(labels=c(tool_names, "not significant"),
              col=c(col[2:length(col)], 'lightgrey'), pch=16, pt.cex=2)
    title = "Development (before*-after birth). Tool comparison (Base condition: before birth)"
    mtext(side=3, text = title, outer=TRUE, cex= title_cex, line=1)  }
  thresholds_text = setTitles(thresholds)
  mtext(side=1,  text = thresholds_text, outer=TRUE, cex= title_cex, line=2)
}


plotResultsRepot = function(outputs.prepr.list, tumor=FALSE, file='', thresholds,
                            metrics_png='metrics_plot.png', fdr_png='fdr_plot.png'){
  col = c('sajr.norm.tumor' = '#979A9A',
          sajr = "#984EA3",
          dje = "orange3",
          diego = "#5DADE2",
          'dje&sajr' = "#FF9900",
          'diego&sajr' = '#E8FF00',
          'diego&dje' = '#F781BF',
          'diego&dje&sajr' = "#4DAF4A")
  
  col.metrics.if = !grepl("FDR|gene|id", colnames(outputs.prepr.list[[1]]$all.jxns.info))
  col.fdr.if = grepl("FDR", colnames(outputs.prepr.list[[1]]$all.jxns.info))
  
  png(metrics_png, width = 25, height = 35, units = "cm", res = 700)
  plotGraphs(outputs.prepr.list=outputs.prepr.list, 
             cols.tf=col.metrics.if, tumor=tumor, thresholds=thresholds, col=col, file=file)
  dev.off()
  
  png(fdr_png, width = 25, height = 35, units = "cm", res = 700)
  plotGraphs(outputs.prepr.list=outputs.prepr.list, 
             cols.tf=col.fdr.if, tumor=tumor, log='xy', show_all_xtick_labels=TRUE,  
             add_regression_curve=FALSE, thresholds=thresholds, col=col, file=file)
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


