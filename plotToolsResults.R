# install.packages("ggVennDiagram")
library(ggplot2)
library(ggVennDiagram)
library(grid)
library(gridExtra)
library(eulerr)

getNrowsNcols = function(outputs_tissue, tumor){
  nrow = length(outputs_tissue)
  cols = colnames(outputs_tissue[[1]]$all.jxns.info)
  cols = length(cols)-4  # Remove first two columns
  if(tumor){
    cols=cols-2
  }
  list(nrow=nrow, ncol=cols)
}

setPlotParameters <- function(bottom_page_margin=2, left_page_margin=1, top_page_margin=2, right_page_margin=0,
                              bottom_plot_margin=2.4, left_plot_margin=0, top_plot_margin=0, right_plot_margin=0,
                              title_axis_distance=2, axis_label_distance=0.7, axis_line_distance=0,
                              font_size = 0.7, tick_length=-0.4,
                              nrow, ncol) {
  par(tcl=tick_length, # -0.2 specifies the length of the ticks as a fraction of the height of a line of text. A negative value means the ticks will point inwards towards the plot, while a positive value would make them point outwards.
      mar = (c(bottom_plot_margin, left_plot_margin, top_plot_margin, right_plot_margin) + 0.1),  # mar = c(bottom, left, top, right)
      oma = (c(bottom_page_margin, left_page_margin, top_page_margin, right_page_margin) + 0.1),  # Bottom, left, top, right
      mgp = c(title_axis_distance, axis_label_distance, axis_line_distance),  # Title, label, and line distances
      xpd=TRUE,
      pty = "s", # square plots
      cex.axis = font_size)

  # Create the plot matrix
  total_numb_of_plots = nrow*ncol
  layout_matrix = matrix(1:total_numb_of_plots, nrow = nrow, ncol = ncol, byrow = TRUE)
  layout(mat = layout_matrix)
}

col = c(sajr = "#984EA3",
        dje = "#A65628",
        diego = "#377EB8",
        'dje&sajr' = "#FF7F00",
        'diego&sajr' = '#FFFF33',
        'diego&dje' = '#F781BF',
        'diego&dje&sajr' = "#4DAF4A")


makeDotplots = function(tf, all.jxns, intersections, tissue, log = '', grid, tumor){
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
      
      if (par("mfg")[1]==grid[[1]]){
        axis(1, at=ticks_dict[[par.1]], labels = TRUE)
        axis(2, at=ticks_dict[[par.2]], labels = TRUE)
      }

    } else {
      axis(1, labels = TRUE)
      axis(2, labels = TRUE)
    }

    if (par("mfg")[1]==grid[[1]]) {
      mtext(side=1, text = par.1, line = 2, cex= 0.7)
    }
    if (par("mfg")[2]==1) {
      mtext(side=2, text = tissue, line = 3.5, cex= 1)
    }

    # points
    for (ids in names(intersections)){
      i = which(all.jxns$junction_id_sajr %in% intersections[[ids]])
      points(x[i,], col=col[ids], pch = 16)

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

plotGraphs = function(all.jxns, intersections, tissue, grid, tumor){
  col.fdr.if = grepl("FDR", colnames(all.jxns))
  col.metrics.if = !grepl("FDR|gene|id", colnames(all.jxns))  
  makeDotplots(col.metrics.if, all.jxns, intersections, tissue=tissue, grid=grid, tumor=tumor)
  makeDotplots(col.fdr.if, all.jxns, intersections, tissue=tissue, log='xy', grid=grid,tumor=tumor)

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
    p[[tissue]] = ggVennDiagram(outputs_tissue[[tissue]]$sign.jxns.info.list$all.single.tool) +
      labs(title = paste0(tissue, '( #jxns = ', n_jxns, ' )' )) +
      theme(legend.position = "none") 
  }
  rc = ceiling(sqrt(length(outputs_tissue)))
  do.call(grid.arrange, c(p, ncol = rc, top = title, bottom = thresholds_text))
}



plotFisher = function(outputs_tissue){
  all.single.tool = lapply(outputs_tissue$sign.jxns.info.list$all.single.tool, function(x) length(x))
  pairs = lapply(outputs_tissue$sign.jxns.info.list$intersections$only.pair.tools, function(x) length(x))
  
  data = matrix(c(all.single.tool$diego-pairs$'diego&dje', 
                  pairs$'diego&dje', 
                  all.single.tool$dje-pairs$'diego&dje',
                  pairs$'diego&dje'), nrow = 2, ncol=2, byrow = TRUE)
  
  print(data)
  
  # Perform Fisher's exact test
  result <- fisher.test(data)

  # Print the results
  print(result)
  
}


plotResultsRepot = function(outputs_tissue, tumor=FALSE, file='', thresholds){
  grid = getNrowsNcols(outputs_tissue,tumor)
  setPlotParameters(nrow = grid[[1]], ncol = grid[[2]])
  for (tissue in names(outputs_tissue)){
    intersect = Reduce(append, outputs_tissue[[tissue]]$sign.jxns.info.list$intersections)
    intersect = intersect[c('sajr',
                            'dje',
                            'diego',
                            'dje&sajr',
                            'diego&sajr',
                            'diego&dje',
                            'diego&dje&sajr')]
    plotGraphs(outputs_tissue[[tissue]]$all.jxns.info,
               intersect, tissue, grid, tumor)
  }
  thresholds_text = paste0('Filtration thresholds. logFC >= ', thresholds$logfc_threshold,
                           'dPSI >= ', thresholds$dpsi_threshold,
                           'abundance change >= ', abund_change_threshold, 
                           'FDR >= ', fdr_threshold)
  if (tumor==TRUE){  
    title=paste0(file, " and development. Tool comparison")
    mtext(side=3, text = title, outer=TRUE, cex= 0.7, line=1)}
  else{  
    title = "Development. Tool comparison"
    mtext(side=3, text = title, outer=TRUE, cex= 0.9, line=1)  }
  mtext(side=1,  text = thresholds_text, outer=TRUE, cex= 0.7, line=1)  
  
  plotVennDiagram(outputs_tissue, title, thresholds_text)
  if (tumor==FALSE){
    plotEulerDiagram(outputs_tissue, title = title, thresholds_text = thresholds_text)
  }
  plotFisher(outputs_tissue[[tissue]])
}
