
# install.packages("ggVennDiagram")
library(ggplot2)
library(ggVennDiagram)
library(grid)
library(gridExtra)
library(eulerr)

getNrowsNcols = function(outputs_tissue){
  nrow = length(names(outputs_tissue))
  cols = colnames(outputs_tissue[[1]])
  cols = length(cols[-c(1,2, length(cols)-1,length(cols))])  # Remove first two columns
  list(nrow=nrow, ncol=cols)
}

setPlotParameters <- function(bottom_page_margin=2, left_page_margin=1, top_page_margin=1.5, right_page_margin=0,
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


makeDotplots = function(tf, all.jxns, tissue, log = ''){
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
  
  #comb = combn(all.jxns[,tf], 2, simplify = FALSE)
  # Specify the target column
  target_col <- grep("tum", colnames(all.jxns[,tf]), value = TRUE)
  # Get all other column names
  other_cols <- setdiff(colnames(all.jxns[,tf]), target_col)
  # Generate combinations
  
  print(target_col)
  print(other_cols)
  
  for (i in other_cols){
    
  }
  lapply(other_cols, function(col) {
    print(col)
    print(target_col)
    # lapply(comb, function(x) {
    x = all.jxns[,c(col, target_col)]
    par.1 = col
    par.2 = target_col
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
      lm_model <- lm(x[, par.2] ~ x[, par.1], data = x[complete.cases(x), ])
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
    
    # Calculate Spearman correlation
    a = x[complete.cases(x),]
    
    result = cor.test(a[,par.1], a[,par.2], method = "spearman")
    corr.coef = round(result$estimate, digits=2)
    mtext(paste0('rho = ', corr.coef), side=3,
          col = "black", cex=0.5)
    
    
  })
}

plotGraphs = function(all.jxns, tissue){
  print(head(all.jxns))
  col.fdr.if_norm2tum = grepl("FDR", colnames(all.jxns)) & !grepl("_gtex2tum", colnames(all.jxns))
  col.fdr.if_gtex2tum = grepl("FDR", colnames(all.jxns)) & !grepl("_norm2tum", colnames(all.jxns))
  
  col.metrics.if_norm2tum = !grepl("FDR|gene|id|dPSI_gtex2tum", colnames(all.jxns))
  col.metrics.if_gtex2tum = !grepl("FDR|gene|id|dPSI_norm2tum", colnames(all.jxns))
  
  print(col.metrics.if_norm2tum)
  
  makeDotplots(col.metrics.if_norm2tum, all.jxns, tissue=tissue)
  #makeDotplots(col.metrics.if_gtex2tum, all.jxns, tissue=tissue)
  
  makeDotplots(col.fdr.if_norm2tum, all.jxns, tissue=tissue, log='xy')
  #makeDotplots(col.fdr.if_gtex2tum, all.jxns, tissue=tissue, log='xy')
  
}

# grid = getNrowsNcols(all.jxns.info.dev.cans.list)
# colnames(all.jxns.info.dev.cans.list[[1]])
setPlotParameters(nrow = 6, ncol = 6)
for (tissue in names(all.jxns.info.dev.cans.list)){
  print(tissue)
  plotGraphs(all.jxns.info.dev.cans.list[[tissue]], tissue)
}
