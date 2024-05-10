#source('ages_comparison.R')
#install.packages("RColorBrewer")
#dev.off()
library(RColorBrewer)
# display.brewer.pal(n = 8, name = "Spectral")
# set2_colors <- brewer.pal(n = 8, name = "Spectral")
# set2_colors

# легенда
# корреляционная кривая ?
# группы в djexpress?
# дубликаты значимые только в djexpress? раздваиваются?
# corr coef для p value
# дать название всему графику
# 
# # сделать таблицу значимых во всех методах
# # сделать реверсед таблицу

setPlotParameters <- function(bottom_page_margin=0, left_page_margin=2, top_page_margin=1, right_page_margin=0,
                              bottom_plot_margin=1, left_plot_margin=3, top_plot_margin=2, right_plot_margin=1,
                              title_axis_distance=2, axis_label_distance=1, axis_line_distance=0,
                              font_size = 0.7, tick_length=-0.2,
                              outputs_tissue) {
  par(tcl=-tick_length, # -0.2 specifies the length of the ticks as a fraction of the height of a line of text. A negative value means the ticks will point inwards towards the plot, while a positive value would make them point outwards.
      mar = (c(bottom_plot_margin, left_plot_margin, top_plot_margin, right_plot_margin) + 0.1),  # mar = c(bottom, left, top, right)
      oma = (c(bottom_page_margin, left_page_margin, top_page_margin, right_page_margin) + 0.1),  # Bottom, left, top, right
      mgp = c(title_axis_distance, axis_label_distance, axis_line_distance),  # Title, label, and line distances
      xpd=TRUE,
      pty = "s", # square plots
      cex.axis = font_size)

  # Create the plot matrix
  nrow = length(outputs_tissue)+1
  cols = colnames(outputs_tissue[[1]]$all.jxns.info)
  ncol = length(cols[-c(1,2, length(cols)-1,length(cols))])+1  # Remove first two columns
  total_numb_of_plots = nrow*ncol
  layout_matrix = matrix(1:total_numb_of_plots, nrow = nrow, ncol = ncol, byrow = TRUE)
  layout(mat = layout_matrix)
}

setPlotParameters(outputs_tissue = outputs_tissue)


library(eulerr)

plotVienn = function(intersections){
  layout.show()
  data = c(sapply(intersections, function(x) length(x)))
  
  # Example data
  fit <- euler(data)
  # Plot with numbers in intersections
  p = plot(fit, quantities = TRUE)
}

plot_graphs = function(all.jxns, intersections){
  x_lim_dict = list('dPSI.sajr' = c(-1,1), 'logFC.dje' = c(-4,4), 'abund_change.diego' = c(-3,3))
  x_ticks_dict = list('dPSI.sajr' = seq(-1,1,by=0.5), 'logFC.dje'=seq(-4,4,by=2), 'abund_change.diego' = seq(-3,3,by=1.5))
  
  col.fdr.if = grepl("FDR", colnames(outputs_tissue[[1]]$all.jxns))
  col.metrics.if = !grepl("FDR|gene|id", colnames(outputs_tissue[[1]]$all.jxns))
  
  # 1. Using negation with grep
  fdr_comb = combn(all.jxns[,col.fdr.if], 2, simplify = FALSE)
  metrics_comb = combn(all.jxns[,col.metrics.if], 2, simplify = FALSE)
  

  p = plotVienn(intersections)
  show(p)
  lapply(metrics_comb, function(x) {
    plot(x)
    for (ids in names(intersections)){
      col = c(diego.unique = "#377EB8",
              dje.unique = "#A65628",
              sajr.unique = "#984EA3",
              diego_and_dje = '#F781BF',
              diego_and_sajr = '#FFFF33',
              dje_and_sajr = "#FF7F00",
              all.tools = "#4DAF4A")
      i = which(all.jxns$junction_id_sajr %in% intersections[[ids]])
      points(x[i,], col=col[ids], pch = 16)
      }
    }
  )


  lapply(fdr_comb, function(x) {
    plot(x, log='xy')
    for (ids in names(intersections)){
      col = c(diego.unique = "#377EB8",
              dje.unique = "#A65628",
              sajr.unique = "#984EA3",
              diego_and_dje = '#F781BF',
              diego_and_sajr = '#FFFF33',
              dje_and_sajr = "#FF7F00",
              all.tools = "#4DAF4A")
      i = which(all.jxns$junction_id_sajr %in% intersections[[ids]])
      points(x[i,], col=col[ids], pch = 16)
    }
  })
  
  
}

# all = outputs_tissue[[1]]$all.jxns.info
# inter = outputs_tissue[[1]]$sign.jxns.info.list$intersections[[1]]
# 
# which(all$junction_id_sajr %in% inter)



for (tissue.output in outputs_tissue){
  plot_graphs(tissue.output$all.jxns.info, tissue.output$sign.jxns.info.list$intersections)

}

