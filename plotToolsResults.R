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

library(eulerr)

setPlotParameters <- function(bottom_page_margin=0, left_page_margin=2, top_page_margin=1, right_page_margin=0,
                              bottom_plot_margin=1, left_plot_margin=1, top_plot_margin=2, right_plot_margin=1,
                              title_axis_distance=2, axis_label_distance=1, axis_line_distance=0,
                              font_size = 0.7, tick_length=-0.4,
                              outputs_tissue) {
  par(tcl=tick_length, # -0.2 specifies the length of the ticks as a fraction of the height of a line of text. A negative value means the ticks will point inwards towards the plot, while a positive value would make them point outwards.
      mar = (c(bottom_plot_margin, left_plot_margin, top_plot_margin, right_plot_margin) + 0.1),  # mar = c(bottom, left, top, right)
      oma = (c(bottom_page_margin, left_page_margin, top_page_margin, right_page_margin) + 0.1),  # Bottom, left, top, right
      mgp = c(title_axis_distance, axis_label_distance, axis_line_distance),  # Title, label, and line distances
      xpd=TRUE,
      pty = "s", # square plots
      cex.axis = font_size)

  # Create the plot matrix
  nrow = length(outputs_tissue)+1
  cols = colnames(outputs_tissue[[1]]$all.jxns.info)
  ncol = length(cols[-c(1,2, length(cols)-1,length(cols))])  # Remove first two columns
  total_numb_of_plots = nrow*ncol
  layout_matrix = matrix(1:total_numb_of_plots, nrow = nrow, ncol = ncol, byrow = TRUE)
  layout(mat = layout_matrix)
}

setPlotParameters(outputs_tissue = outputs_tissue)

makeDotplots = function(tf, all.jxns, intersections){
  x_lim_dict = list('dPSI_sajr' = c(-1,1), 'logFC_dje' = c(-4,4), 'abund_change_diego' = c(-3,3))
  x_ticks_dict = list('dPSI_sajr' = seq(-1,1,by=0.5), 'logFC_dje'=seq(-4,4,by=2), 'abund_change_diego' = seq(-3,3,by=1.5))
  
  comb = combn(all.jxns[,tf], 2, simplify = FALSE)
  
  lapply(comb, function(x) {
    print(colnames(x))
    print(x_lim_dict[colnames(x)[1]])
    plot(x)
         # xlim = x_lim_dict[colnames(x)[1]],
         # ylim = x_lim_dict[colnames(x)[2]])
    for (ids in names(intersections)){
      col = c(diego = "#377EB8",
              dje = "#A65628",
              sajr = "#984EA3",
              'diego&dje' = '#F781BF',
              'diego&sajr' = '#FFFF33',
              'dje&sajr' = "#FF7F00",
              'diego&dje&sajr' = "#4DAF4A")
      i = which(all.jxns$junction_id_sajr %in% intersections[[ids]])
      points(x[i,], col=col[ids], pch = 16)
    }
  }
  )
}

plot_graphs = function(all.jxns, intersections){
  col.fdr.if = grepl("FDR", colnames(all.jxns))
  col.metrics.if = !grepl("FDR|gene|id", colnames(all.jxns))
  
  makeDotplots(col.metrics.if, all.jxns, intersections)
  makeDotplots(col.fdr.if, all.jxns, intersections)

}


plotVienn = function(intersections){
  data = sapply(intersections, function(x) length(x))
  fit = euler(data)
  plot(fit,
       quantities = TRUE)
}


for (tissue.output in outputs_tissue){
  plot_graphs(tissue.output$all.jxns.info, tissue.output$sign.jxns.info.list$intersections)
}
# 
# plots = list()
# for (tissue in names(outputs_tissue)){
#   plots[[tissue]] = plotVienn(outputs_tissue[[tissue]]$sign.jxns.info.list$intersections)
# }
# 
# library(gridExtra)
# 
# grid.arrange(plots[[1]], plots[[2]], ncol=6)

