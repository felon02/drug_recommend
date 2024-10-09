# script_KMplot
# 2023.8.23

library(survminer)

# KM function  ------------------------------------------------------------

        
km_plot <- function(fit, df, ylab='Probability', xlab='Month', title='', pval_coord, 
                    legend_coord, legend_title, pat = colorJama){
  #
  fig_x <- ggsurvplot(fit, data = df,
                      ylab = ylab,
                      xlab = xlab,
                      size = 0.5,
                      censor.size = 2.5,
                      title = title,
                      pval = TRUE, pval.method = TRUE,    # Add p-value &  method name
                      pval.size = 2, pval.method.size=2,
                      pval.coord = pval_coord,
                      pval.method.coord=c(pval_coord[1]*0.3, pval_coord[2]),
                      surv.median.line = "hv",            # Add median survival lines
                      legend.title = legend_title,               # Change legend titles
                      legend = legend_coord,
                      # legend.labs = legend_labs,  # Change legend labels
                      palette = pat,                    # Use JCO journal color palette
                      font.main = c(8),
                      font.x = c(8),
                      font.y = c(8),
                      font.tickslab = c(6),
                      font.legend = c(6),
                      ggtheme = theme_survminer(base_size = 8, font.main = c(8)),
                      risk.table = TRUE,                  # Add No at risk table
                      #cumevents = TRUE,                   # Add cumulative No of events table
                      tables.height = 0.25,               # Specify tables height
                      tables.theme = theme_cleantable(font.main = 8),  # Clean theme for tables
                      tables.y.text = FALSE,              # Hide tables y axis text
                      fontsize = 2.5,                       # fontsize for the risk table and the cumulative events table.
                      #font.family = "Arial",              # font.family for the risk table and the cumulative events table.
  ) #+ guides(color=guide_legend(override.aes=list(fill=NA)))
  #
  return(fig_x)
}

# ?arrange_ggsurvplots

#
km_plot_arranged <- function(fit, df, ylab='Probability', xlab='Month', pval_coord, title='',
                             legend_coord, legend_title, pat = colorJama, km_table_heights=c(1.5, 0.5)){
  # fig
  fig_x <- km_plot(fit, df, ylab=ylab, xlab=xlab, pval_coord, title=title,
                   legend_coord, legend_title, pat = pat)
  #
  res = ggarrange(fig_x$plot+theme(legend.background=element_blank()), 
                  fig_x$table, ncol = 1, nrow = 2, heights = km_table_heights)
  #
  return(res)
} 

#
km_plot_oneCurve_figx <- function(fit, df, ylab='Probability', xlab='Month', title='',
                                  pat = colorJama, font_size = 6){
  fig_x = ggsurvplot(fit, 
                     data = df,
                     ylab = ylab,
                     xlab = xlab,
                     size = 0.5,
                     censor.size = 2.5,
                     title = title,
                     surv.median.line = "hv",            # Add median survival lines
                     palette = pat,                    # Use JCO journal color palette
                     font.main = c(font_size),
                     font.x = c(font_size),
                     font.y = c(font_size),
                     font.tickslab = c(font_size),
                     font.legend = c(font_size),
                     ggtheme = theme_survminer(base_size = font_size,
                                               # legend = "none" # useness
                                               font.main = c(font_size)),
                     risk.table = TRUE,                  # Add No at risk table
                     #cumevents = TRUE,                   # Add cumulative No of events table
                     tables.height = 0.25,               # Specify tables height
                     tables.theme = theme_cleantable(font.main = font_size),  # Clean theme for tables
                     tables.y.text = FALSE,              # Hide tables y axis text
                     fontsize = 2.5) 
  return(fig_x)
}
km_plot_oneCurve <- function(fit, df, ylab='Probability', xlab='Month', title='',
                             pat = colorJama, km_table_heights=c(1.5, 0.5)){
  # fig
  fig_x = km_plot_oneCurve_figx(fit, df, ylab=ylab, xlab=xlab, title=title, pat = pat)
  #
  res = ggarrange(fig_x$plot+theme(legend.position = "none"), 
                  fig_x$table, ncol = 1, nrow = 2, heights = km_table_heights)
  #
  return(res)
} 