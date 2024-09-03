## function 1 ==================================================================================
customize_labels <- function (p, font.title = NULL,
                              font.subtitle = NULL, font.caption = NULL,
                              font.x = NULL, font.y = NULL, font.xtickslab = NULL, font.ytickslab = NULL)
{
  original.p <- p
  if(is.ggplot(original.p)) list.plots <- list(original.p)
  else if(is.list(original.p)) list.plots <- original.p
  else stop("Can't handle an object of class ", class (original.p))
  .set_font <- function(font){
    font <- ggpubr:::.parse_font(font)
    ggtext::element_markdown (size = font$size, face = font$face, colour = font$color)
  }
  for(i in 1:length(list.plots)){
    p <- list.plots[[i]]
    if(is.ggplot(p)){
      if (!is.null(font.title)) p <- p + theme(plot.title = .set_font(font.title))
      if (!is.null(font.subtitle)) p <- p + theme(plot.subtitle = .set_font(font.subtitle))
      if (!is.null(font.caption)) p <- p + theme(plot.caption = .set_font(font.caption))
      if (!is.null(font.x)) p <- p + theme(axis.title.x = .set_font(font.x))
      if (!is.null(font.y)) p <- p + theme(axis.title.y = .set_font(font.y))
      if (!is.null(font.xtickslab)) p <- p + theme(axis.text.x = .set_font(font.xtickslab))
      if (!is.null(font.ytickslab)) p <- p + theme(axis.text.y = .set_font(font.ytickslab))
      list.plots[[i]] <- p
    }
  }
  if(is.ggplot(original.p)) list.plots[[1]]
  else list.plots
}


## function 2 ================================================================================== 
ggsurvplot_customised <- function(dataframe, #lung.surv
                                  pval_coord, #c(800, 0.75)
                                  break_time_by, #500
                                  conf_int, #FALSE or TRUE
                                  title, #Kaplan-Meier Curve for Lung Cancer 
                                  legend_title, #"chr 2q status"
                                  legend_labels, # c("amp", "wt")
                                  subtitle, 
                                  xlim = NULL){ #"OS, hh_ova (n = 115)" c("amp", "wt")
  ggsurv <- ggsurvplot(fit, data = dataframe, 
                       risk.table = TRUE, 
                       pval = TRUE, 
                       conf.int = conf_int, 
                       conf.int.style = "ribbon",
                       xlab = "Days", 
                       #pval.method = TRUE, 
                       break.time.by = break_time_by, 
                       ggtheme = theme_light(), 
                       risk.table.y.text.col = TRUE,
                       risk.table.height = 0.2, 
                       risk.table.y.text = FALSE, 
                       ncensor.plot = TRUE,
                       ncensor.plot.height = 0.20,  
                       #surv.median.line = "hv", 
                       legend.labs = legend_labels, 
                       test.for.trend = FALSE, 
                       title = title, 
                       subtitle = subtitle, 
                       legend.title = legend_title, 
                       legend = "top", 
                       surv.scale = "percent", 
                       pval.coord = pval_coord, 
                       pval.size = 3, 
                       xlim = xlim, 
                       fontsize = 3, 
                       size = 0.5, 
                       font.family = "Arial", 
                       font.legend = c(9), 
                       censor.size = 1, 
                       censor.shape = 124)
  ggsurv <- customize_labels(ggsurv, 
                             font.title = c(9, "bold"), 
                             font.subtitle = c(9, "italic", "darkgrey"),
                             font.x = c(9), 
                             font.y = c(9), 
                             font.xtickslab = c(8))
  return(ggsurv)
} 


