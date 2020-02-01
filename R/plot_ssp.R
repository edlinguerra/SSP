## plot_ssp: Function to plot MultSE and sampling effort relationship of simulated data
#'@importFrom ggplot2 ggplot
#'@importFrom ggplot2 aes
#'@importFrom ggplot2 aes_
#'@importFrom ggplot2 geom_point
#'@importFrom ggplot2 geom_errorbar
#'@importFrom ggplot2 theme_bw
#'@importFrom ggplot2 ylab
#'@importFrom ggplot2 xlab
#'@importFrom ggplot2 rel
#'@importFrom ggplot2 scale_x_continuous
#'@importFrom ggplot2 scale_y_continuous
#'@importFrom ggplot2 theme
#'@importFrom ggplot2 element_text
#'@importFrom ggplot2 element_blank
#'@importFrom ggplot2 element_rect
#'@importFrom ggplot2 element_line
#'@importFrom ggplot2 annotate
#'@importFrom ggplot2 geom_text
#'@importFrom ggplot2 facet_grid
#'@importFrom ggplot2 geom_rect
#'@export

plot_ssp<-function(xx, opt, multi.site){
  if (multi.site == FALSE){
    ul<-round(max(xx$upper) + 0.1*max(xx$upper), 2)

    fig<-ggplot(xx, aes_(x=~samples, y=~mean))+
      geom_point(size = 0.5)+
      geom_errorbar(aes_(ymin=~lower, ymax=~upper), size=0.1, width=.2)+
      theme_bw(base_size=16) +
      ylab ("Multivariate pseudo SE")+
      xlab("Sampling effort (n)")+
      scale_y_continuous(breaks=seq(0.0, ul, 0.025))+
      scale_x_continuous(breaks=seq(min(xx$samples), max(xx$samples), 1))+
      theme(axis.text.x = element_text(colour="black", size=rel(0.7)),
            axis.text.y = element_text(colour="black", size=rel(0.7)),
            axis.title.x = element_text(colour="black", size=rel(0.9)),
            axis.title.y = element_text(colour="black", size=rel(0.9)),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(size=0.4),
            axis.ticks= element_line(size=0.2))+
      annotate("rect", xmin=opt[2], xmax=opt[3], ymin=min(xx$lower),
              ymax=max(xx$upper), alpha=.5, fill="grey10")+
      annotate("rect", xmin=opt[1], xmax=opt[2], ymin=min(xx$lower),
               ymax=max(xx$upper), alpha=.5, fill="grey50")+
      geom_text(aes_(x = ~samples, y = ~upper + 0.01, label = ~cum), na.rm = TRUE)
    return(fig)
  } else {

    shade.opt<-data.frame("xmin" = c(opt[1,2], opt[2,2]),
                          "xmax" = c(opt[1,3], opt[2,3]),
                          "ymin" = c(min(xx$lower), min(xx$lower)),
                          "ymax" = c(max(xx$lower), max(xx$lower)),
                          "sv" = c("sites", "samples"))

    shade.sub<-data.frame("xmin" = c(opt[1,1], opt[2,1]),
                          "xmax" = c(opt[1,2], opt[2,2]),
                          "ymin" = c(min(xx$lower), min(xx$lower)),
                          "ymax" = c(max(xx$lower), max(xx$lower)),
                          "sv" = c("sites", "samples"))

    ul<-round(max(xx$upper) + 0.1*max(xx$upper), 2)

    my_breaks<-function(x) {
      y<-seq(min(x), max(x), 1)
      y<-round(y,0)
    }

    fig<-ggplot(xx, aes_(x=~samples, y=~mean))+
      geom_point()+
      geom_errorbar(aes_(ymin=~lower, ymax=~upper), width=.1)+
      facet_grid(.~sv, scales = "free_x")+
      theme_bw(base_size=16) +
      ylab ("Multivariate pseudo SE")+
      xlab("Sampling effort")+
      scale_y_continuous(breaks=seq(0.0, ul, 0.025))+
      scale_x_continuous(breaks = my_breaks)+
      theme(axis.text.x = element_text(colour="black", size=rel(0.7)),
            axis.text.y = element_text(colour="black", size=rel(0.7)),
            axis.title.x = element_text(colour="black", size=rel(0.9)),
            axis.title.y = element_text(colour="black", size=rel(0.9)),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_rect(size=0.4),
            axis.ticks= element_line(size=0.2))+
            geom_rect(data = shade.opt, aes_(x = NULL,y = NULL,
                                                    xmin=~xmin, xmax=~xmax, ymin=~ymin, ymax=~ymax), alpha=0.5, fill="grey10")+
            geom_rect(data = shade.sub, aes_(x = NULL,y = NULL,
                                       xmin=~xmin, xmax=~xmax, ymin=~ymin, ymax=~ymax), alpha=0.5, fill="grey50")+
      geom_text(aes_(x = ~samples, y = ~upper + 0.01, label = ~cum), na.rm = TRUE)

    return(fig)
   }
}
