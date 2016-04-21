fte_theme <- function() {
    ## inspiré de http://minimaxir.com/2015/02/ggplot-tutorial/
    ## Generate the colors for the chart procedurally with RColorBrewer
    palette <- RColorBrewer::brewer.pal("Greys", n=9)
    color.background = palette[2]
    color.grid.major = palette[3]
    color.axis.text = palette[6]
    color.axis.title = palette[7]
    color.title = palette[9]
    ## Begin construction of chart
    theme_bw(base_size=9) +
        ## Set the entire chart region to a light gray color
        theme(
            text = element_text(family = "Ubuntu Light"),
            panel.background=element_rect(
                fill=alpha(color.background, 1), color=color.background),
            plot.background=element_rect(
                fill=alpha(color.background, 0.5), color=color.background),
            panel.border=element_rect(color=color.background),
            panel.grid.major=element_line(
                color=color.grid.major,size=.25, linetype = "dotted"),
            panel.grid.minor=element_blank(),
            axis.ticks=element_blank(),
            strip.text.y = element_text(angle = 0, size = 8),
            legend.position="none", # Format the legend, but hide by default
            legend.background = element_blank(),
            legend.text = element_text(size=9,color=color.axis.title),
            legend.key = element_rect(
                fill = scales::alpha("gray", 0),
                colour = scales::alpha("gray", 0)),
            ## Set title and axis labels, and format these and tick marks
            plot.title = element_text(
                color = color.title, size = 14, face = "bold", vjust = 1.25,
                hjust = 0, family = "Ubuntu"),
            axis.text.x = element_text(size = 7,color = color.axis.text),
            axis.text.y = element_text(size = 7,color = color.axis.text),
            axis.title.x = element_text(
                size = 8, color = color.axis.title, vjust = 0, hjust = 0.9),
            axis.title.y = element_text(
                size = 8, color = color.axis.title, vjust = 0.9, angle = 0 ),
            ## Plot margins
            plot.margin = unit(c(0.35, 0.2, 0.3, 0.35), "cm")
        )
}
#
theme_set(theme_bw() + fte_theme())
#
legend_position <- function(x=NULL, y=NULL) {
    custom_legend <- theme(
        legend.margin = unit(-0.3,"lines"),
        legend.key = element_rect(fill = scales::alpha("gray", 0),
                                  colour = scales::alpha("gray", 0)))
    if (is.null(x) & is.null(y)) {
        custom_legend + theme(legend.position = "bottom")
    } else {
        custom_legend + theme(legend.position = c(x, y))
    }
}
