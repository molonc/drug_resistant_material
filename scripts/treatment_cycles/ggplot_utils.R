library(grid)
library(gridExtra)



### To use this theme with ggplot
## 1. source this file in your script
## 2. make the ggplot object
## 3. ggplot_object + thesis_theme
## 4. save as pdf with the final height and width sizes.

thesis_theme <- cowplot::theme_cowplot() +
  theme(
    text = element_text(size = 8),
    axis.title.x = element_text(size=8),  
    axis.title.y = element_text(size=8),
    axis.text.x = element_text(size=7),  
    axis.text.y = element_text(size=7),
    plot.title = element_text(size=10, face="bold", hjust=0),
    strip.text.x = element_text(size=9),
    strip.text.y = element_text(size=9),
    legend.title=element_text(size=7), 
    legend.text=element_text(size=7),
    legend.spacing.x = unit(0.1, 'mm'),
    legend.spacing.y = unit(0.1, 'mm'),
    legend.key.height=unit(1,"line")
  )


### To save the legend in a different file:
## 1. mylegend <- get_legend(ggplot_object) 
## 2. as_ggplot(mylegend)
## 3. ggsave("mylegend.pdf")

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}