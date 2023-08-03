library(dplyr)
library(RColorBrewer)
library(ggplot2)

library(extrafont)
# font_import(prompt=F, paths ='/usr/share/fonts/truetype/myfonts/') # import Helvetica font
fonts()
my_font <- "Helvetica"


# dat: data.frame with columns: (id, s)
# clone_dic: a data.frame with columns: (id, K, letters, pretty_names, colour)
plot_box_posterior_v2 <- function(dat, color_fn, tag, plottitle, ymax=NULL, ymin=NULL,  
                               clone_axis_dir = 'bottom', sort_by = 'median', add_linebreak = FALSE) {
  color_df <- data.table::fread(color_fn)
  # color_df <- data.table::fread(paste0(output_dir,'colorcode_total.csv'))
  color_df <- color_df %>%
    dplyr::filter(datatag==tag) %>%
    dplyr::select(-datatag)
  # clone_dic$K <- as.character(clone_dic$K)
  
  ## Loading clone_dic here
  clone_dic <- dplyr::left_join(dat, color_df, by=c('letters'='clone_id')) %>%
               dplyr::group_by(letters, colour) %>%
               dplyr::summarise(nb=n())
  # dat <- dat[!is.na(dat$id), ]
  print(clone_dic)
  
  # Sort clones by their posterior mean fitness values
  if (sort_by == 'median') {
    s_means <- dat %>% dplyr::group_by(letters) %>% dplyr::summarise(mean_s = median(s)) %>% dplyr::arrange(mean_s)
  } else {
    s_means <- dat %>% dplyr::group_by(letters) %>% dplyr::summarise(mean_s = mean(s)) %>% dplyr::arrange(mean_s)
  }
  
  dat$letters <- factor(x = dat$letters, levels = s_means$letters, ordered = T)
  
  
  myColors <- clone_dic$colour
  names(myColors) <- clone_dic$letters
  # clones_ls <- unique(dat$letters)
  # myColors <- brewer.pal(length(clones_ls), "Set3")
  # names(myColors) <- clones_ls
  
  # dat$pretty_names <- factor(dat$pretty_names, levels = names(myColors), ordered = T)
  # stopifnot(all(levels(dat$pretty_names) == as.character(names(myColors))))
  
  # change to s + 1
  dat$s <- dat$s + 1
  
  # Set the relevant s-values, length of the whiskers  
  # Put the min and max in range of the whiskers...
  whisker_coef <- 1.5
  
  # lims <- dat %>% dplyr::group_by(K) %>% dplyr::summarise(iqr = IQR(s), median = median(s), q1 = quantile(s, .25), q3 = quantile(s, .75), bottom = min(s), top = max(s)) 
  lims <- dat %>% dplyr::summarise(iqr = IQR(s), median = median(s), q1 = quantile(s, .25), q3 = quantile(s, .75), bottom = min(s), top = max(s)) 
  if(is.null(ymin)){
    ymin <- min(lims$bottom)
  }
  if(is.null(ymax)){
    ymax <- max(lims$top)
  }
  p <- 
    dat %>% 
    ggplot(aes(letters, s, color = letters, fill = letters)) +
    geom_boxplot(outlier.shape = NA, coef = whisker_coef) + 
    stat_boxplot(geom = 'errorbar', width = 0.3, size = .5, coef = whisker_coef) + 
    scale_x_discrete(limits = levels(dat$letters), breaks = levels(dat$letters), labels = levels(dat$letters), position = clone_axis_dir) + 
    coord_flip(ylim = c(ymin, ymax)) + #c(min(lims$bottom), max(lims$top))
    geom_hline(yintercept = 1.0, linetype = 'longdash', size = .5) +
    scale_color_manual(values = myColors) + 
    scale_fill_manual(values = myColors) + 
    # xlab('') + 
    # ylab(sprintf('1 + s (relative to clone %s)', ref_clone_letter)) +
    # ylab('1 + s ') +
    labs(x='Clones',y='1 + s', title=plottitle) + 
    # theme_light(base_size = 20) + 
    theme_bw()+
    theme(axis.text.x = element_text(color="black",size=10, hjust = 0.5, family=my_font),
          # axis.text.x = element_blank(),
          axis.text.y = element_text(color="black",size=11, hjust = 0.5, family=my_font, face = 'bold'),
          axis.title.x = element_text(color="black",size=11, hjust = 0.5, family=my_font),
          axis.title.y = element_text(color="black",size=11, hjust = 0.5, family=my_font),
          plot.title = element_text(color="black",size=12, hjust = 0, family=my_font, face = 'bold'),
          # panel.grid.major = element_blank(), 
          strip.background = element_rect(fill = 'white', colour = 'white'),
          panel.grid.minor = element_blank(),
          plot.margin = unit(c(0, 0, 0, 0), "null"),
          panel.spacing = unit(c(0, 0, 0, 0), "null"),
          axis.ticks.y = element_blank(),
          legend.position = "none") #+#, angle = 90) +
    # fitclone_get_theme_no_grid() + 
    # theme()
  
  if (add_linebreak) {
    p <- p + ylab(sprintf('1 + s\n(rel. to clone %s)', ref_clone_letter))
  }
  
  # p
  # if (show_title) {
  #   p <- p + ggtitle(paste0(datatag), subtitle = get_ID_from_edge_list_path(edge_list_path))
  # } else {
  #   p <- p + theme(legend.position = "none")
  # }
  # p <- p + theme(legend.position = "none")
  data <- ggplot_build(p)$data[[1]]
  p <- p + geom_segment(inherit.aes = FALSE, data=data, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", size=.5)
  
  return(p)
}

## Original version from Sohrab
# dat: data.frame with columns: (id, s)
# clone_dic: a data.frame with columns: (id, K, letters, pretty_names, colour)
plot_box_posterior_original_version <- function(dat, clone_dic, show_title=TRUE, ref_clone_letter = '', clone_axis_dir = 'bottom', sort_by = 'median', add_linebreak = FALSE) {
  
  clone_dic$K <- as.character(clone_dic$K)
  dat <- dplyr::right_join(dat, clone_dic)
  dat <- dat[!is.na(dat$id), ]
  
  
  # Sort clones by their posterior mean fitness values
  if (sort_by == 'median') {
    s_means <- dat %>% dplyr::group_by(letters) %>% dplyr::summarise(mean_s = median(s)) %>% dplyr::arrange(mean_s)
  } else {
    s_means <- dat %>% dplyr::group_by(letters) %>% dplyr::summarise(mean_s = mean(s)) %>% dplyr::arrange(mean_s)
  }
  
  dat$letters <- factor(x = dat$letters, levels = s_means$letters, ordered = T)
  
  
  myColors <- clone_dic$colour
  names(myColors) <- clone_dic$letters
  
  dat$pretty_names <- factor(dat$pretty_names, levels = names(myColors), ordered = T)
  stopifnot(all(levels(dat$pretty_names) == as.character(names(myColors))))
  
  # change to s + 1
  dat$s <- dat$s + 1
  
  # Set the relevant s-values, length of the whiskers  
  # Put the min and max in range of the whiskers...
  whisker_coef <- 1.5
  
  lims <- dat %>% dplyr::group_by(K) %>% dplyr::summarise(iqr = IQR(s), median = median(s), q1 = quantile(s, .25), q3 = quantile(s, .75), bottom = min(s), top = max(s)) 
  
  p <- 
    dat %>% 
    ggplot(aes(letters, s, color = letters, fill = letters)) +
    geom_boxplot(outlier.shape = NA, coef = whisker_coef) + 
    stat_boxplot(geom = 'errorbar', width = 0.3, size = .5, coef = whisker_coef) + 
    scale_x_discrete(limits = levels(dat$letters), breaks = levels(dat$letters), labels = levels(dat$letters), position = clone_axis_dir) + 
    coord_flip(ylim = c(min(lims$bottom), max(lims$top))) + 
    geom_hline(yintercept = 1.0, linetype = 'longdash', size = .5) +
    scale_color_manual(values = myColors) + 
    scale_fill_manual(values = myColors) + 
    xlab('') + ylab(sprintf('1 + s (relative to clone %s)', ref_clone_letter)) +
    theme_light(base_size = 20) + 
    fitclone_get_theme_no_grid() + 
    theme(axis.ticks.y = element_blank())
  
  if (add_linebreak) {
    p <- p + ylab(sprintf('1 + s\n(rel. to clone %s)', ref_clone_letter))
  }
  
  
  if (show_title) {
    p <- p + ggtitle(paste0(datatag), subtitle = get_ID_from_edge_list_path(edge_list_path))
  } else {
    p <- p + theme(legend.position = "none")
  }
  
  data <- ggplot_build(p)$data[[1]]
  p <- p + geom_segment(inherit.aes = FALSE, data=data, aes(x=xmin, xend=xmax, y=middle, yend=middle), color="white", size=.5)
  
  return(p)
}