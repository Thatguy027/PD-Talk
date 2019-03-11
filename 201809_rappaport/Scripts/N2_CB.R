library(tidyverse)
library(extrafont)

try(setwd(dirname(rstudioapi::getActiveDocumentContext()$path)))

# import fonts
#font_import()

# strain colors
strain_colors <- c("N2" = "#F9A227", "CB4856" = "#2790F9")

# color blind friendly colors
# orange, sky blue, bluish green, yellow, blue, vermillion, reddish purple, gray
# from - https://serialmentor.com/dataviz/color-pitfalls.html
col_blind_colors <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999")

# colors
axis_color <- "#000F08"
highlight_color <- "#D7263D"
background_color <- "#FCFFFC"

# font
number_font <- "Itim"
axes_text_size <- 20
axes_title_font <- "Montserrat ExtraBold"
axes_title_size <- 20
title_size <- 20

base_theme <- theme(
  line = element_line(colour = axis_color, size = 0.5, linetype = 1), 
  rect = element_rect(fill = background_color, colour = axis_color, size = 0.5, linetype = 1), 
  text = element_text(family = axes_title_font, size = axes_text_size), 
  
  axis.text = element_text(family = number_font,size = rel(0.8), colour = "grey30", margin = unit(0.1, "cm")),
  axis.text.x = element_text(vjust = 1), 
  axis.text.y = element_text(hjust = 1), 
  axis.ticks = element_line(colour = "gray90"), 
  axis.title.x = element_text(margin = margin(t = 10, r = 0, b = 0, l = 0)), 
  axis.title.y = element_text(margin = margin(t = 0, r = 10, b = 0, l = 0), angle = 90), 
  axis.ticks.length = unit(0.15, "cm"),
  
  strip.text = element_text(size = rel(0.8)), 
  strip.background = element_rect(fill = background_color, colour = NA, size = 0.5, linetype = 1),
  
  plot.background = element_rect(fill = background_color, color = NA),
  
  legend.background = element_rect(fill=background_color,colour = background_color), 
  legend.spacing = unit(0.2, "cm"), 
  legend.key = element_rect(fill = background_color, colour = NA), 
  legend.key.size = unit(1, "lines"), 
  legend.key.height = NULL, 
  legend.key.width = NULL, 
  legend.text = element_text(size = rel(0.6)), 
  legend.text.align = NULL, 
  legend.title = element_text(size = rel(0.6), hjust = 0), 
  legend.title.align = NULL, 
  legend.position = "right", 
  legend.direction = NULL, 
  legend.justification = "center", 
  legend.box = NULL, 
  
  panel.background = element_rect(fill = background_color, colour = NA),  
  panel.grid.major = element_line(colour = "gray90"), 
  panel.grid.minor = element_blank(), 
  panel.spacing = unit(1, "lines"), 
  panel.margin.x = NULL, 
  panel.margin.y = NULL)




