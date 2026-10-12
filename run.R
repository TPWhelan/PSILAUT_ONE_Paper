library(reshape2)
library(tidyverse)
library(shadowtext)


source('parseData.R')

#insert name of data folder below 

filelist <- dir(path = 'data', 
                full.names = TRUE) %>% 
  as_tibble()

names <- read.csv(file = 'schaefer_100_labels.csv', header = FALSE, fileEncoding = 'UTF-8-BOM')

filelist %>% 
  dplyr::mutate(data = purrr::map(value, 
                                  ~parseData(.x, names[,1]))) %>% 
  tidyr::unnest(data) %>% 
  dplyr::group_by(Var1, Var2) %>% 
  dplyr::summarise(bigMean = mean(mean_z)) %>% 
  
  
  ggplot2::ggplot(aes(x = as.factor(Var1), 
                      y = as.factor(Var2))) + 
  geom_tile(aes(fill = bigMean)) + 
  geom_text(aes(label = format(round(bigMean, 
                                           digits = 2), 
                                     nsmall = 2)), 
                  colour = I("black")) + 
  scale_fill_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 11,
                                                              name = 'RdBu')), 
                       limit = c(-1, 1), 
                       na.value = 'transparent', 
                       name = 'Mean Z-score') + 
  coord_flip() + 
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 0),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.4, 0.7),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))


