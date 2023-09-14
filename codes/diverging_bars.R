# Diverging bars plot
# Reference:
# https://www.r-bloggers.com/2018/05/diverging-bar-charts-plotting-variance-with-ggplot2/
# Written by Julia N. Usuda 
# https://github.com/JNUsuda/

# An R script for a diverging bars (butterfly) plot.
# For only total DEGs (separated by up and downregulation), 
#  use first excel sheet and the first ggplot().  
# For total DEGs and subset DEGs  (separated by up and downregulation), 
#  use second excel sheet and the secong ggplot().
# Adjust settings (plot size, filename) before saving.

library(tidyverse)
library(ggplot2)
library(readxl) # to read from excel
library(reshape) # to melt dataframe
library(svglite) # to export as .svg

# import table
## !!! for only DEGs:
tabela <- read_excel("./tabela_diverging_bars.xlsx", sheet = 1)
## !!! for DEGs and subset:
tabela <- read_excel("./tabela_diverging_bars.xlsx", sheet = 2)

tabela1 <- reshape::melt(as.data.frame(tabela), id.vars = c("Names", "Classification"), variable_name = "Type")
tabela1 <- tidyr::separate(tabela1, Type, sep = "_", into = c("Type", "UPDOWN"))

# set order 
neworder <- rev(unique(tabela$Names))
tabela_final <- dplyr::arrange(transform(tabela1, Names=factor(Names,levels=neworder)), Names)

# !!! plot for only total DEGs:
## skip to bottom to save
ggplot(data = tabela_final) +
  #downregulated total DEGs, and labels
  geom_bar(data = subset(tabela_final, (UPDOWN == "DOWN")), 
           position = position_dodge(), stat = "identity", 
           aes(fill = paste(UPDOWN), y = -value, x = Names),
           alpha = 1, width = 0.9) +
  geom_text(data = subset(tabela_final, (UPDOWN == "DOWN")),
            aes(y = -value, x = Names, label = value, hjust = 1.2),
            size = 3) +
  #upregulated total DEGs, and labels
  geom_bar(data = subset(tabela_final, (UPDOWN == "UP")), 
           position = position_dodge(), stat = "identity", 
           aes(fill = paste(UPDOWN), y = value, x = Names), 
           alpha = 1, width = 0.9)+
  geom_text(data = subset(tabela_final, (UPDOWN == "UP")),
            aes(y = value, x = Names, label = value, hjust = - 0.2),
            size = 3) +
  theme(legend.position = "bottom", 
        legend.title =  element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(color = "black"),
        panel.background = element_rect( fill = "white", color = "grey"),
        panel.grid = element_blank(),
        strip.text.y = element_text(vjust = 1, hjust = 0.5), 
        strip.background = element_rect(colour="grey", fill="gray95")) +
  labs(x = element_blank(), y = element_blank()) +
  coord_flip() +
  facet_grid(rows = vars(Classification), scale = "free_y", space = "free_y") +
  scale_fill_manual(values=c("#999999", "#85bb41"),
                    labels = c("Down-regulated DEGs", "Up-regulated DEGs")) 


# !!! plot por DEGs and DEGs subset:
ggplot(data = tabela_final) +
  #downregulated total DEGs and downregulated subset DEGs, and labels
  geom_bar(data = subset(tabela_final, (UPDOWN == "DOWN")), 
           position = position_dodge(), stat = "identity", 
           aes(fill = paste(UPDOWN, Type), y = -value, x = Names),
           alpha = 1, width = 0.9)+
  geom_text(data = subset(tabela_final, (UPDOWN == "DOWN")),
            aes(y = -value, x = Names, label = value, hjust = 1.1, 
                vjust =  ifelse(test = (Type == "subset"), yes = -0.5, no = 1.3)),
            size = 3) +
  #upregulated total DEGs and upregulated subset DEGs, and labels
  geom_bar(data = subset(tabela_final, (UPDOWN == "UP")), 
           position = position_dodge(), stat = "identity", 
           aes(fill = paste(UPDOWN, Type), y = value, x = Names), 
           alpha = 1, width = 0.9)+
  geom_text(data = subset(tabela_final, (UPDOWN == "UP")),
            aes(y = value, x = Names, label = value, hjust = -0.1, 
                vjust =  ifelse(test = (Type == "subset"), yes = -0.5, no = 1.3)),
            size = 3) +
  theme(legend.position = "bottom", 
        legend.title =  element_blank(),
        axis.ticks = element_blank(),
        axis.text.x = element_blank(), 
        axis.text.y = element_text(color = "black"),
        panel.background = element_rect( fill = "white", color = "grey"),
        panel.grid = element_blank(),
        strip.text.y = element_text(vjust = 1, hjust = 0.5), 
        strip.background = element_rect(colour="grey", fill="gray95")) +
  labs(x = element_blank(), y = element_blank()) +
  coord_flip() +
  facet_grid(rows = vars(Classification), scale = "free_y", space = "free_y") +
  scale_fill_manual(values=c("#cccccc","#999999", "#bddb98", "#85bb41"),
                    labels = c("Down-regulated DEGs", "Down-regulated subset DEGs", 
                               "Up-regulated DEGs", "Up-regulated subset DEGs")) +
  guides(fill = guide_legend(nrow = 2, byrow = F)) 


# !!! set width and height before exporting 
# !!! WARNING: overwrites files with same filename
#png:
ggsave(filename = "bars_DEGs.png", path = "./", width = 2900, height = 2200, 
       units = "px", dpi = 300, device = "png", scale = 1)
#svg:
ggsave(filename = "bars_DEGs.svg", path = "./", width = 2900, height = 2200, 
       units = "px", device = "svg", scale = 1)
