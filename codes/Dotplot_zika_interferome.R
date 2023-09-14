install.packages("tidyverse")
library(tidyverse)
library(viridis)
library(ggplot2)
library(readxl)
library(xlsx)

outputs1 <- list.files(path = "./data/enrichr_results/", full.names = T, recursive = F)
GO_BPs <- grep(x = basename(outputs1), pattern = "^GO")
outputs_BPs <- outputs1[GO_BPs]
#folders <- list.dirs(path = outputs2, full.names = T, recursive = F)
react <- grep(x = basename(outputs1), pattern = "^React")
outputs_react <- outputs1[react]

# export tables in the list as sheets in a xlsx
lapply(outputs1, function(i){
  tabela <- read.table(file = i, header = TRUE, sep="\t", dec=".", quote = "")
  exclude <- grep(x = colnames(tabela), pattern = "Old")
  tabela2 <- tabela[-(exclude)]
  tabela2 <- tabela2 %>% dplyr::mutate(TermID = (gsub(x = Term, pattern = ".*\\(|\\)", "")))
  tabela2 <- tabela2 %>% dplyr::relocate(TermID, .after = Term)
  tabela3 <- subset(tabela2, tabela2$Adjusted.P.value < 0.05)
  xlname <- gsub(x = i, pattern = ".*/|\\.txt", replacement = "")
  write.xlsx(x = tabela3, file = "./data/enrichr_results/Enrichr_filtered.xlsx", 
             sheetName = xlname, row.names = F, append = TRUE)
})


tabela <- read.delim(file = outputs_BPs[1])



#tabela <- read_xlsx(path = "data/dotplot_example.xlsx")
tabela$p.adjust <- as.numeric(tabela$p.adjust)
df <- tabela

df$Description <- factor(df$Description, levels = df$Description[order(df$GeneRatio)])


# create the plot
df <- df %>% filter(p.adjust < 0.05)
df <- df[1:30,]

ggplot(df, aes(x = GeneRatio, y = Description, color=p.adjust)) +
  geom_point(mapping = aes(size = p.adjust)) +
 # scale_x_continuous(limits = c(0, 70)) +
  theme_bw(14) +
  theme(axis.text.y = element_text(size = rel(1.0)),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size = rel(.75)),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(size = 0.5),
        panel.grid.minor.x = element_blank()) +
  xlab ("Gene Ratio") +
  ylab("GO category") +
  ggtitle("DotPlot by GO")+
  scale_color_gradient(low="green", high="royalblue")

