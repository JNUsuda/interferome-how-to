#http://yulab-smu.top/clusterProfiler-book/index.html

#Import your dataset
#Col1 name as "ID" (Genes ID in Entrez form). Use transgeneID or Biomart for conversion
#Col2 "LogFC"(LogFoldChange)

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("ReactomePA")
BiocManager::install("DOSE")
BiocManager::install("enrichplot")
BiocManager::install("ggnewscale")
install.packages("rlang")

#Access libraries
library(DOSE)
library(viridis)
library(ggnewscale)
library(enrichplot)
library(ReactomePA)
library(clusterProfiler)
library(dplyr)
library(ggplot2)

#COLOR PALLETES WEBSITE: http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually
## assume that 1st column is ID
## 2nd column is fold change

d <- geneList
geneList <- d[,2]
names(geneList) <- as.character(d[,1])
genevector <- names(geneList)
nomes <- d$ID
valores <- d$LogFC
names(valores) <- nomes
geneList <- sort(geneList, decreasing = TRUE)


edo <- enrichPathway(genevector)
edo <- pairwise_termsim(edo, method = "JC", semData = NULL, showCategory = 200)
edo2 <- gsePathway(geneList, pvalueCutoff = 0.2)
egoBP <- enrichGO(genevector, ont="BP",  OrgDb = "org.Hs.eg.db")
egoBP <- pairwise_termsim(egoBP, method = "JC", semData = NULL, showCategory = 200)
egoBP2 <- gseGO(geneList, ont="BP",  OrgDb = "org.Hs.eg.db")
egoCC <- enrichGO(genevector, ont="CC",  OrgDb = "org.Hs.eg.db")
egoCC <- pairwise_termsim(egoCC, method = "JC", semData = NULL, showCategory = 200)
egoCC2 <- gseGO(geneList, ont="CC",  OrgDb = "org.Hs.eg.db")
egoMF <- enrichGO(genevector, ont="MF",  OrgDb = "org.Hs.eg.db")
egoMF <- pairwise_termsim(egoMF, method = "JC", semData = NULL, showCategory = 200)
egoMF2 <- gseGO(geneList, ont="MF",  OrgDb = "org.Hs.eg.db")

edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
egoxCC <- setReadable(egoCC, 'org.Hs.eg.db', 'ENTREZID')
egoxBP <- setReadable(egoBP, 'org.Hs.eg.db', 'ENTREZID')
egoxMF <- setReadable(egoMF, 'org.Hs.eg.db', 'ENTREZID')


###############################################################
## PLOTS (ORA, Pahtways)

#Dot plot of enriched terms.
dotplot(edo, showCategory=15) + ggtitle("dotplot for ORA") + scale_color_viridis()
# to change the color scheme replace scale_color_viridis() for (gradient as: low midle and high): scale_color_gradientn(colours = c("#89DC8C","#34A593", "#0B1199"))

# Bar plot of enriched terms.
edo@result <- edo@result %>% arrange(desc(Count)) 

col_fun = RColorBrewer::brewer.pal(n =9, name = "RdYlBu") 
barplot(edo, showCategory=15, fill="p.adjust")+ 
  ggtitle("barplot for ORA")+ 
  scale_fill_viridis()

edo@result <- edo@result %>% arrange(desc(p.adjust))

#Ridgeline plot - Histogram
ridgeplot(edo2, fill="pvalue")+ 
#fill= "p.adjust" if necessary
xlab("logFC")

#Network plot of enriched terms.
p1 <- cnetplot(edox, foldChange=valores)
p2 <- cnetplot(edox, foldChange=valores, circular = TRUE, colorEdge = TRUE, showCategory = 15) + scale_color_viridis()
cowplot::plot_grid(p2, ncol=1, rel_widths=c(1, 1))

#Labelling nodes by selected subset
p3 <- cnetplot(edox, node_label="none")
p4 <- cnetplot(edox, node_label="gene") 
p5 <- cnetplot(edox, node_label="category") 
p6 <- cnetplot(edox, node_label="all") 
cnetplot(edox, foldChange=valores)+scale_color_viridis(option = "plasma")
cowplot::plot_grid(p6, ncol=1) + xlab("logFC")


#Heatmap plot of enriched terms
p7 <- heatplot(edox)
p8 <- heatplot(edox, foldChange=valores)+ scale_color_viridis()
cowplot::plot_grid(p8, ncol=1)


#Plot for results obtained from hypergeometric test and gene set enrichment analysis.
p9 <- emapplot(edo)
p10 <- emapplot(edo, cex_category=1.5)
p11 <- emapplot(edo,layout="circle")
p12 <- emapplot(edo, cex_category=1.5,layout="circle", showCategory = 15, cex_line = 0.4)+ scale_color_viridis()
cowplot::plot_grid(p12, ncol=1)

#Network using specific pathways
View(edox@result)
#choose your pathways
goid <- c("R-HSA-162909", "R-HSA-162906", "R-HSA-1169410", "R-HSA-877300")
edox@result <- edox@result[goid, ]
cnetplot(edox, foldChange=valores, circular = TRUE, colorEdge = TRUE, showCategory = 6) + scale_color_viridis()

#After running this, you should run the enrichment again
edo <- enrichPathway(genevector)
edo <- pairwise_termsim(edo, method = "JC", semData = NULL, showCategory = 200)
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')

###############################################################

## PLOTS (BP)

#Dot plot of enriched terms.
dotplot(egoBP, showCategory=20) + ggtitle("dotplot for BP")+ scale_color_viridis()
View(egoBP)

# Bar plot of enriched terms.
egoBP@result <- egoBP@result %>% arrange(desc(Count))

barplot(egoBP, showCategory=13, color="p.adjust")+ 
  ggtitle("barplot for BP")+ 
  scale_fill_gradient(low="peachpuff", high="indianred3")

egoBP@result <- egoBP@result %>% arrange(desc(p.adjust))

#Ridgeline plot - Histogram
ridgeplot(egoBP2, fill="pvalue")+
#fill= "p.adjust" if necessary
  xlab("logFC")
  

#Network plot of enriched terms.
p1 <- cnetplot(egoxBP, foldChange=valores)
p2 <- cnetplot(egoxBP, foldChange=valores, circular = TRUE, colorEdge = TRUE, showCategory = 5)+ scale_color_viridis()
cowplot::plot_grid(p2, ncol=1, rel_widths=c(1, 1))

#Labelling nodes by selected subset
p3 <- cnetplot(egoxBP, node_label="none")
p4 <- cnetplot(egoxBP, node_label="gene") 
p5 <- cnetplot(egoxBP, node_label="category") 
p6 <- cnetplot(egoxBP, node_label="all") 
cowplot::plot_grid(p6, ncol=1)

#Heatmap plot of enriched terms
p7 <- heatplot(egoxBP)
p8 <- heatplot(egoxBP, foldChange=valores)+ scale_color_viridis()
cowplot::plot_grid(p8, ncol=1)

#Plot for results obtained from hypergeometric test and gene set enrichment analysis.
p9 <- emapplot(egoBP)
p10 <- emapplot(egoBP, cex_category=1.5)
p11 <- emapplot(egoBP,layout="kk")
p12 <- emapplot(egoBP, cex_category=1.5,layout="kk")+ scale_color_viridis() 
cowplot::plot_grid(p12, ncol=1)

#Network using specific pathways.Change according the refs below
View(egoxBP@result)

#choose your pathways. change according the refs below
goidBP <- c("GO:0006260","GO:0140014","GO:0071103","GO:0007059", "GO:0000075", "GO:0000082","GO:0000070","GO:0006302","GO:0006310")
egoxBP@result <- egoxBP@result[goidBP, ]
cnetplot(egoxBP, foldChange=valores, circular = TRUE, colorEdge = TRUE)+ scale_color_viridis()
emapplot(egoxBP, cex_category=1.5,layout="kk")+ scale_color_viridis() 
cnetplot(egoxBP, foldChange=valores, colorEdge = TRUE)+ scale_color_viridis() 
  
  
barplot(egoxBP, showCategory=15, color="p.adjust")+ 
  ggtitle("barplot for BP")+ 
  scale_color_viridis()
dotplot(egoxBP, showCategory=20) + ggtitle("dotplot for BP")+ scale_color_viridis()
#Ridgeline plot - Histogram
ridgeplot(egoBP2, fill="pvalue")+
  #fill= "p.adjust" if necessary
  xlab("logFC")
heatplot(egoxBP, foldChange=valores)+ scale_color_viridis()


#After running this, you should run the enrichment again
egoBP <- enrichGO(genevector, ont="BP",  OrgDb = "org.Hs.eg.db")
egoBP <- pairwise_termsim(egoBP, method = "JC", semData = NULL, showCategory = 200)
egoxBP <- setReadable(egoBP, 'org.Hs.eg.db', 'ENTREZID')

###########################################################
## PLOTS (CC)

#Dot plot of enriched terms.
dotplot(egoCC, showCategory=20) + ggtitle("dotplot for CC") + scale_color_viridis()

# Bar plot of enriched terms.
egoCC@result <- egoCC@result %>% arrange(desc(Count))

barplot(egoCC, showCategory=15, color="p.adjust")+ 
  ggtitle("barplot for CC")+ 
  scale_color_viridis()

egoCC@result <- egoCC@result %>% arrange(desc(p.adjust))

#Ridgeline plot - Histogram
ridgeplot(egoCC2, fill="pvalue")+ 
#fill= "p.adjust" if necessary
xlab("logFC")

#Network plot of enriched terms.
p1 <- cnetplot(egoxCC, foldChange=valores)
p2 <- cnetplot(egoxCC, foldChange=valores, circular = TRUE, colorEdge = TRUE, showCategory = 5) + scale_color_viridis()
cowplot::plot_grid(p2, ncol=1, rel_widths=c(1, 1))

#Labelling nodes by selected subset
p3 <- cnetplot(egoxCC, node_label="none")
p4 <- cnetplot(egoxCC, node_label="gene") 
p5 <- cnetplot(egoxCC, node_label="category") 
p6 <- cnetplot(egoxCC, node_label="all") 
cowplot::plot_grid(p6, ncol=1)

#Heatmap plot of enriched terms
p7 <- heatplot(egoxCC)
p8 <- heatplot(egoxCC, foldChange=valores)+ scale_color_viridis()
cowplot::plot_grid(p8, ncol=1)

#Plot for results obtained from hypergeometric test and gene set enrichment analysis.
p9 <- emapplot(egoCC)
p10 <- emapplot(egoCC, cex_category=1.5)
p11 <- emapplot(egoCC,layout="kk")
p12 <- emapplot(egoCC, cex_category=1.5,layout="kk")  + scale_color_viridis()
cowplot::plot_grid(p12, ncol=1)


#Network using specific pathways
View(egoxCC@result)
#choose your pathways
goidCC <- c("GO:1990204", "GO:0098563", "GO:0030665", "GO:0044815", "GO:0031985")
egoxCC@result <- egoxCC@result[goidCC, ]
cnetplot(egoxCC, foldChange=valores, circular = TRUE, colorEdge = TRUE) + scale_color_viridis()

#After running this, you should run the enrichment again
egoCC <- enrichGO(genevector, ont="CC",  OrgDb = "org.Hs.eg.db")
egoCC <- pairwise_termsim(egoCC, method = "JC", semData = NULL, showCategory = 200)
egoxCC <- setReadable(egoCC, 'org.Hs.eg.db', 'ENTREZID')


###############################################################

## PLOTS (MF)

#Dot plot of enriched terms.
dotplot(egoMF, showCategory=20) + ggtitle("dotplot for MF")+ scale_color_viridis()

# Bar plot of enriched terms.
#x - counts
egoMF@result <- egoMF@result %>% arrange(desc(Count))

barplot(egoMF, showCategory=15, color="p.adjust")+ 
  ggtitle("barplot for MF")+ 
  scale_color_viridis()

egoMF@result <- egoMF@result %>% arrange(desc(p.adjust))

#Ridgeline plot - Histogram
ridgeplot(egoMF2, fill="pvalue")+ 
#fill= "p.adjust" if necessary
xlab("logFC")

#Network plot of enriched terms.
p1 <- cnetplot(egoxMF, foldChange=valores)
p2 <- cnetplot(egoxMF, foldChange=valores, circular = TRUE, colorEdge = TRUE, showCategory = 5)+ scale_color_viridis()
cowplot::plot_grid(p2, ncol=1, rel_widths=c(1, 1))

#Labelling nodes by selected subset
p3 <- cnetplot(egoxMF, node_label="none")
p4 <- cnetplot(egoxMF, node_label="gene") 
p5 <- cnetplot(egoxMF, node_label="category") 
p6 <- cnetplot(egoxMF, node_label="all") 
cowplot::plot_grid(p6, ncol=1)

#Heatmap plot of enriched terms
p7 <- heatplot(egoxMF)
p8 <- heatplot(egoxMF, foldChange=valores)+ scale_color_viridis()
cowplot::plot_grid(p8, ncol=1)

#Plot for results obtained from hypergeometric test and gene set enrichment analysis.
p9 <- emapplot(egoMF)
p10 <- emapplot(egoMF, cex_category=1.5)
p11 <- emapplot(egoMF,layout="kk")
p12 <- emapplot(egoMF, cex_category=1.5,layout="kk")+ scale_color_viridis() 
cowplot::plot_grid(p12, ncol=1)


#Network using specific pathways. change according the refs below
View(egoxMF@result)
#choose your pathways
goidMF <- c("GO:0016887", "GO:0003678", "GO:0003777", "GO:0003688", "GO:0003774")
egoxMF@result <- egoxMF@result[goidMF, ]
cnetplot(egoxMF, foldChange=valores, circular = TRUE, colorEdge = TRUE)+ scale_color_viridis()

#After running this, you should run the enrichment again
egoMF <- enrichGO(genevector, ont="MF",  OrgDb = "org.Hs.eg.db")
egoMF <- pairwise_termsim(egoMF, method = "JC", semData = NULL, showCategory = 200)

egoxMF <- setReadable(egoMF, 'org.Hs.eg.db', 'ENTREZID')

#######################################################################
#http://yulab-smu.top/clusterProfiler-book/chapter4.html#gsencg-fuction
#Density plots for diseases association with the list of DEGs

#gseDO means GSE Diseases Ontology
y <- gseDO(geneList,
           nPerm         = 100,
           minGSSize     = 120,
           pvalueCutoff  = 1,
           verbose       = FALSE)
ridgeplot(y, fill= "pvalue")

#gseDGN means GSE Diseases Ontology

dgn <- gseDGN(geneList,
              nPerm         = 100,
              minGSSize     = 120,
              pvalueCutoff  = 1,
              verbose       = FALSE)

ridgeplot(dgn, fill= "pvalue")

#gseNCG means GS Network of Cancer Gene

ncg <- gseNCG(geneList,
              nPerm         = 100,
              minGSSize     = 120,
              pvalueCutoff  = 1,
              verbose       = FALSE)

ridgeplot(ncg, fill= "pvalue")


#######################################################################
#Export tables

install.packages("writexl")
library(writexl)

tabelaORA <- edox@result
tabelaCC <- egoxCC@result
tabelaBP <- egoxBP@result
tabelaMF <- egoxMF@result


write_xlsx(tabelaORA, "tabelaORA.xlsx")
write_xlsx(tabelaBP, "tabelaBP-upacute.xlsx")
write_xlsx(tabelaCC, "tabelaCC.xlsx")
write_xlsx(tabelaMF, "tabelaMF.xlsx")

