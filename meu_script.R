# DEG FILTERING ----------------------------------------------------------------

# project organization:
# folders: codes, data, and outputs 

# import data
original_name <- "./data/GSE25001_Whole/GSE25001.100_4to8d_vs_34_convalescent.tsv"
table <- read.table(file = original_name, header = TRUE, sep = "\t", dec = ".", quote = "")

# subset by adjusted p-value < 0.05
p_filtered <- subset(x = table, table$adj.P.Val < 0.05)
# subset by log FC < -1 OR > 1
p_FC_filtered <- subset(x = p_filtered, 
                        subset = (p_filtered$logFC < -1 | p_filtered$logFC > 1))
# remove duplicates by gene symbol, after filtering
DEGs <- p_FC_filtered[!duplicated(p_FC_filtered$Gene.symbol), ]

# remove empty values
if (any(nchar(DEGs$Gene.symbol) == 0) == TRUE){
  to_remove <- which(nchar(DEGs$Gene.symbol) == 0)
  DEGs <- DEGs[-to_remove,]}
if (any(DEGs$Gene.symbol == " ") == TRUE){
  to_remove <- which(DEGs$Gene.symbol == " ")
  DEGs <- DEGs[-to_remove,]}

# subset upregulated and downregulated genes
upreg <- subset(x = DEGs, subset = (DEGs$logFC > 1))
downreg <- subset(x = DEGs, subset = (DEGs$logFC < -1))

## Export results
# Export tables as .tsv:
# all filtered DEGs
write.table(DEGs, file = (paste0(out_path, "/DEG_", new_folder, ".tsv")), 
            sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = TRUE)
# upregulated DEGs
write.table(upreg, file = (paste0(out_path, "/UP_", new_folder, ".tsv")), 
            sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = TRUE)
# downregulated DEGs
write.table(downreg, file = (paste0(out_path, "/DOWN_", new_folder, ".tsv")), 
            sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = TRUE)

# export excel spreadsheet with the all, upregulated and downregulated DEGs as 'tabs'
  write.xlsx(x = DEGs, file = (paste0(out_path, "/excel_", new_folder, ".xlsx")), 
             sheetName = "DEG", row.names = FALSE)
  write.xlsx(x = upreg, file = (paste0(out_path, "/excel_", new_folder, ".xlsx")), 
             sheetName = "UP", row.names = FALSE, append = TRUE)
  write.xlsx(x = downreg, file = (paste0(out_path, "/excel_", new_folder, ".xlsx")), 
             sheetName = "DOWN", row.names = FALSE, append = TRUE)
