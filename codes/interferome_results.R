## RESULTS SUMMARY TABLE -------------------------------------------------------

# Summarizes the data downloaded from the interferome
# from the data search results (experimental data), as it is more accurate.

# list paths
outputs1 <- list.dirs(path = "./output", full.names = T, recursive = F)
GSEs <- grep(x = basename(outputs1), pattern = "^GSE")
outputs2 <- outputs1[GSEs]
folders <- list.dirs(path = outputs2, full.names = T, recursive = F)

# load the comparisons table
COMPARISONS <- read.table(file = "./output/tabela_comparisons.txt", 
                          header = T, sep = "\t", dec = ".", quote = "")
comparisons_table <- COMPARISONS

# list all 'data search results' files
DataSearchRes <- list.files(path = paste0("./output"), pattern = 
                              ".+_DataSearchResults\\.txt", full.names = T, recursive = T)

# get number of unique genes from 'Data search results'
DataN <- lapply(DataSearchRes, function(i){
  table <- read.table(file = i, header = TRUE, sep = "\t", dec = ".", quote = "", skip = 19)
  table1 <- table[!duplicated(table$Gene.Name), ]
  nrow(table1)})
# add file names and manipulate
names(DataN) <- gsub(x = DataSearchRes, pattern = "/[[:alnum:]]+_[[:alnum:]]+\\.txt|\\./output/", replacement = "") 
DataNdf <- melt(DataN)
DataNdf <- dplyr::mutate(DataNdf, folder1 = (gsub(x = DataNdf$L1, pattern = "/.*", "")))
DataNdf <- dplyr::mutate(DataNdf, temp1 = (gsub(x = DataNdf$L1, pattern = "^[[:alnum:]_]+/", "")))
DataNdf <- dplyr::mutate(DataNdf, folder2 = (gsub(x = DataNdf$temp1, pattern = "/.*", "")))
DataNdf <- dplyr::mutate(DataNdf, UPDOWN = (gsub(x = DataNdf$L1, pattern = ".*_", "")))
exclude <- grep(colnames(DataNdf), pattern = "temp.")
DataNdf <- DataNdf[,-exclude]


# get number of UP genes and merge tables
UPDataN <- subset(x = DataNdf, DataNdf$UPDOWN == "UP")
comparisons_table <- merge(x = comparisons_table, y = UPDataN[,c(1,4)], by = "folder2", all.x = TRUE)
colnames(comparisons_table)[colnames(comparisons_table)== "value"] <- "UP_interferome"
# get number of DOWN genes and merge tables
DOWNDataN <- subset(x = DataNdf, DataNdf$UPDOWN == "DOWN")
comparisons_table <- merge(x = comparisons_table, y = DOWNDataN[,c(1,4)], by = "folder2", all.x = TRUE)
colnames(comparisons_table)[colnames(comparisons_table)== "value"] <- "DOWN_interferome"

### number of unique genes by IFN type from 'data search results"

# for each IFN type,
for (IFN_type in c("I","II","III")) {
  # for each table
  # subset by select IFN type; remove duplicates; get number of genes;
  DataNType <- lapply(DataSearchRes, function(i){
    table <- read.table(file = i, header = TRUE, sep = "\t", dec = ".", quote = "", skip = 19)
    table1 <- subset(x = table, table$Inteferome.Type == IFN_type)
    table2 <- table1[!duplicated(table1$Gene.Name), ]
    nrow(table2)})
  # turn to df, manipulate names
  names(DataNType) <- gsub(x = dirname(DataSearchRes), pattern = "\\./output/", replacement = "") 
  DataNTypedf <- melt(DataNType)
  DataNTypedf <- dplyr::mutate(DataNTypedf, folder1 = (gsub(x = DataNTypedf$L1, pattern = "/.*", "")))
  DataNTypedf <- dplyr::mutate(DataNTypedf, temp1 = (gsub(x = DataNTypedf$L1, pattern = "^[[:alnum:]_]+/", "")))
  DataNTypedf <- dplyr::mutate(DataNTypedf, folder2 = (gsub(x = DataNTypedf$temp1, pattern = "/.*", "")))
  DataNTypedf <- dplyr::mutate(DataNTypedf, UPDOWN = (gsub(x = DataNTypedf$L1, pattern = ".*_", "")))
  exclude <- grep(colnames(DataNTypedf), pattern = "temp.")
  DataNTypedf <- DataNTypedf[,-exclude]
  
  # get number of UP unique genes by selected IFN type
  UPDataNType <- subset(x= DataNTypedf, DataNTypedf$UPDOWN == "UP")
  comparisons_table <- merge(x=comparisons_table, y =UPDataNType[,c(1,4)], by = "folder2", all.x = TRUE)
  colnames(comparisons_table)[colnames(comparisons_table)== "value"] <- paste0("UP_Type_",IFN_type)
  
  # get number of DOWN unique genes by selected IFN type
  DOWNDataNType <- subset(x= DataNTypedf, DataNTypedf$UPDOWN == "DOWN")
  comparisons_table <- merge(x=comparisons_table, y =DOWNDataNType[,c(1,4)], by = "folder2", all.x = TRUE)
  colnames(comparisons_table)[colnames(comparisons_table)== "value"] <- paste0("DOWN_Type_",IFN_type)
}

# arrange columns order
# turn NAs to 0
comparisons_table <- comparisons_table %>% dplyr::relocate(DEG, UP, .after = g2sev)
comparisons_table <- comparisons_table %>% dplyr::relocate(DOWN, DOWN_interferome, 
                                                           DOWN_Type_I, DOWN_Type_II ,DOWN_Type_III, .after = UP_Type_III)
comparisons_table <- comparisons_table %>% dplyr::relocate(folder2, .after = studytitle)
comparisons_table[is.na(comparisons_table)] <- 0

# write table with interferome results summarized
write.table(comparisons_table, file = "./output/tabela_interferome_datasearchres.tsv", 
            sep="\t", dec=".", quote = FALSE, row.names = FALSE, col.names = TRUE)
