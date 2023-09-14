
library(biomaRt)

mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")

searchAttributes(mart = ensembl, pattern = "hgnc")
table=read_excel("transcripts_up.xlsx")
table
transcript_ids <- table$...1

res <- getBM(attributes = c('ensembl_transcript_id_version',
                            'external_gene_name',
                            'entrezgene_id'),
             filters = 'ensembl_transcript_id_version',
             values = transcript_ids,
             mart = mart)

c <- merge(table, res, by.x = "...1", by.y = "ensembl_transcript_id_version")
d <- c[order(c$padj),]
e <- d[!duplicated(d$entrezgene_id),,drop=FALSE]
f <- e[order(e$entrezgene_id),]
g <- f[!duplicated(f$external_gene_name),,drop=FALSE]

write.xlsx(as.data.frame(g),"transcripts_converted_up.xlsx", row.names=TRUE)