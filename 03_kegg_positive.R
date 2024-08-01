## KEGG functional enrichment on over expressed metabolic markers

list.files()


res<-read.table("limmasig278.tsv",h=T,sep="\t",row.names=1)


library(dplyr)

res%>%filter(logFC>=1)->pos

vector<-row.names(pos)

library(clusterProfiler)
library(org.Hs.eg.db)
library(pathview)
library(DOSE)
library(enrichplot)
gene_list <- bitr(vector, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
entrez_ids <- gene_list$ENTREZID


kegg_enrich <- enrichKEGG(gene = entrez_ids, 
                          organism = 'hsa', 
                          pvalueCutoff = 0.05)

# Convert results back to gene symbols
kegg_result <- as.data.frame(kegg_enrich)
kegg_result$geneID <- sapply(kegg_result$geneID, function(x) {
    genes <- unlist(strsplit(x, "/"))
    symbols <- gene_list$SYMBOL[match(genes, gene_list$ENTREZID)]
    paste(symbols, collapse = "/")
})


kegg_result%>%filter(category=="Metabolism")->kegg_result

kegg_enrich <- pairwise_termsim(kegg_enrich)

dotplot(kegg_enrich) + ggtitle("KEGG Pathway Enrichment")
emapplot(kegg_enrich)
barplot(kegg_enrich, showCategory = 15)



# Pathway visualization
pathview(gene.data = entrez_ids, 
         pathway.id = "hsa00650",  # Example pathway ID
         species = "hsa")

kegg<-as.data.frame(kegg_result)
library(writexl)
write_xlsx(kegg,"kegg_activated_hb.xlsx")

