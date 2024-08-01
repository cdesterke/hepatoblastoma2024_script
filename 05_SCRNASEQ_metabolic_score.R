## single metabolic score

load("data.rda")


library(Seurat)
library(viridis)
library(RColorBrewer)
library(ggplot2)



## compute metabolic score in seurat 
meta.score <- list(c('GLUL','ATP10A','UGT3A2','UGT3A2','CEL','ACSL4','ATP1A2',
	'OXCT1','NT5DC2','ALDH1L2','SQLE','PDE4C','PRKAA2','DNMT3B','GAL3ST4',
	'FKBP10','CHST10','PYCR1','HSDL1','SOAT2','CKB','BCAT1','SULT1C2','GPAM',
	'SOD3','GPX8','PDE5A','ENO2','ATP1B3','SCD','SCD','HDC','ASNS','PAPSS1',
	'PLA2G7','PKM','HMGCS1','GPX7','GSTP1','GSTP1','PFKFB4','P4HA2','ENTPD1',
	'CHST15','ACLY','CHST3','B3GALT2','FDFT1','MTHFD1L','CA5B','PFKM','ISYNA1',
	'HK2','MTR','NUDT4','MSRB3','DDAH2','LPCAT2','CS'))

data<- AddModuleScore(object = data,
  features = meta.score, name = 'metabolism.score',ctrl=59)

FeaturePlot(data,"metabolism.score1",reduction = "umap",pt.size=0.5,min.cutoff = "q9") + 
  scale_colour_gradientn(colours = rev(brewer.pal(n = 10, name = "Spectral"))) + 
ggtitle("metabolism.score1")

VlnPlot(data, features = c("metabolism.score1"), slot = "data", log = TRUE,pt.size=1,split.by="Cell.class_concise",col=cols25())

length(meta.score[[1]])


## subset tumor cell  and heptocytes

small<-subset(data,idents=c("Tumor Cell","Hep"))
library(pals)
VlnPlot(small, features = c("DNMT3B"), slot = "data", log = TRUE,pt.size=1,split.by="Cell.class_concise",col=cols25())


## export meta data of small object to perform machine learning
x<-small[[]]

x$Cell.class_concise<-as.factor(x$Cell.class_concise)

library(pROC)
library(Epi)

       
ROC(form = Cell.class_concise~metabolism.score1, plot="ROC", data=x)


genes_of_interest<-meta.score[[1]]


count_matrix <- GetAssayData(object = small, slot = "data")

# Ensure the genes of interest are present in the count matrix
genes_in_matrix <- rownames(count_matrix)
missing_genes <- setdiff(genes_of_interest, genes_in_matrix)
if(length(missing_genes) > 0){
  stop(paste("The following genes are not found in the count matrix:", paste(missing_genes, collapse = ", ")))
}

# Extract expression data for the genes of interest
gene_expression_data <- as.data.frame(t(count_matrix[genes_of_interest, ]))

# Add gene expression data to Seurat metadata
for(gene in genes_of_interest){
  small <- AddMetaData(object = small, metadata = gene_expression_data[, gene], col.name = gene)
}

# Check the Seurat object metadata to confirm addition
df<-small@meta.data

write.csv(df,file="metabolism_meta_data.csv",row.names=T)





