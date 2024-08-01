## prepare single cell RNA-seq object

library(Seurat)
library(hdf5r)
library(zellkonverter)
library(SingleCellExperiment)


# Read the .h5ad file into an AnnData object
sce <- readH5AD("GSE180665_hb_integrated_normalized_annotated_harmony.h5ad")

names(assays(sce))=c("counts")

assay(sce,"logcounts")<-assay(sce,"counts")

data<-as.seurat(sce)

save(data,file="hepatoblastoma.rda")
load("hepatoblastoma.rda")

str(data)
x<-data[[]]

## Seurat read and add original meta data annotation
meta<-read.csv("meta.csv",h=T,row.names=1)
all(row.names(meta)==row.names(x))
data <- AddMetaData(object = data, metadata = meta)

## unintegrate seurat pipeline

data<-FindVariableFeatures(data)
data<-ScaleData(data)
data<-RunPCA(data)
ElbowPlot(data,ndims = 50)
data<-RunUMAP(data,dims=1:30)


## set Seurat identities
Idents(object = data) <- 'Cell.class_concise'
Idents(object = data)

library(pals)
DimPlot(data,reduction="umap",group.by ="sample",pt.size=1,cols=cols25(),label=T)



### integration
data[["originalexp"]] <- split(data[["originalexp"]], f = data$sample)

data <- IntegrateLayers(
  object = data, method = HarmonyIntegration,
  orig.reduction = "pca", new.reduction = "harmony",
  verbose = FALSE
)


DimPlot(data, reduction = "harmony",group.by="sample",cols=cols25(),label=F)
DimPlot(data, reduction = "harmony",group.by="Cell.class_concise",cols=cols25(),label=F)
save(all,file="harmony.rda")



all$class<-x$class
DimPlot(data, reduction = "umap",group.by="Cell.class_concise")


data <- JoinLayers(data)
data

data<-RunPCA(data)
ElbowPlot(data,ndims = 50)


## run umap post harmony integration
data<-RunUMAP(data,dims=1:30, reduction = "harmony")

save(data,file="all.rda")



DimPlot(data,reduction="umap",group.by ="sample_group",pt.size=0.5,cols=cols25(),label=F)

DimPlot(data,reduction="umap",group.by ="RNA_snn_res.0.5",pt.size=0.5,cols=cols25(),label=T)
save(data,file="harmony.rda")

