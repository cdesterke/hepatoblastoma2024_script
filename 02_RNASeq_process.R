## RNAseq preprocess

list.files()

## load pheno data
pheno<-read.csv("pheno.csv",h=T)

load  counts data
data<-read.table("counts.txt",sep="\t",h=T)

## annotate with biomart 110 dataset
bm110%>%dplyr::rename(Geneid="hs.id")->bm110

bm110%>%inner_join(data,by="Geneid")->df

df$Geneid<-NULL
df$mm.id<-NULL
df$mm.gene<-NULL
df%>%dplyr::rename(gene="hs.gene")->data

## prepare count dataframe with transpipe R-package
## library(devtools)
## devtools::install_github("cdesterke/transpipe")

library(transpipe)
ok<-filtermatrix(data)

## edgeR normalization
library(edgeR)

d0 <- DGEList(ok)

d0<-calcNormFactors(d0,method="TMM")

## edgeR voom transformation
y<-voom(d0, plot=T)
data<-y$E

data<-as.data.frame(data)

save(data,file="voom.rda")

write.csv(data,file="TNMVOOM.csv",row.names=T)


colramp = colorRampPalette(c(3,"darkblue",5))(24)
plot(density(data[,1]),col=colramp[1],lwd=3,ylim=c(0,.2))
	for(i in 1:24){lines(density(data[,i]),lwd=1,col=colramp[i])}

## remove low expressed genes
edata = data[rowMeans(data) > 0, ]
plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.25))
	for(i in 1:24){lines(density(edata[,i]),lwd=1,col=colramp[i])}

save(edata,file="TNMVOOMfiltrated12386.rda")
write.csv(edata,file="TNMVOOMfiltrated12386.csv",row.names=T)


## quantile normalization
library(preprocessCore)


norm = normalize.quantiles(as.matrix(edata))

plot(density(norm[,1]),col=colramp[1],lwd=3,ylim=c(0,.25))
	for(i in 1:24){lines(density(norm[,i]),lwd=1,col=colramp[i])}



row.names(norm)<-row.names(edata)
colnames(norm)<-colnames(edata)

save(norm,file="quantile.rda")


## remove cell line samples from dataset
pheno%>%filter(tissue!="Hepatoblastoma cell line")->pheno

row.names(pheno)<-pheno$id

small<-norm[,row.names(pheno)]


data<-small[row.names(small)%in%vector,]


## DEG with transpipe
res<-deg(data,pheno$sample,control="normal")
vollimma(res,nb=500,fc=0.5,p=0.05,size=4,alpha=1)
sig<-filtresig(res)
process<-reducedf(sig,data,n=287)

## perform PCA without ellipses
pcatrans(process,pheno,group="tissue",pal="Set1",alpha=0.7,names=F)

bestheat(process,pheno,font=10,rownames=F,scale="row")

pheno$id<-NULL
pheno$gsm<-NULL


write.table(res,file="limma.tsv",row.names=T,sep="\t")

write.table(sig,file="limmasig278.tsv",row.names=T,sep="\t")

## export results with metabolic annotations
head(sig)

sig$hs.gene<-row.names(sig)

all%>%right_join(sig,by="hs.gene")->final
final%>%arrange(desc(logFC))->final

library(writexl)

write_xlsx(final,"metabolism_Heptatoblastoma.xlsx")


