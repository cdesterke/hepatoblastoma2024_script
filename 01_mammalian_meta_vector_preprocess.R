list.files()

## load mammalian metabolism database
meta<-read.csv("metabolism.csv",h=T)

## library(devtools)
## devtools::install_github("cdesterke/geneconverter")

library(geneconverter)

## load biomart data from geneconverter r-package
load("bm110.rda")

library(dplyr)

bm110%>%dplyr::rename(mm.gene="gene.mm")->bm110

bm110%>%right_join(meta,by="mm.gene")->all
save(all,file="humanannotation.rda")
head(all)

dim(all)
all%>%distinct(hs.gene)->genes

vector<-as.vector(genes$hs.gene)

vector<-as.vector(na.omit(vector))

save(vector,file="meta_vector.rda")



## RNAseq preprocess

list.files()

pheno<-read.csv("pheno.csv",h=T)

data<-read.table("counts.txt",sep="\t",h=T)

bm110%>%dplyr::rename(Geneid="hs.id")->bm110

bm110%>%inner_join(data,by="Geneid")->df

df$Geneid<-NULL
df$mm.id<-NULL
df$mm.gene<-NULL
df%>%dplyr::rename(gene="hs.gene")->data

library(transpipe)
ok<-filtermatrix(data)


library(edgeR)

d0 <- DGEList(ok)

d0<-calcNormFactors(d0,method="TMM")


y<-voom(d0, plot=T)
data<-y$E

data<-as.data.frame(data)

save(data,file="voom.rda")

write.csv(data,file="TNMVOOM.csv",row.names=T)


colramp = colorRampPalette(c(3,"darkblue",5))(24)
plot(density(data[,1]),col=colramp[1],lwd=3,ylim=c(0,.2))
	for(i in 1:24){lines(density(data[,i]),lwd=1,col=colramp[i])}


edata = data[rowMeans(data) > 0, ]
plot(density(edata[,1]),col=colramp[1],lwd=3,ylim=c(0,.25))
	for(i in 1:24){lines(density(edata[,i]),lwd=1,col=colramp[i])}

save(edata,file="TNMVOOMfiltrated12386.rda")
write.csv(edata,file="TNMVOOMfiltrated12386.csv",row.names=T)



library(preprocessCore)


norm = normalize.quantiles(as.matrix(edata))

plot(density(norm[,1]),col=colramp[1],lwd=3,ylim=c(0,.25))
	for(i in 1:24){lines(density(norm[,i]),lwd=1,col=colramp[i])}



row.names(norm)<-row.names(edata)
colnames(norm)<-colnames(edata)

save(norm,file="quantile.rda")



pheno%>%filter(tissue!="Hepatoblastoma cell line")->pheno

row.names(pheno)<-pheno$id

small<-norm[,row.names(pheno)]


data<-small[row.names(small)%in%vector,]

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


head(sig)

sig$hs.gene<-row.names(sig)

all%>%right_join(sig,by="hs.gene")->final
final%>%arrange(desc(logFC))->final

library(writexl)

write_xlsx(final,"metabolism_Heptatoblastoma.xlsx")


