memory.limit(999999)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(viridis)
library(ggsci)
all=readRDS("all.rds")
head(all@meta.data)
all@meta.data[,10]=as.numeric(all@meta.data[,10])
all@meta.data[,10]=all@meta.data[,10]-1
ll=names(table(all@meta.data$seurat_clusters))
ll[11]="SOX9+"
all@meta.data[which(all@meta.data$seurat_clusters=="Pulp"),"seurat_clusters"]="SOX9+"
all@meta.data$seurat_clusters=factor(all@meta.data$seurat_clusters,
levels=ll)

DimPlot(all, reduction = "umap",label=T)+theme(legend.position="none")
p_UMAP=DimPlot(all, reduction = "umap",group.by="seurat_clusters",label=TRUE)+
p_UMAP+theme(
legend.position="none",axis.text=element_blank(),strip.background = element_blank(),
axis.title=element_blank(),axis.ticks=element_blank())+scale_color_brewer(palette="Set3")
axis.ticks=element_blank(),axis.text=element_blank(),
tiff("allUMAP.tiff",res=300,units="in",height=4,width=4)
p_UMAP
dev.off()

pdf("allUMAP.pdf",height=4,width=4)
p_UMAP
dev.off()

FeaturePlot(all,reduction = "umap", features = c("CGF1","ENAM","DMP1","AMBN"),
cols = viridis(11)[c(1,6,11)])+theme(legend.position="NULL")

all@meta.data[which(all@meta.data[,10] %in% c(1,2,4,5)),10]="Tcell"
all@meta.data[which(all@meta.data[,10] %in% c(0,3)),10]="Neutrophil"
all@meta.data[which(all@meta.data[,10] %in% c("Bcell")),10]="B"
all@meta.data[which(all@meta.data[,10] %in% c(17)),10]="pDC"
all@meta.data[which(all@meta.data[,10] %in% c(6)),10]="DC"
all@meta.data[which(all@meta.data[,10] %in% c(12)),10]="Macrophage"
all@meta.data[which(all@meta.data[,10] %in% c("Monocyte")),10]="Mono"
all@meta.data[which(all@meta.data[,10] %in% c(14)),10]="Osteoclast"
all@meta.data[which(all@meta.data[,10] %in% c(8)),10]="Endothelium"
all@meta.data[which(all@meta.data[,10] %in% c("Fibroblast")),10]="Osteoblast"
all@meta.data[which(all@meta.data[,10] %in% c(9)),10]="Pulp"
all@meta.data[which(all@meta.data[,10] %in% c(16)),10]="Pericyte"


draw=t(as.matrix(all[["RNA"]][c("CD14","CD1C","LILRA4","S100A9","CCL3","CD3E",
"CD79A","FCN1","ACP5","VWF","SPARC","SOX9","RGS5"),]))
y=Embeddings(all[["umap"]])
draw=data.frame(draw,x=y[rownames(draw),1],y=y[rownames(draw),2])
library(tidyr)
dmm=gather(draw,marker,value,CD14:RGS5, factor_key=TRUE)
mm=ggplot(dmm,aes(x,y,color=value))+
geom_point(size=0.1)+facet_wrap(~marker, ncol = 7)+
scale_color_gradientn(breaks=c(0,0.2,0.4,0.6,0.8,1),
colors=brewer.pal(11,"Spectral")[c(11,9,7,5,3,1)])+theme_minimal()+theme(
axis.text=element_blank(),strip.text=element_text(
size=10),strip.background=element_blank(),legend.position = "none",
axis.title=element_blank(),axis.ticks=element_blank())
tiff("marker.tiff",res=300,units="in",height=4,width=14)
mm
dev.off()


get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}
dat=Embeddings(all[["umap"]])[rownames(all@meta.data)[which(all@meta.data[,1]=="PatientB")],]
colnames(dat)=c("x","y")
dat=data.frame(dat)
dat$density <- get_density(dat$x, dat$y, n = 100)
dat$density=max(dat$density)-dat$density
dab=ggplot(dat) + geom_point(aes(x, y, color = density),size=3.5) + 
scale_color_viridis(option="A") +theme_minimal()

head(all@meta.data)
dbar=all@meta.data[,c(1,10)]
colnames(dbar)=c("id","cluster")
dbar$id=as.character(dbar$id)
dbar$id=ifelse(dbar$id=="PatientA","adult","child")
dbar$id=factor(dbar$id,levels=c("child","adult"))

abar=ggplot(dbar,aes(id,fill=cluster))+geom_bar(position="fill")
abar=abar+scale_fill_brewer(palette="Set3")+theme_minimal()
abar=abar+theme(axis.ticks=element_blank(),axis.title=element_blank(),
axis.text=element_text(size=10,color="black"))
abar=abar+theme(legend.position="none",axis.text.y=element_blank(),
axis.title.y=element_text(size=10))+ylab("proportion")


#pdf("test.pdf",height=3,width=8)
f1_mid=abar+mm+plot_layout(width=c(1,7))
#dev.off()
f1_up=plot_spacer()+p_UMAP+plot_layout(width=c(5,3))



#########################################################################
TC=rownames(all@meta.data)[which(all@meta.data[,10] %in% c("Tcell"))]
TC=SubsetData(all,cells=TC)

TC<- FindVariableFeatures(TC, selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
TC<- RunPCA(TC, npcs = 30, verbose = FALSE)
ElbowPlot(TC)
TC<- RunUMAP(TC, reduction = "pca", dims = 1:13)
TC<- FindNeighbors(TC, dims = 1:13)
TC<- FindClusters(TC, resolution = 0.6)
TC=SubsetData(TC,cells=rownames(TC@meta.data)[which(TC@meta.data[,10]!=8)])
DimPlot(TC, reduction = "umap",label=T)+theme(legend.position="none")
FeaturePlot(TC,reduction = "umap", features = c("CD103","CD69","CCR7","CD45RA"),
cols = viridis(11)[c(1,6,11)])+theme(legend.position="NULL")

cluster1.markers=FindMarkers(TC, ident.1 =6,only.pos = TRUE,min.pct = 0.25)
head(cluster1.markers, n = 15)


#0 CD8 cytotoxic GZMH+NKG7
#1 CD8 memory IL7R+LTB
#2 NK NKG7+GNLY
#3 CD4 Th17 CCR6+CXCR6
#4 CD4 memory CCR7+SELL
#5 CD8 activation CRTAM+CCR5
#6 NK naive SELL
#7 MKI67 proliferation
TC@meta.data[,10]=as.numeric(TC@meta.data[,10])-1
TC@meta.data[which(TC@meta.data[,10]==0),10]="cytotoxic CD8T"
TC@meta.data[which(TC@meta.data[,10]==1),10]="memory CD8T"
TC@meta.data[which(TC@meta.data[,10]==2),10]="NK"
TC@meta.data[which(TC@meta.data[,10]==3),10]="Th17"
TC@meta.data[which(TC@meta.data[,10]=="naive CD4T"),10]="memory CD4T"
TC@meta.data[which(TC@meta.data[,10]==5),10]="activated CD8T"
TC@meta.data[which(TC@meta.data[,10]==6),10]="naive NK"
TC@meta.data[which(TC@meta.data[,10]==7),10]="proliferation"

T_UMAP=DimPlot(TC,pt.size=1, reduction = "umap",group.by="seurat_clusters",label=T)+
theme(legend.position="none",axis.title=element_blank(),
axis.ticks=element_blank())+scale_color_brewer(palette="Set3")
tiff("TUMAP.tiff",res=300,units="in",height=5,width=5)
T_UMAP
dev.off()

draw=t(as.matrix(TC[["RNA"]][c("CD4","CD8A","NKG7","GZMH","IL7R","LTB",
"GNLY","CCR6","CCR7","SELL","CRTAM","CXCR6","MKI67"),]))
y=Embeddings(TC[["umap"]])
draw=data.frame(draw,x=y[rownames(draw),1],y=y[rownames(draw),2])
library(tidyr)
dmm=gather(draw,marker,value,CD4:MKI67, factor_key=TRUE)
Tmm=ggplot(dmm,aes(x,y,color=value))+
geom_point(size=0.5)+facet_wrap(~marker, ncol = 7)+
scale_color_gradientn(
colors=brewer.pal(11,"Spectral")[c(11,9,7,5,3,1)])+theme_minimal()+theme(
axis.text=element_blank(),strip.text=element_text(
size=6),strip.background=element_blank(),
axis.title=element_blank(),axis.ticks=element_blank())
tiff("Tmarker.tiff",res=300,units="in",height=2,width=7)
Tmm
dev.off()
legend.position="none",
hvg<- FindVariableFeatures(all, selection.method = "vst", 
        nfeatures = 5000, verbose = FALSE)@assays$RNA@var.features

dbar=TC@meta.data[,c(1,10)]
colnames(dbar)=c("id","cluster")
dbar$id=as.character(dbar$id)
dbar$id=ifelse(dbar$id=="PatientA","adult","child")
dbar$id=factor(dbar$id,levels=c("child","adult"))

tbar=ggplot(dbar,aes(id,fill=cluster))+geom_bar(position="fill")
tbar=tbar+scale_fill_brewer(palette="Set3")+theme_minimal()
tbar=tbar+theme(axis.ticks=element_blank(),axis.title=element_blank(),
axis.text=element_text(size=10,color="black"),legend.position="none")
tbar

ints=rownames(TC@meta.data)[which(TC@meta.data[,10] %in% c("memory CD4T","memory CD8T"))]
dvln=data.frame(value=t(as.matrix(TC[["RNA"]][c("CD69"),ints])),
patient=ifelse(TC@meta.data[ints,1]=="PatientA","adult","child"),
cluster=TC@meta.data[ints,10])
dvln$patient=factor(dvln$patient,levels=c("child","adult"))
dvln$cluster=factor(dvln$cluster,levels=unique(TC@meta.data[,10]))


xx=c("memory CD8T:p=1.1e-9","memory CD4T:p=4.1e-5")
names(xx)=c("memory CD8T","memory CD4T")
pvln=ggplot(dvln,aes(patient,CD69,color=cluster,fill=cluster))+geom_violin(alpha=0.5)
pvln=pvln+geom_jitter()+theme_minimal()+facet_wrap(.~cluster,ncol=1,
labeller = labeller(cluster=xx))+
theme(axis.ticks=element_blank(),axis.title.x=element_blank(),
axis.text=element_text(size=10,color="black"),legend.position="none")
pvln=pvln+scale_color_manual(values=c("#FB8072","#BEBADA"))
pvln=pvln+scale_fill_manual(values=c("#FB8072","#BEBADA"))
pvln
t.test(CD69~patient,data=dvln[which(dvln[,3]=="memory CD8T"),])

T_UMAP+tbar+pvln+plot_layout(width=c(1.1,0.3,0.6))

metaT=TC@meta.data[,10,drop=F]

metaT=metaT[which(metaT[,1]!="proliferation"),,drop=F]
metaO=all@meta.data[which(all@meta.data[,10]!="Tcell"),10,drop=F]
metaT=rbind(metaT,metaO)
d=(all[["RNA"]]@counts)[hvg,rownames(metaT)]
write.table(d,"cdbcount.txt",quote=F,sep="\t")
write.table(metaT,"cdbmetat.txt",quote=F,sep="\t")

###############################################################
head(all@meta.data)

pulp=rownames(all@meta.data)[which(all@meta.data[,10] %in% c("Pulp"))]
pulp=SubsetData(all,cells=pulp)

pulp<- FindVariableFeatures(pulp, selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
pulp<- RunPCA(pulp, npcs = 30, verbose = FALSE)
ElbowPlot(pulp)
pulp<- RunUMAP(pulp, reduction = "pca", dims = 1:6)
pulp<- FindNeighbors(pulp, dims = 1:6)
pulp<- FindClusters(pulp, resolution = 0.1)

pulpUMAP=DimPlot(pulp, reduction = "umap",group.by="seurat_clusters",
label=TRUE)+theme(
legend.position="none",axis.text=element_blank(),strip.background = element_blank(),
axis.title=element_blank(),axis.ticks=element_blank())


FeaturePlot(pulp,reduction = "umap", features = c("CD13","CD24","CLU","AMBN"),
cols = viridis(11)[c(1,6,11)])+theme(legend.position="NULL")

cluster1.markers=FindMarkers(pulp, ident.1 =0,only.pos = TRUE,min.pct = 0.25)
head(cluster1.markers, n = 25)

draw=t(as.matrix(pulp[["RNA"]][c("ID2","CD24","CLU","AMBN"),]))
y=Embeddings(pulp[["umap"]])
draw=data.frame(draw,x=y[rownames(draw),1],y=y[rownames(draw),2])
draw[,1:4]=apply(draw[,1:4],2,function(x){x/max(x)})
library(tidyr)
dmm=gather(draw,marker,value,ID2:AMBN, factor_key=TRUE)
Pmm=ggplot(dmm,aes(x,y,color=value))+
geom_point(size=1)+facet_wrap(~marker, ncol = 2)+
scale_color_gradientn(breaks=c(0,0.2,0.4,0.6,0.8,1),
colors=brewer.pal(11,"Spectral")[c(11,9,7,5,3,1)])+theme_classic()+theme(
legend.position="none",axis.text=element_blank(),strip.background = element_blank(),
axis.title=element_blank(),axis.ticks=element_blank())
tiff("Pmarker.tiff",res=300,units="in",height=2.7,width=3)
Pmm
dev.off()

f1e=pulpUMAP/Pmm



pulp@meta.data[,10]=as.numeric(pulp@meta.data[,10])-1
pulp@meta.data[which(pulp@meta.data[,10]==0),10]="APSC"
pulp@meta.data[which(pulp@meta.data[,10]==1),10]="ameloblast"
metaP=pulp@meta.data[,10,drop=F]

###############################################################
head(all@meta.data)

ob=rownames(all@meta.data)[which(all@meta.data[,10] %in% c("Osteoblast"))]
ob=SubsetData(all,cells=ob)

ob<- FindVariableFeatures(ob, selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
ob<- RunPCA(ob, npcs = 30, verbose = FALSE)
ElbowPlot(ob)
ob<- RunUMAP(ob, reduction = "pca", dims = 1:11)
ob<- FindNeighbors(ob, dims = 1:11)
ob<- FindClusters(ob, resolution = 0.1)
ob=SubsetData(ob,cells=rownames(ob@meta.data)[which(ob@meta.data[,10]!=2)])
obUMAP=DimPlot(ob, reduction = "umap",group.by="seurat_clusters",label=TRUE)+theme(
legend.position="none",axis.text=element_blank(),strip.background = element_blank(),
axis.title=element_blank(),axis.ticks=element_blank())

FeaturePlot(pulp,reduction = "umap", features = c("DLK1","VIM","GJA1","SPARC"),
cols = viridis(11)[c(1,6,11)])+theme(legend.position="NULL")

cluster1.markers=FindMarkers(pulp, ident.1 =0,only.pos = TRUE,min.pct = 0.25)
head(cluster1.markers, n = 20)

draw=t(as.matrix(ob[["RNA"]][c("SPON1","VIM","GJA1","SPARC"),]))
y=Embeddings(ob[["umap"]])
draw=data.frame(draw,x=y[rownames(draw),1],y=y[rownames(draw),2])
library(tidyr)
dmm=gather(draw,marker,value,SPON1:SPARC, factor_key=TRUE)
Omm=ggplot(dmm,aes(x,y,color=value))+
geom_point(size=0.5)+facet_wrap(~marker, ncol = 2)+
scale_color_gradientn(breaks=c(0,0.2,0.4,0.6,0.8,1),
colors=brewer.pal(11,"Spectral")[c(11,9,7,5,3,1)])+theme_classic()+theme(
legend.position="none",axis.text=element_blank(),strip.background = element_blank(),
axis.title=element_blank(),axis.ticks=element_blank())
tiff("Omarker.tiff",res=300,units="in",height=3,width=3)
Omm
dev.off()

f1f=obUMAP/Omm
design="
AB##
CD##"
f1_down=(f1e|f1f|plot_spacer())+plot_layout(width=c(1,1,2))
f1=(f1_up/plot_spacer()/f1_mid/plot_spacer()/f1_down)+plot_layout(height=c(3,0.3,2.4,0.3,4))
pdf("Figure 1.pdf",height=10,width=8)
f1
dev.off()

tiff("f1e.tiff",res=300,units="in",height=5,width=5)
f1e|f1f
dev.off()


ob@meta.data[,10]=as.numeric(ob@meta.data[,10])-1
ob@meta.data[which(ob@meta.data[,10]=="DOsteoblast"),10]="dOsteoblast"
ob@meta.data[which(ob@meta.data[,10]=="dOsteoblast"),10]="iOsteoblast"
metaOB=ob@meta.data[,10,drop=F]

metaO=all@meta.data[-which(all@meta.data[,10] %in% c("Osteoblast","Pulp")),10,drop=F]
metaO=rbind(metaO,metaP)
metaO=rbind(metaO,metaOB)
meta=rbind(metaT,metaO[which(metaO[,1]!="Tcell"),1,drop=F])
write.table(meta,"cdbmetat.txt",quote=F,sep="\t")
d=(all[["RNA"]]@counts)[hvg,rownames(meta)]
write.table(d,"cdbcount.txt",quote=F,sep="\t")

#############################
N=rownames(all@meta.data)[which(all@meta.data[,10] %in% c("Neutrophil"))]
N=SubsetData(all,cells=N)

N<- FindVariableFeatures(N, selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
N<- RunPCA(N, npcs = 30, verbose = FALSE)
ElbowPlot(N)
N<- RunUMAP(N, reduction = "pca", dims = 1:15)
N<- FindNeighbors(N, dims = 1:15)
N<- FindClusters(N, resolution = 0.3)
N=SubsetData(N,cells=rownames(N@meta.data)[which(N@meta.data[,10]!=5)])
DimPlot(N, reduction = "umap",group.by="seurat_clusters",label=TRUE)+
theme(legend.position="none")
FeaturePlot(N,reduction = "umap", features = c("S100P","MKI67","PRRG4","P2RY13"),
cols = viridis(11)[c(1,6,11)])+theme(legend.position="NULL")

cluster1.markers=FindMarkers(N, ident.1 =0,only.pos = TRUE,min.pct = 0.25)
head(cluster1.markers, n = 15)


#0 unclassifed
#1 S100P+
#2 P2RY13+
#3 antiviral MX1
#4 PRRG4+
#5 
#6 innate PGLYRP1
#7 inhibitory SPLI+PI

N@meta.data[,10]=as.numeric(N@meta.data[,10])-1
N@meta.data[which(N@meta.data[,12]==0),10]="unclassifed"
N@meta.data[which(N@meta.data[,12]==1),10]="S100P+"
N@meta.data[which(N@meta.data[,12]==2),10]="P2RY13+"
N@meta.data[which(N@meta.data[,12]==3),10]="antiviral"
N@meta.data[which(N@meta.data[,12]==4),10]="PRRG4+"

N@meta.data[which(N@meta.data[,12]==6),10]="innate"
N@meta.data[which(N@meta.data[,12]==7),10]="inhibitory"


metaN=N@meta.data[,10,drop=F]


metaO=metaO[which(metaO[,1]!="Neutrophil"),1,drop=F]
metaN=rbind(metaN,metaO)
d=(all[["RNA"]]@counts)[hvg,rownames(metaN)]
write.table(d,"Ncount.txt",quote=F,sep="\t")
write.table(metaN,"metaN.txt",quote=F)
write.csv(metaN,"metaN.csv")

draw=t(as.matrix(N[["RNA"]][c("S100P","P2RY13","MX1","PRRG4","PGLYRP1","SLPI"),]))
y=Embeddings(N[["umap"]])
draw=data.frame(draw,x=y[rownames(draw),1],y=y[rownames(draw),2])
library(tidyr)
dmm=gather(draw,marker,value,S100P:SLPI, factor_key=TRUE)
Nmm=ggplot(dmm,aes(x,y,color=value))+
geom_point(size=0.5)+facet_wrap(~marker, ncol = 2)+
scale_color_gradientn(
colors=brewer.pal(11,"Spectral")[c(11,9,7,5,3,1)])+theme_minimal()+theme(
axis.text=element_blank(),strip.text=element_text(
size=8),strip.background=element_blank(),legend.key.size = unit(0.3, 'cm'),
axis.title=element_blank(),axis.ticks=element_blank())
tiff("Tmarker.tiff",res=300,units="in",height=2,width=7)
Nmm
dev.off()

dbar=N@meta.data[,c(1,10)]
colnames(dbar)=c("id","cluster")
dbar$id=as.character(dbar$id)
dbar$id=ifelse(dbar$id=="PatientA","adult","child")
dbar$id=factor(dbar$id,levels=c("child","adult"))

nbar=ggplot(dbar,aes(id,fill=cluster))+geom_bar(position="fill")
nbar=nbar+scale_fill_brewer(palette="Set3")+theme_minimal()
nbar=nbar+theme(axis.ticks=element_blank(),axis.title=element_blank(),
axis.text=element_text(size=10,color="black"),legend.position="none")
nbar

N_UMAP=DimPlot(N,pt.size=1, reduction = "umap",group.by="seurat_clusters",label=TRUE)+
theme(legend.position="none",axis.title=element_blank(),
axis.ticks=element_blank())+scale_color_brewer(palette="Set3")

library(patchwork)
p_f3_one=N_UMAP+nbar+Nmm+plot_layout(width=c(0.9,0.3,0.8))
p_f3_two=plot_spacer()+p_endN_bub
p_f3_three=plot_spacer()+p_BD_bub

pdf("Figure3.pdf",height=10,width=8)
p_f3_one/p_f3_two/p_f3_three
dev.off()

######################################################
SC=SubsetData(all,cells=APSC)
Idents(SC)="orig.ident"

head(SC@meta.data)
table(SC@meta.data$Phase)

ggplot(SC@meta.data,aes(S.Score,G2M.Score,color=Phase,shape=orig.ident))+geom_point()+
theme_minimal()

SC<- FindVariableFeatures(SC, selection.method = "vst", 
        nfeatures = 2000, verbose = FALSE)
hvg=SC@assays$RNA@var.features

library(monocle)
fd=data.frame(gene_short_name=rownames(SC[["RNA"]]))
rownames(fd)=fd[,1]
pd <- new("AnnotatedDataFrame", data = SC@meta.data)
fd <- new("AnnotatedDataFrame", data = fd)
HSMM <- newCellDataSet(as.matrix(SC[["RNA"]][rownames(fd),rownames(pd)]),
    phenoData = pd, featureData = fd,expressionFamily=uninormal())
diff_test_res <- differentialGeneTest(HSMM[hvg,],
              fullModelFormulaStr = "~Phase")
diff_test_res=diff_test_res[order(diff_test_res[,3]),]
ordering_genes <- row.names (subset(diff_test_res, pval < 0.05))
HSMM <- setOrderingFilter(HSMM, ordering_genes)
HSMM <- reduceDimension(HSMM, max_components = 2,norm_method="none",
    method = 'DDRTree')
HSMM <- orderCells(HSMM)

pData(HSMM)[which(pData(HSMM)$S.Score>0.2),"S.Score"]=0.2
pData(HSMM)[which(pData(HSMM)$G2M.Score>0.2),"G2M.Score"]=0.2
pData(HSMM)[which(pData(HSMM)$S.Score<(-0.1)),"S.Score"]=-0.1
pData(HSMM)[which(pData(HSMM)$G2M.Score<(-0.1)),"G2M.Score"]=-0.1
pData(HSMM)[which(pData(HSMM)$Phase=="S"),"Phase"]="G1/S"
pData(HSMM)[which(pData(HSMM)$Phase=="G1"),"Phase"]="G0"
pData(HSMM)[which(pData(HSMM)$Phase=="G2M"),"Phase"]="G2/M"

p_sscore=plot_cell_trajectory(HSMM, color_by = "S.Score")+annotate(geom="text",
x=-0.1,y=-1.5,label="G1/S Score")+scale_color_gradientn(breaks=c(-0.1,-0.05,0,0.05,0.1,0.15),
colors=brewer.pal(11,"Spectral")[c(11,9,7,5,3,1)])+theme_minimal()+ylim(c(-2,1))+
theme(legend.position="none",axis.title=element_blank(),axis.text.x=element_blank())


p_gscore=plot_cell_trajectory(HSMM, color_by = "G2M.Score")+annotate(geom="text",
x=-0.1,y=-1.5,label="G2/M Score")+scale_color_gradientn(breaks=c(-0.1,-0.05,0,0.05,0.1,0.15),
colors=brewer.pal(11,"Spectral")[c(11,9,7,5,3,1)])+theme_minimal()+ylim(c(-2,1))+
theme(legend.position="none",axis.title=element_blank(),axis.text.x=element_blank())

p_monocle=plot_cell_trajectory(HSMM, color_by = "Phase")+theme_minimal()+
theme(legend.position=c(0.4,0.3),legend.spacing.y=unit(0,"mm"),
legend.key.size = unit(1, "mm"))

lb=rownames(pData(HSMM))[which(pData(HSMM)$State!=1)]
rb=rownames(pData(HSMM))[which(pData(HSMM)$State!=2)]

diff_lb_time <- differentialGeneTest(HSMM[hvg,lb],
fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_lb_time=diff_lb_time[order(diff_lb_time[,3]),]
rownames(diff_lb_time)[which(diff_lb_time$pval<0.01)]

diff_rb_time <- differentialGeneTest(HSMM[hvg,rb],
fullModelFormulaStr = "~sm.ns(Pseudotime)")

diff_rb_time=diff_rb_time[order(diff_rb_time[,3]),]
rownames(diff_rb_time)[which(diff_rb_time$pval<0.01)]

head(diff_lb_time)

head(pData(HSMM))
pData(HSMM)$bin=c()
pData(HSMM)=pData(HSMM)[order(pData(HSMM)$Pseudotime),]
lbin=ceiling((1:94)/10)
lbin[91:94]=9
rbin=ceiling((1:33)/10)
rbin[31:33]=3
pData(HSMM)[which(pData(HSMM)$State==2),"bin"]=lbin+10
pData(HSMM)[which(pData(HSMM)$State==3),"bin"]=10
pData(HSMM)[which(pData(HSMM)$State==1),"bin"]=rbin

lbgene=rownames(diff_lb_time)[which(diff_lb_time$pval<0.01)]
rbgene=rownames(diff_rb_time)[which(diff_rb_time$pval<0.01)]

d_pthm=exprs(HSMM)["DKK3",lb]

d_time=data.frame(exp=pthm,time=pData(HSMM)[names(pthm),12],phase=
pData(HSMM)[names(pthm),7])
p_apsc_c=ggplot(d_time,aes(time,exp,color=phase,group=1))+geom_point()+geom_smooth(
data=d_time,se=F,color="black")+theme_minimal()+
xlab("pseudotime")+ylab("DKK3")+theme(legend.position="none")

[which(d_time[,1]!=0),]
method="lm",
table(pData(HSMM)$bin)
f=function(x){
id=rownames(pData(HSMM))[which(pData(HSMM)$bin==x)]
int=apply(d_pthm[,id],1,mean)
return(int)
}
dpthm=sapply(1:19,f)
colnames(dpthm)=1:19
dpthm=dpthm[,10:19]
ord=order(unlist(apply(dpthm,1,function(x){which(x==max(x))[1]})))
dpthm=dpthm[ord,]
dpthm[]=apply(dpthm,1,function(x){x/max(x)})

cort=apply(dpthm,1,function(x){cor(x,10:19,method="spearman")})
cort=cort[order(cort)]
arg=data.frame(diff_lb_time[names(cort),],cort=cort)
write.csv(arg,"Table SARG.csv")


gene.d =bitr(names(cort)[which(cort<0)], fromType = "SYMBOL", 
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db)
gl=as.character(gene.d[,2])
dGO=enrichGO(gene         =gl,
                OrgDb         = org.Hs.eg.db,
                
                ont           = "BP",
                pAdjustMethod = "BH",
              
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)
rGO=data.frame(simplify(dGO))
rGO[1:30,c(2,6)]

source("https://gist.githubusercontent.com/jokergoo/bfb115200df256eeacb7af302d4e508e/raw/14f315c7418f3458d932ad749850fd515dec413b/word_cloud_grob.R")
getR=function(x){
x1=unlist(strsplit(as.character(x[1]),"/",fixed=TRUE))
x1=as.numeric(x1[1])/as.numeric(x1[2])
x2=unlist(strsplit(as.character(x[2]),"/",fixed=TRUE))
x2=as.numeric(x2[1])/as.numeric(x2[2])
r=x1/x2
return(r)
}

library(stringr)
words=gGO[c("GO:0009612","GO:0007162","GO:0016049","GO:0010463",
"GO:0043687","GO:0048511"),c(2)]
pval=-log10(gGO[c("GO:0009612","GO:0007162","GO:0016049","GO:0010463",
"GO:0043687","GO:0048511"),c(6)])

r=apply(gGO[c("GO:0009612","GO:0007162","GO:0016049","GO:0010463",
"GO:0043687","GO:0048511"),3:4],1,getR)
words=str_wrap(words, width = 35)

col_fun= colorRamp2(c(3,6),brewer.pal(11,"RdBu")[c(11,1)])

gb_g = word_cloud_grob(words, fontsize = 2.5*sqrt(r), max_width = unit(45, "mm"),
col=col_fun(pval))
grid.newpage()
grid.draw(gb_g)
grid.rect(width = 1.5*grobWidth(gb_g), height = 2*grobHeight(gb_g), gp = gpar(fill = NA))

#
write.csv(rGO,"rGO.csv")
write.csv(gGO,"gGO.csv")

words=rGO[c("GO:0061448","GO:0006027","GO:0030111","GO:0042533",
"GO:0030856","GO:0001667"),c(2)]
pval=-log10(rGO[c("GO:0061448","GO:0006027","GO:0030111","GO:0042533",
"GO:0030856","GO:0001667"),c(6)])

r=apply(rGO[c("GO:0061448","GO:0006027","GO:0030111","GO:0042533",
"GO:0030856","GO:0001667"),3:4],1,getR)
words=str_wrap(words, width = 35)

col_fun= colorRamp2(c(2,7),brewer.pal(11,"RdBu")[c(11,1)])

gb_r = word_cloud_grob(words, fontsize = 2.5*sqrt(r), max_width = unit(45, "mm"),
col=col_fun(pval))
grid.newpage()
grid.draw(gb_r)
grid.rect(width = 1.5*grobWidth(gb_r), height = 2*grobHeight(gb_r), gp = gpar(fill = NA))

breaks=c(rep("descend",82),rep("ascend",62))
breaks=factor(breaks,levels=c("descend","ascend"))

panel_fun = function(index, nm) {
     grid.rect(gp = gpar(fill = "#EEEEEE", col = NA))
     if(nm=="descend"){
    grid.draw(gb_r)
    }else{grid.draw(gb_g)
     }
 }
ra=rowAnnotation(word_cloud = anno_link(align_to =breaks, which = "row", 
        panel_fun = panel_fun, size =1.2*grobHeight(gb_r), 
        width = 1.2*grobWidth(gb_r), 
        link_gp = gpar(fill = "#EEEEEE", col = NA)
    ))

hm=Heatmap(dpthm[names(cort)[order(cort)],],col=col_hm,heatmap_legend_param = list(at = c(0,0.5,1),
title="scaled expression",direction = "horizontal",title_position="leftcenter"),show_row_names=F,
cluster_rows = FALSE,row_split=breaks,cluster_columns=FALSE,show_column_names=F,
right_annotation=ra)

pthm = grid::grid.grabExpr(draw(hm,heatmap_legend_side = "bottom"))


#######################
library(patchwork)
pf6a=p_monocle/p_sscore/p_gscore+plot_layout(height=c(2,1,1))
pf6a=pf6a|plot_spacer()
pdf("figure 6a.pdf",width=8,height=5)
pf6a+plot_layout(width=c(3,5))
dev.off()

pdf("figure 6b.pdf",width=5,height=5)
draw(hm,heatmap_legend_side = "bottom")
dev.off()

#############################################
library(nichenetr)
rm(lr_network)
rm(ligand_target_matrix)
rm(weighted_networks)
rm(weighted_networks_lr)
rm(all)

Idents(all)="seurat_clusters"
all@meta.data[APSC,"seurat_clusters"]="APSC"
all@meta.data[ame,"seurat_clusters"]="amelo"
receiver = "APSC"

expressed_genes_receiver = get_expressed_genes(receiver, all, pct = 0.25)

background_expressed_genes = expressed_genes_receiver %>% .[. %in% rownames(ligand_target_matrix)]
bg=background_expressed_genes

sender_celltypes=levels(Idents(all))[-8]

list_expressed_genes_sender = sender_celltypes %>% unique() %>% lapply(get_expressed_genes, all, 0.10) # lapply to get the expressed genes of every sender cell type separately here
expressed_genes_sender = list_expressed_genes_sender %>% unlist() %>% unique()

geneset_oi = names(cort)
geneset_oi = geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]

ligands = lr_network %>% pull(from) %>% unique()
receptors = lr_network %>% pull(to) %>% unique()

expressed_ligands = intersect(ligands,expressed_genes_sender)
expressed_receptors = intersect(receptors,expressed_genes_receiver)

potential_ligands = lr_network %>% filter(from %in% expressed_ligands & to %in% expressed_receptors) %>% pull(from) %>% unique()

ligand_activities = predict_ligand_activities(geneset = geneset_oi, background_expressed_genes = background_expressed_genes, 
ligand_target_matrix = ligand_target_matrix, potential_ligands = potential_ligands)
ligand_activities = ligand_activities %>% arrange(-pearson) %>% mutate(rank = rank(desc(pearson)))

best_upstream_ligands = ligand_activities %>% top_n(20, pearson) %>% arrange(-pearson) %>% pull(test_ligand) %>% unique()

active_ligand_target_links_df = best_upstream_ligands %>% lapply(get_weighted_ligand_target_links,geneset = geneset_oi, ligand_target_matrix = ligand_target_matrix, n = 400) %>% bind_rows() %>% drop_na()

active_ligand_target_links = prepare_ligand_target_visualization(ligand_target_df = active_ligand_target_links_df, ligand_target_matrix = ligand_target_matrix, cutoff = 0.33)

order_ligands = intersect(best_upstream_ligands, colnames(active_ligand_target_links)) %>% rev() %>% make.names()
order_targets = active_ligand_target_links_df$target %>% unique() %>% intersect(rownames(active_ligand_target_links)) %>% make.names()
rownames(active_ligand_target_links) = rownames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23
colnames(active_ligand_target_links) = colnames(active_ligand_target_links) %>% make.names() # make.names() for heatmap visualization of genes like H2-T23

vis_ligand_target = active_ligand_target_links[order_targets,order_ligands] %>% t()
p_ligand_target_network = vis_ligand_target %>% make_heatmap_ggplot("Prioritized ligands","Predicted target genes", color = "purple",legend_position = "top", x_axis_position = "top",legend_title = "Regulatory potential")  + theme(axis.text.x = element_text(face = "italic")) + scale_fill_gradient2(low = "whitesmoke",  high = "purple", breaks = c(0,0.0045,0.0090))
p_ligand_target_network

dimnames(active_ligand_target_links)

lr_network_top = lr_network %>% filter(from %in% best_upstream_ligands & to %in% expressed_receptors) %>% distinct(from,to)
best_upstream_receptors = lr_network_top %>% pull(to) %>% unique()

lr_network_top_df_large = weighted_networks_lr %>% filter(from %in% best_upstream_ligands & to %in% best_upstream_receptors)

lr_network_top_df = lr_network_top_df_large %>% spread("from","weight",fill = 0)
lr_network_top_matrix = lr_network_top_df %>% dplyr::select(-to) %>% as.matrix() %>% magrittr::set_rownames(lr_network_top_df$to)

dist_receptors = dist(lr_network_top_matrix, method = "binary")
hclust_receptors = hclust(dist_receptors, method = "ward.D2")
order_receptors = hclust_receptors$labels[hclust_receptors$order]
    
dist_ligands = dist(lr_network_top_matrix %>% t(), method = "binary")
hclust_ligands = hclust(dist_ligands, method = "ward.D2")
order_ligands_receptor = hclust_ligands$labels[hclust_ligands$order]

order_receptors = order_receptors %>% intersect(rownames(lr_network_top_matrix))
order_ligands_receptor = order_ligands_receptor %>% intersect(colnames(lr_network_top_matrix))

lr_network_top_matrix = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(lr_network_top_matrix) = order_receptors %>% make.names()
colnames(lr_network_top_matrix) = order_ligands_receptor %>% make.names()

lr_network_top_matrix=data.frame(lr_network_top_matrix,receptor=rownames(lr_network_top_matrix))
vis_ligand_receptor_network = lr_network_top_matrix[order_receptors, order_ligands_receptor]
rownames(vis_ligand_receptor_network) = order_receptors %>% make.names()
colnames(vis_ligand_receptor_network) = order_ligands_receptor %>% make.names()

p_ligand_receptor_network = vis_ligand_receptor_network %>% make_heatmap_ggplot("Ligands","Receptors", color = "mediumvioletred", x_axis_position = "top",legend_title = "Prior interaction potential")
p_ligand_receptor_network=p_ligand_receptor_network+theme(legend.position="none")
p_ligand_target_network=p_ligand_target_network+theme(legend.position="none")



library(patchwork)

p_apsc=((p_apsc_a/p_apsc_b/p_apsc_c)|p_ligand_receptor_network)/p_ligand_target_network
save(p_apsc,file="p_apsc.RData")
pdf("Figure SAPSC.pdf",height=10,width=8)
p_apsc
dev.off()





