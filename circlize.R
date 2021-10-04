library(nichenetr)
library(Seurat)
library(ggplot2)
library(tidyr)
library(dplyr)
library(reshape2)
library(circlize)
library(RColorBrewer)
library(ggplot2)
active_ligand_target_links=active_ligand_target_links[,intersect(colnames(active_ligand_target_links),rownames(all[["RNA"]]))]
scaled=t(as.matrix(all[["RNA"]][rownames(active_ligand_target_links),]))
scaled[]=apply(scaled,2,function(x){x/max(x)})
hm_ligand=c()
for (i in unique(all@meta.data$celltype)){
int=scaled[rownames(all@meta.data)[which(all@meta.data$celltype==i)],]
int=apply(int,2,mean)
hm_ligand=rbind(hm_ligand,int)
}
rownames(hm_ligand)=unique(all@meta.data$celltype)
ligand_group=apply(hm_ligand,2,function(x){which(x==max(x))})
ligand_group=rownames(hm_ligand)[ligand_group]
names(ligand_group)=colnames(hm_ligand)

which(active_ligand_target_links[ligands_all,]!=0)]
dim(lr_network_top_matrix)
dim(active_ligand_target_links)
active_ligand_target_links=active_ligand_target_links[,-ncol(active_ligand_target_links)]
head(active_ligand_target_links)
rownames(active_ligand_target_links)

ligand_activity=data.frame(ligand_activities)[1:20,]
gr=data.frame(weighted_networks$gr)
gr=gr[which(gr[,1] %in% bg & gr[,2] %in% bg),]

circlize_data=c()
for(ligands_all in ligand_activity[,1]){
ligand_pearson=ligand_activity[which(ligand_activity[,1]==ligands_all),4]
targets_all = colnames(active_ligand_target_links)[
which(active_ligand_target_links[ligands_all,]!=0)]

receptor=lr_network_top_matrix[which(lr_network_top_matrix[,ligands_all]!=0),ligands_all,drop=F]
receptor=receptor[order(receptor[,1],decreasing=TRUE),,drop=F]

int_ligand=c()
for (r in rownames(receptor)){
receptor_activity=receptor[r,1]
tier1=gr[which(gr[,1] %in% r),2]
tier2=gr[which(gr[,1] %in% tier1 & gr[,3]>0.05),2]

target_link=intersect(targets_all,unique(c(tier1,tier2)))
if(length(target_link)>0){
int=data.frame(target=target_link,ligand=ligands_all,receptor=r,
receptor_activity=receptor_activity,ligand_pearson=ligand_pearson,
weight=t(active_ligand_target_links[ligands_all,target_link])[,1])
int_ligand=rbind(int_ligand,int)
}
}
circlize_data=rbind(circlize_data,int_ligand)
}
table(circlize_data[,3])
circlize_data=circlize_data[order(circlize_data[,4],decreasing=TRUE),]
target_grou=circlize_data[!duplicated(circlize_data[,1]),c(1,3)]
target_group=target_grou[,2]
names(target_group)=target_grou[,1]

target_group[which(target_group=="BMPR2")]="BMPR1A"
target_group[which(target_group=="TNFRSF1B")]="TNFRSF1A"

removed_rec=names(table(target_group))[which(table(target_group)<3)]
target_group=target_group[-which(target_group %in% removed_rec)]
x=rep("other receptor",66)
names(x)=setdiff(colnames(active_ligand_target_links),names(target_group))
target_group=c(target_group,x)
ligand_group=ligand_group[order(ligand_group)]
group=c(ligand_group,target_group)
group1=c(ligand_group,target_group)
group=factor(group,levels=c("Mono","Macrophage","Tcell","Pulp","other cell",
"TGFBR1","FGFR3","TGFBR2","RARG","TNFRSF1A","LTBR","BMPR1A"))

allc=c(sender_celltypes,"Osteoblast")
allc=allc[order(allc)]
ligand_col_map=data.frame(cell=allc,color=brewer.pal(12,"Set3"))
rownames(ligand_col_map)=ligand_col_map[,1]

active_ligand_target_links=t(active_ligand_target_links)[,names(ligand_group)]

group[which(group %in% c("B","Endothelium","Pericyte","Osteoblast"))]="other cell"

circos.clear()

group=group[which(group!="other receptor")]
d_circos=as.matrix(t(active_ligand_target_links))
d_circos=d_circos[which(rownames(d_circos) %in% names(group)),
which(colnames(d_circos) %in% names(group))]

cdm=chordDiagram(d_circos,group=group,
annotationTrack = c("text"),grid.col="grey",row.col=ligand_col_map[ligand_group,2])

#
circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index, 
        facing = "clockwise", niceFacing = TRUE,adj=c(-0.2,0.5),cex=0.8)
}, bg.border = NA)
#
head(cdm)
cdm[which(cdm[,2]=="H2AFX"),]
allgene=unique(c(cdm[,1],cdm[,2]))

scaled=t(as.matrix(all[["RNA"]][allgene,]))

hm_allgene=c()
for (i in c("Mono","Macrophage","Tcell","Pulp","B","Endothelium","Osteoblast","Pericyte")){
int=scaled[rownames(all@meta.data)[which(all@meta.data$seurat_clusters==i)],]
int=apply(int,2,mean)
hm_allgene=rbind(hm_allgene,int)
}
rownames(hm_allgene)=c("Mono","Macrophage","Tcell","Pulp","B","Endothelium","Osteoblast","Pericyte")
hm_allgene[]=apply(hm_allgene,2,function(x){x/max(x)})
celllabel=data.frame(gene=rownames(hm_allgene),symbol=c("T","Mono","Mac",
"Pulp","B","End","Osteoblast","Pericyte"))
rownames(celllabel)=celllabel[,1]

brewer.pal(11,"RdBu")
col_fun = colorRamp2(c(0,1), c("white","red"))
col_fun_weight=colorRamp2(c(0.1,0.7), c("#DADAEB","#54278F"))
circos.clear()
pdf("circosOB.pdf",height=8,width=8)
circos.par("cell.padding"=c(0,0,0,0),start.degree=120)

chordDiagram(d_circos,group=group,transparency =0,
annotationTrack = c(),grid.col="white",row.col=ligand_col_map[ligand_group,2],
preAllocateTracks = list(list(track.index = 1,ylim=c(0,12),track.height=cm_h(2.4),bg.border ="white")))

for(i in 1:8){
cell=rownames(hm_allgene)[i]
#label_pos=i-0.5

#circos.text(sector.index="FASLG",track.index =1,-0.01,
#            label_pos, labels = celllabel[cell,2],facing = "inside", 
#            adj = c(1.2,0),cex=0.6)
y1=i-1
y2=i
for(gene in names(ligand_group)){
if(gene %in% cdm[,1]){
		x=max(cdm[which(cdm[,1]==gene),"x1"])
	}else{
		x=max(cdm[which(cdm[,2]==gene),"x2"])
	}

circos.rect(x,y1,0,y2, 
col = col_fun(hm_allgene[cell,gene]), 
border = "black",
sector.index = gene, track.index = 1)
}
}
circos.yaxis(side="left", at=seq(0.5,7.5,length.out=8),labels=celllabel[c(2,3,1,4:8),2],
sector.index="ITGAM", track.index=1,tick=FALSE,labels.cex=0.6,lwd=0,col="white")

#circos.text(sector.index="ALPL",track.index =1,-0.01,
#            0, labels = "receptor\nweight",facing = "inside", 
#           adj = c(1,0),cex=0.5)

y1=0
y2=1
for(gene in setdiff(colnames(d_circos),names(ligand_group))){
if(gene %in% cdm[,1]){
		x=max(cdm[which(cdm[,1]==gene),"x1"])
	}else{
		x=max(cdm[which(cdm[,2]==gene),"x2"])
	}
rec=target_group[gene]
rec=max(lr_network_top_matrix[rec,1:20])
col=col_fun_weight(rec)
circos.rect(x,y1,0,y2, col = col, border = col,
sector.index = gene, track.index = 1)
}


circos.track(track.index = 1, panel.fun = function(x, y) {
    gene=CELL_META$sector.index
    if(gene %in% names(ligand_group)){
		yaxis=9.1
		adj=0.3
            col="black"
            cex=0.6}else{
		yaxis=1.1
		adj=0.05
		col=ifelse(cort[gene]<0,"#00BFC4","#F8766D")
            cex=0.5}
    circos.text(CELL_META$xcenter,yaxis , gene, col=col,
        facing = "clockwise", niceFacing = TRUE,adj=c(adj,0.5),cex=cex)
}, bg.border = NA)

for(i in levels(group)[6:12]){
genes=names(group)[which(group==i)]
if(i %in% unique(target_group)){genes=setdiff(genes,names(ligand_group))}
lg=max(nchar(genes))

highlight.sector(genes, track.index = 1,
cex = 0.6, text.col = "black", niceFacing = TRUE,col = "#FF000000",
text = i, text.vjust=0.1,text.hjust=0.5)
}

dev.off()
-1.2*lg



lgd1=Legend(col_fun = col_fun, title = "mean\nexpression", at = c(0,0.5,1),
direction="horizontal",title_position ="topcenter")
lgd2=Legend(col_fun = col_fun_weight, title = "receptor\nweight", at = c(0.1,0.4,0.7),
direction="vertical",title_position ="leftcenter-rot")


pd = packLegend(lgd1,lgd2, direction = "vertical", 
row_gap = unit(5, "mm"))
p_legend=grid::grid.grabExpr(draw(pd))

scaled=t(as.matrix(TC[["RNA"]][c("FASLG","IFNG","TGFB1"),]))
hm_TC=c()
for (i in unique(TC@meta.data$seurat_clusters)[1:7]){
int=scaled[rownames(TC@meta.data)[which(TC@meta.data$seurat_clusters==i)],]
int=apply(int,2,mean)
hm_TC=rbind(hm_TC,int)
}
rownames(hm_TC)=unique(TC@meta.data$seurat_clusters)[1:7]
hm_TC[]=apply(hm_TC,2,function(x){x/max(x)})
HMtc=Heatmap(t(hm_TC),col=col_fun,cluster_rows = FALSE,cluster_columns=FALSE,
row_names_side="left",row_names_gp=gpar(size=5),column_names_gp=gpar(size=5),
column_names_rot=-45)

pdf("obcircos_hm.pdf",height=3,width=7.5)
draw(HMtc,annotation_legend_list=pd)
dev.off()

scaled=t(as.matrix(pulp[["RNA"]][c("BMP2","SFRP2","BMP7"),]))
hm_pulp=c()
for (i in unique(pulp@meta.data$seurat_clusters)){
int=scaled[rownames(pulp@meta.data)[which(pulp@meta.data$seurat_clusters==i)],]
int=apply(int,2,mean)
hm_pulp=rbind(hm_pulp,int)
}
rownames(hm_pulp)=unique(pulp@meta.data$seurat_clusters)
hm_pulp[]=apply(hm_pulp,2,function(x){x/max(x)})
HMpulp=Heatmap(t(hm_pulp),col=col_fun,cluster_rows = FALSE,cluster_columns=FALSE,
row_names_side="left",row_names_gp=gpar(size=5),column_names_gp=gpar(size=5),
column_names_rot=-45)

pdf("obcircos_hm.pdf",height=3,width=7.5)
draw(HMtc,annotation_legend_list=pd)
dev.off()
########################################################
library(data.table)
grn=fread("All.grn_output.tsv",data.table=FALSE)
aucell=fread("All.aucell_output.tsv",data.table=FALSE)
head(grn)

lr_network_top_matrix[1:3,1:3]

lr_niche=gather(lr_network_top_matrix[names(table(target_group)),],ligand,weight,1:20)
lr_niche=lr_niche[which(lr_niche[,3]!=0),]
cdb_index=c(1,2,grep("|iOsteoblast",colnames(sigmean),fixed=TRUE),
grep("|dOsteoblast",colnames(sigmean),fixed=TRUE))

re="TNFRSF1A"
sigmean[which(sigmean[,1]==re),1:3]
lr_niche[which(lr_niche[,1]==re),]


for(i in 1:nrow(sigmean)){
if(sigmean[i,4]){
a=sigmean[i,1]
b=sigmean[i,2]
sigmean[i,1]=b
sigmean[i,2]=a
}
}
sigmean=sigmean[,-3]
sigmean=sigmean[,-3]
sigmean[1:3,1:5]
#################################################
head(circlize_data)
head(gr)
intersect(grn[,1],gr[which(gr[,1] %in% names(table(target_group))),2])
table(grn[,1])
grn=grn[which(grn[,1] %in% bg1 & grn[,2] %in% bg),]
nichenet_scenic=c()
for(i in names(table(target_group))){
allds=gr[which(gr[,1]==i),2]
if(length(intersect(grn[,1],allds))>0){
int=circlize_data[which(circlize_data[,3]==i),]
targets=int[,1]
tf=intersect(grn[,1],allds)
inttf=grn[which(grn[,1] %in% tf & grn[,2] %in% targets &grn[,3]>1),]
intrt=gr[which(gr[,1]==i & gr [,2] %in% tf),]
rownames(intrt)=intrt[,2]
inttf$receptor=i
inttf$rtweight=intrt[inttf[,1],3]
nichenet_scenic=rbind(nichenet_scenic,inttf)
}
}
head(nichenet_scenic)

###########################################################
table(nichenet_scenic[,2])
library(Biobase)
rownames(pData(HSMM))
regulon=sub("(+)","",aucell[,1],fixed=TRUE)
rownames(aucell)=regulon
HSMM1=HSMM[regulon,]
exprs(HSMM1)=as.matrix(aucell[,rownames(pData(HSMM))])

diff_tf_time <- differentialGeneTest(HSMM1,
fullModelFormulaStr = "~sm.ns(Pseudotime)")
sigTF=rownames(diff_tf_time)[which(diff_tf_time$qval<0.01)]
############################################################

sigTF=intersect(names(table(nichenet_scenic[,1])),ob_tgene)
nichenet_scenic=nichenet_scenic[which(nichenet_scenic[,1] %in% sigTF),]
write.csv(nichenet_scenic,"dscenic.csv")
dscenic=fread("dscenic.csv",data.table=F)
dscenic=spread(dscenic,TF,importance)
head(dscenic)
rownames(dscenic)=dscenic[,1]
dscenic=dscenic[,-1]
dscenic[]=apply(dscenic,2,function(x){
x[which(is.na(x))]=1
return(x)
})
dscenic=dscenic[which(apply(dscenic,1,max)>10),]
dscenic=log10(dscenic)

dscenic=dscenic[,order(match(colnames(dscenic),rownames(dpthm)))]
dscenic=dscenic[order(match(rownames(dscenic),rownames(dpthm))),]

col_scenic= colorRamp2(c(0,0.4,0.8,1.2,1.6,2),
brewer.pal(11,"Spectral")[c(11,9,7,5,3,1)])

hmscenic=Heatmap(dscenic,col=col_scenic,heatmap_legend_param = list(at = c(0,1,2),
title="TF regulation",direction = "vertical",title_position="leftcenter-rot",title_gp=gpar(fontsize=8)),
show_column_names=F,cluster_columns=FALSE,cluster_rows=FALSE,
row_names_gp=gpar(fontsize=5,col=ifelse(cort[rownames(dscenic)]<0,"#00BFC4","#F8766D")))

lgd1=Legend(col_fun = col_fun, title = "mean expression", at = c(0,0.5,1),
direction="vertical",title_position ="leftcenter-rot")
lgd3=Legend(col_fun = col_fun_weight, title = "link weight", at = c(0.1,0.4,0.7),
direction="vertical",title_position ="leftcenter-rot")
,
merge_legend = TRUE,annotation_legend_list=list(lgd1,lgd3))
hmS=grid::grid.grabExpr(draw(hmscenic,heatmap_legend_side = "left"))
library(patchwork)
pdf("p4d.pdf",height=5,width=3)
plot_spacer()/hmS+plot_layout(height=c(2,3))
dev.off()

lr_net=gather(lr_network_top_matrix[names(table(target_group)),],
ligand,weight,1:20)
head(lr_net)
lr_net=lr_net[which(lr_net[,3]!=0),]
lr_net[which(lr_net[,2]=="BMP7"),2]="BMP7(L)"
lr_net=lr_net[,c(2,1,3)]
colnames(lr_net)=c("from","to","weight")
rt_net=gr[which(gr[,1] %in% names(table(target_group)) &
gr[,2] %in% colnames(dscenic)),]

annotation=rbind(data.frame(node=unique(rt_net[,1]),type="receptor"),
data.frame(node=unique(rt_net[,2]),type="TF"))

annotation=rbind(annotation,
data.frame(node=unique(lr_net[,1]),type="ligand"))

net=rbind(lr_net,rt_net)
write.csv(annotation,"annotation.csv")
write.csv(net,"net.csv")



