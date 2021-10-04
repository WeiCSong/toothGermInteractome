library(data.table)
sigmean=fread("significant_means.txt",data.table=F)
sigmean=sigmean[,c(5:9,11:412)]
colnames(sigmean)

f=function(x){
length(which(!is.na(x)))
}

MAX=function(x){
max(x[which(!is.na(x))])
}

type=function(x){
if(length(which(!is.na(x)))==0){
return(NA)}else{
x[which(is.na(x))]=0
return(sigmean[which(x==max(x)),3])
}
}


N=apply(sigmean[,8:407],2,f)
M=apply(sigmean[,8:407],2,MAX)
M[which(M==-Inf)]=0
T=apply(sigmean[,8:407],2,type)
netT=data.frame(N,M)
head(netT)
Ttype=c("NK","naive NK","cytotoxic CD8T","activated CD8T","memory CD8T",
"memory CD4T","Th17")
write.csv(netT,"netT.csv")
netT=fread("netT.csv",data.table=F)
netT$filter=ifelse(netT$C1==netT$C2|!(netT$C1 %in% Ttype),0,1)

ame=c("NK|ameloblast",
"Th17|ameloblast",
"activated CD8T|ameloblast",
"cytotoxic CD8T|ameloblast",
"memory CD4T|ameloblast",
"memory CD8T|ameloblast",
"naive NK|ameloblast")

ost=c("NK|Osteoclast",
"Th17|Osteoclast",
"activated CD8T|Osteoclast",
"cytotoxic CD8T|Osteoclast",
"memory CD4T|Osteoclast",
"memory CD8T|Osteoclast",
"naive NK|Osteoclast")

amerow=apply(sigmean[,ame],1,function(x){length(which(!is.na(x)))})

ameB=sigmean[which(amerow!=0),ame]
rownames(ameB)=paste(sigmean[which(amerow!=0),1],sigmean[which(amerow!=0),2],sep="|")
ameB[]=apply(ameB,2,function(x){
x[which(is.na(x))]=0
return(x)
})
ameB=log10(ameB+1)

pT=fread("pvalues.txt",data.table=F)
pT=pT[which(pT[,5]!=""&pT[,6]!=""),]
pT=pT[-100,]
rownames(pT)=paste(pT[,5],pT[,6],sep="|")
pT=pT[rownames(ameB),colnames(ameB)]

ameBU$LRpair=rownames(ameBU)
library(tidyr)

ameBU=gather(ameBU,cellPair, strength,1:7)
pT$LRpair=rownames(pT)
pT=gather(pT,cellPair, strength,1:7)
colnames(pT)[3]="pval"
d_ameT_bub=data.frame(ameBU,pval=pT[,3])
library(ggplot2)
d_ameT_bub[,2]=gsub("|ameloblast","",d_ameT_bub[,2])
d_ameT_bub[,2]=gsub("|","",d_ameT_bub[,2],fixed=TRUE)
d_ameT_bub[,4]=1-d_ameT_bub[,4]
p=ggplot(d_ameT_bub,aes(cellPair,LRpair,size=pval,color=strength))+geom_point()
p=p+theme_minimal()+scale_color_gradientn(breaks=c(0,0.3,0.6,0.9,1.2),colors=
brewer.pal(9,"Spectral")[c(9,7,5,3,1)])+scale_size(breaks=c(0,0.5,1),labels=c(1,0.5,0))
p=p+xlab("cell interacted with ameloblast")+theme(axis.title.y=element_blank(),
axis.text.x=element_text(size=8,color="black",angle=90,hjust=1),axis.text.y=element_text(
size=8,color="black"))
p_ameT_bub=p
#p_ostT_bub=p_ostT_bub+xlab("osteoclast partner")
#p_ameT_bub=p_ameT_bub+theme(legend.position = "bottom")


ostrow=apply(sigmean[,ost],1,function(x){length(which(!is.na(x)))})

ostB=sigmean[which(ostrow!=0),ost]
rownames(ostB)=paste(sigmean[which(ostrow!=0),1],sigmean[which(ostrow!=0),2],sep="|")
ostB[]=apply(ostB,2,function(x){
x[which(is.na(x))]=0
return(x)
})
ostB=log10(ostB+1)

pT=fread("pvalues.txt",data.table=F)
pT=pT[which(pT[,5]!=""&pT[,6]!=""),]
pT=pT[-100,]
rownames(pT)=paste(pT[,5],pT[,6],sep="|")
pT=pT[rownames(ostB),colnames(ostB)]

ostBU=ostB
ostBU$LRpair=rownames(ostBU)
library(tidyr)

ostBU=gather(ostBU,cellPair, strength,1:7)
pT$LRpair=rownames(pT)
pT=gather(pT,cellPair, strength,1:7)
colnames(pT)[3]="pval"
d_ostT_bub=data.frame(ostBU,pval=pT[,3])
library(ggplot2)
d_ostT_bub=d_ostT_bub[which(d_ostT_bub[,1]!="THBS1|"),]
d_ostT_bub[,2]=gsub("|Osteoclast","",d_ostT_bub[,2],fixed=TRUE)
d_ostT_bub[,4]=1-d_ostT_bub[,4]
p=ggplot(d_ostT_bub,aes(cellPair,LRpair,size=pval,color=strength))+geom_point()
p=p+theme_minimal()+scale_color_gradientn(breaks=c(0,0.3,0.6,0.9,1.2),colors=
brewer.pal(9,"Spectral")[c(9,7,5,3,1)])+scale_size(breaks=c(0,0.5,1),labels=c(1,0.5,0))
p=p+xlab("cell interacted with osteoclast")+theme(axis.title.y=element_blank(),
axis.text.x=element_text(size=8,color="black",angle=90,hjust=1),axis.text.y=element_text(
size=8,color="black"),legend.position="none")
p_ostT_bub=p
p_ostT_bub=p_ostT_bub+scale_y_manual(labels=c("unclassifed","S100P+","))
library(patchwork)
p_Tbub=p_ameT_bub  +p_ostT_bub + plot_layout(guides = "collect") & 
theme(legend.position = "bottom")

p_f2_up=T_UMAP+tbar+pvln+plot_layout(width=c(3.5,1.5,2))
p_f2_down=plot_spacer()+p_Tbub+plot_layout(width=c(1,2))
p_f2=p_f2_up/Tmm/p_f2_down+plot_layout(height=c(3,2,4))
tiff("F2.tiff",res=300,units="in",height=10,width=8)
p_f2
dev.off()

###########################################################
library(data.table)
sigmean=fread("significant_means.txt",data.table=F)
sigmean=sigmean[,c(5:9,13:412)]
colnames(sigmean)

f=function(x){
length(which(!is.na(x)))
}

MAX=function(x){
max(x[which(!is.na(x))])
}

type=function(x){
if(length(which(!is.na(x)))==0){
return(NA)}else{
x[which(is.na(x))]=0
return(sigmean[which(x==max(x)),3])
}
}


N=apply(sigmean[,8:405],2,f)
M=apply(sigmean[,8:405],2,MAX)
M[which(M==-Inf)]=0
T=apply(sigmean[,8:405],2,type)
netN=data.frame(N,M)
head(netN)
Ntype=c("unclassifed","S100P+","P2RY13+","antiviral","PRRG4+","innate",
"inhibitory")
int=t(data.frame(strsplit(rownames(netN),"|",fixed=TRUE)))
netN=data.frame(netN,int)
colnames(netN)[3:4]=c("C1","C2")
netN$filter=ifelse(netN$C1==netN$C2|!(netN$C1 %in% Ntype),0,1)
write.csv(netN,"netN.csv")


ame=c("unclassifed|Osteoclast",
"S100P+|Osteoclast",
"P2RY13+|Osteoclast",
"antiviral|Osteoclast",
"PRRG4+|Osteoclast",
"innate|Osteoclast",
"inhibitory|Osteoclast")



amerow=apply(sigmean[,ame],1,function(x){length(which(!is.na(x)))})

ameB=sigmean[which(amerow!=0),ame]
rownames(ameB)=paste(sigmean[which(amerow!=0),1],sigmean[which(amerow!=0),2],sep="|")
ameB[]=apply(ameB,2,function(x){
x[which(is.na(x))]=0
return(x)
})
ameB=log10(ameB+1)

pN=fread("pvalues.txt",data.table=F)
pN=pN[which(pN[,5]!=""&pN[,6]!=""),]
pN=pN[-100,]
rownames(pN)=paste(pN[,5],pN[,6],sep="|")
pN=pN[rownames(ameB),colnames(ameB)]

ameB$LRpair=rownames(ameB)
library(tidyr)

ameB=gather(ameB,cellPair, strength,1:7)
pN$LRpair=rownames(pN)
pN=gather(pN,cellPair, strength,1:7)
colnames(pN)[3]="pval"
d_ameN_bub=data.frame(ameB,pval=pN[,3])
library(ggplot2)
d_ameN_bub[,2]=gsub("|ameloblast","",d_ameN_bub[,2],fixed=TRUE)
d_ameN_bub[,4]=1-d_ameN_bub[,4]
p=ggplot(d_ameN_bub,aes(LRpair,cellPair,size=pval,color=strength))+geom_point()
p=p+theme_minimal()+scale_color_gradientn(breaks=c(0,0.3,0.6,0.9,1.2),colors=
brewer.pal(9,"Spectral")[c(9,7,5,3,1)])+scale_size(breaks=c(0,0.5,1),labels=c(1,0.5,0))
p=p+ylab("cell interacted with\nOsteoclast")+theme(axis.title.x=element_blank(),
axis.text.x=element_text(size=8,color="black",hjust=1,angle=45),axis.text.y=element_text(
size=8,color="black"),legend.key.size = unit(0.3, 'cm'))
p_OSCN_bub=p
p_OSCN_bub=p_OSCN_bub+scale_y_discrete(labels=c("unclassified","S100P+","PRRG4+",
"P2RY13+","innate","inhibitory","antiviral")[7:1])+labs(color="")+theme(
legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,0,-10,-10))
#p_ostT_bub=p_ostT_bub+xlab("osteoclast partner")
#p_ameT_bub=p_ameT_bub+theme(legend.position = "bottom")

###########################################################################

Itype=c("B","DC","pDC","Mono","Macrophage")
Dtype=c("ameloblast","APSC","dOsteoblast","Endothelium","iOsteoblast",
"Osteoclast","Pericyte")
int=t(data.frame(strsplit(rownames(netN),"|",fixed=TRUE)))
netN=data.frame(netN,int)
colnames(netN)[3:4]=c("C1","C2")
netN$filter1=ifelse(netN$C1 %in% Itype & netN$C2 %in% Dtype,1,0)
write.csv(netN,"netI.csv")


ID=paste("B",Dtype,sep="|")



amerow=apply(sigmean[,ID],1,function(x){length(which(!is.na(x)))})

ameB=sigmean[which(amerow!=0)[-1],ID]
rownames(ameB)=paste(sigmean[which(amerow!=0)[-1],1],sigmean[which(amerow!=0)[-1],2],sep="|")
ameB[]=apply(ameB,2,function(x){
x[which(is.na(x))]=0
return(x)
})
ameB=log10(ameB+1)

pN=fread("pvalues.txt",data.table=F)
pN=pN[which(pN[,5]!=""&pN[,6]!=""),]
pN=pN[-100,]
rownames(pN)=paste(pN[,5],pN[,6],sep="|")
pN=pN[rownames(ameB),colnames(ameB)]

ameB$LRpair=rownames(ameB)
library(tidyr)

ameB=gather(ameB,cellPair, strength,1:7)
pN$LRpair=rownames(pN)
pN=gather(pN,cellPair, strength,1:7)
colnames(pN)[3]="pval"
d_BD_bub=data.frame(ameB,pval=pN[,3])
library(ggplot2)
d_BD_bub[,2]=gsub("B|","",d_BD_bub[,2],fixed=TRUE)
d_BD_bub[,4]=1-d_BD_bub[,4]
p=ggplot(d_BD_bub,aes(LRpair,cellPair,size=pval,color=strength))+geom_point()
p=p+theme_minimal()+scale_color_gradientn(breaks=c(0,0.3,0.6,0.9,1.2,1.5),colors=
brewer.pal(11,"Spectral")[c(11,9,7,5,3,1)])+scale_size(breaks=c(0,0.5,1),labels=c(1,0.5,0))
p=p+ylab("cell interacted with\nB cell")+theme(axis.title.x=element_blank(),
axis.text.x=element_text(size=8,color="black",hjust=1,angle=45),axis.text.y=element_text(
size=8,color="black"),legend.margin=margin(0,0,0,0),
legend.box.margin=margin(-10,0,-10,-10),legend.key.size = unit(0.3, 'cm'))+labs(color="")
p_BD_bub=p
save(p_BD_bub,file="p_BD_bub.RData")


###########################################################################

library(data.table)
sigmean=fread("significant_means.txt",data.table=F)
sigmean=sigmean[,c(5:9,13:412)]
colnames(sigmean)

f=function(x){
length(which(!is.na(x)))
}

MAX=function(x){
max(x[which(!is.na(x))])
}

type=function(x){
if(length(which(!is.na(x)))==0){
return(NA)}else{
x[which(is.na(x))]=0
return(sigmean[which(x==max(x)),3])
}
}

sigmean[which(sigmean[,1]=="FFAR2"),c(1,2,grep("))]

N=apply(sigmean[,8:405],2,f)
M=apply(sigmean[,8:405],2,MAX)
M[which(M==-Inf)]=0
T=apply(sigmean[,8:405],2,type)
netN=data.frame(N,M)
head(netN)
Ntype=c("unclassifed","S100P+","P2RY13+","antiviral","PRRG4+","innate",
"inhibitory")
int=t(data.frame(strsplit(rownames(netN),"|",fixed=TRUE)))
netN=data.frame(netN,int)
colnames(netN)[3:4]=c("C1","C2")
netN$filter=ifelse(netN$C1==netN$C2|!(netN$C1 %in% Ntype),0,1)
write.csv(netN,"netN.csv")


ame=c("unclassifed|Endothelium",
"S100P+|Endothelium",
"P2RY13+|Endothelium",
"antiviral|Endothelium",
"PRRG4+|Endothelium",
"innate|Endothelium",
"inhibitory|Endothelium")



amerow=apply(sigmean[,ame],1,function(x){length(which(!is.na(x)))})

ameB=sigmean[which(amerow!=0),ame]
rownames(ameB)=paste(sigmean[which(amerow!=0),1],sigmean[which(amerow!=0),2],sep="|")
ameB[]=apply(ameB,2,function(x){
x[which(is.na(x))]=0
return(x)
})
ameB=log10(ameB+1)

pN=fread("pvalues.txt",data.table=F)
pN=pN[which(pN[,5]!=""&pN[,6]!=""),]
pN=pN[-100,]
rownames(pN)=paste(pN[,5],pN[,6],sep="|")
pN=pN[rownames(ameB),colnames(ameB)]

ameB$LRpair=rownames(ameB)
library(tidyr)

ameB=gather(ameB,cellPair, strength,1:7)
pN$LRpair=rownames(pN)
pN=gather(pN,cellPair, strength,1:7)
colnames(pN)[3]="pval"
d_endN_bub=data.frame(ameB,pval=pN[,3])
library(ggplot2)
d_endN_bub[,2]=gsub("|Endothelium","",d_endN_bub[,2],fixed=TRUE)
d_endN_bub[,4]=1-d_endN_bub[,4]
p=ggplot(d_endN_bub,aes(LRpair,cellPair,size=pval,color=strength))+geom_point()
p=p+theme_minimal()+scale_color_gradientn(breaks=c(0,0.3,0.6,0.9,1.2,1.5),colors=
brewer.pal(11,"Spectral")[c(11,9,7,5,3,1)])+scale_size(breaks=c(0,0.5,1),labels=c(1,0.5,0))
p=p+ylab("cell interacted with\nEndothelium")+theme(axis.title.x=element_blank(),
axis.text.x=element_text(size=8,color="black",hjust=1,angle=45),axis.text.y=element_text(
size=8,color="black"),legend.key.size = unit(0.3, 'cm'))
p_endN_bub=p_endN_bub+labs(color="")+theme(
legend.margin=margin(0,0,0,0),legend.box.margin=margin(-10,0,-10,-10))
#p_ostT_bub=p_ostT_bub+xlab("osteoclast partner")

##################################################
rec=c("MET","PLXND1","LTBR","FLT1","KDR")

gr=gr[which(gr[,1] %in% rec),]
head(gr)
table(gr[,1])

TC=rownames(all@meta.data)[which(all@meta.data[,10] %in% c("Endothelium"))]
all=SubsetData(all,cells=TC)
exp=all[["RNA"]][bg,]
exp=as(exp,"matrix")
dim(exp)
exp[1:3,1:3]
COR=cor(t(exp))
COR=COR[,rec]
head(COR)

receptor="MET"
dgene=gr[which(gr[,1]==receptor),2]
dgene[(which(COR[dgene,receptor]>0.1))]

gene.d =bitr(dgene[(which(COR[dgene,receptor]>0))], fromType = "SYMBOL", 
        toType = "ENTREZID",
        OrgDb = org.Hs.eg.db)
gl=as.character(gene.d[,2])
dGO=enrichGO(gene         =gl,
                OrgDb         = org.Hs.eg.db,
                
                ont           = "BP",
                pAdjustMethod = "BH",
              
                pvalueCutoff  = 0.05,
                qvalueCutoff  = 0.05)
metGO=data.frame(simplify(dGO))
ltbrGO
kdrGO
flt1GO
metGO[1:30,c(2,6)]

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
words=ltbrGO[c("GO:0120192","GO:0045333","GO:0071496"),c(2)]
pval=-log10(ltbrGO[c("GO:0030111","GO:0060485","GO:0001667","GO:0030509",
"GO:0050673","GO:0016049","GO:0090287"),c(6)])

r=apply(ltbrGO[c("GO:0030111","GO:0060485","GO:0001667","GO:0030509",
"GO:0050673","GO:0016049","GO:0090287"),3:4],1,getR)
words=str_wrap(words, width = 40)


col_fun= colorRamp2(c(3,6),brewer.pal(11,"RdBu")[c(11,1)])

gb_i = word_cloud_grob(words, fontsize = 2.2*r, max_width = unit(50, "mm"),
col=col_fun(pval))
grid.newpage()
grid.draw(gb_i)
grid.rect(width = 1.5*grobWidth(gb), height = 2*grobHeight(gb), gp = gpar(fill = NA))


