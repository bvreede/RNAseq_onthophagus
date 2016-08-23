##redoing names for large table

data.large<-read.table(file.choose(),header=T)

colnames(data.large)<-c("gene","smDsx_Male_Brain_S97","smDsx_Male_Brain_S98",
"smDsx_Male_Brain_S99",
"smDsx_Male_Brain_S100",
"smDsx_Male_Brain_S101",
"smDsx_Male_Brain_S102",

"smCntl_Male_Brain_S103",
"smCntl_Male_Brain_S104",
"smCntl_Male_Brain_S105",
"smCntl_Male_Brain_S106",
"smCntl_Male_Brain_S107",
"smCntl_Male_Brain_S108",

"smDsx_Male_TH_S109",
"smDsx_Male_TH_S110",
"smDsx_Male_TH_S111",
"smDsx_Male_TH_S112",
"smDsx_Male_TH_S113",
"smDsx_Male_TH_S114",

"smCntl_Male_TH_S115",
"smCntl_Male_TH_S116",
"smCntl_Male_TH_S117",
"smCntl_Male_TH_S118",
"smCntl_Male_TH_S119",
"smCntl_Male_TH_S120",

"smDsx_Male_HH_S121",
"smDsx_Male_HH_S122",
"smDsx_Male_HH_S123",
"smDsx_Male_HH_S124",
"smDsx_Male_HH_S125",
"smDsx_Male_HH_S126",

"smCntl_Male_HH_S127",
"smCntl_Male_HH_S128",
"smCntl_Male_HH_S129",
"smCntl_Male_HH_S130",
"smCntl_Male_HH_S131",
"smCntl_Male_HH_S132",

"smDsx_Male_Gen_S133",
"smDsx_Male_Gen_S134",
"smDsx_Male_Gen_S135",
"smDsx_Male_Gen_S136",
"smDsx_Male_Gen_S137",
"smDsx_Male_Gen_S138",

"smCntl_Male_Gen_S139",
"smCntl_Male_Gen_S140",
"smCntl_Male_Gen_S141",
"smCntl_Male_Gen_S142",
"smCntl_Male_Gen_S143",
"smCntl_Male_Gen_S144")

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")

data.large<-read.table(file.choose(),header=T)

data.merge<-merge(data.large,data,by="ENSEMBL_ID")
rownames(data.merge)<-data.merge$ENSEMBL_ID
data.merge<-data.merge[,2:145]

x<-colnames(data.large)
y<-colnames(data)
both.names<-as.vector(c(x[2:97],y[2:49]))
both.treatment<-as.vector(c(rep("Cntl",48),rep("Dsx",48),rep("Dsx",6),rep("Cntl",6),rep("Dsx",6),rep("Cntl",6),rep("Dsx",6),rep("Cntl",6),rep("Dsx",6),rep("Cntl",6)))
both.size<-as.vector(c(rep("Large",96),rep("Small",48)))
both.sex<-as.vector(c(rep("Male",24),rep("Female",24),rep("Male",24),rep("Female",24),rep("Male",48)))
both.tissue<-
as.vector(
c(rep("Brain",6),rep("TH",6),rep("HH",6),rep("Gen",6),rep("Brain",6),rep("TH",6),rep("HH",6),rep("Gen",6),rep("Brain",6),rep("TH",6),rep("HH",6),rep("Gen",6),rep("Brain",6),rep("TH",6),rep("HH",6),rep("Gen",6),rep("Brain",12),rep("TH",12),rep("HH",12),rep("Gen",12))
)

metadata<-cbind(both.names,both.treatment,both.size,both.sex,both.tissue)
coldata<-metadata
colnames(coldata)<-c("Sample","Treatment","Size","Sex","Tissue")
rownames(coldata)<-as.vector(c(x[2:97],y[2:49]))
coldata<-coldata[,2:5]
coldata<-as.data.frame(coldata)

ddsFullCountTable<-DESeqDataSetFromMatrix(
countData=data.merge,
colData=coldata,
design= ~Treatment)

dds<-ddsFullCountTable
dds <- dds[ rowSums(counts(dds)) > 1, ]
dds<-DESeq(dds)

###Control
##Male
#HH
coldata.MHH.Cntl<-coldata[coldata$Sex=="Male"&coldata$Tissue=="HH"&coldata$Treatment=="Cntl",]

data.M<-data.merge[,grep("_Male",colnames(data.merge))]
data.MHH<-data.M[,grep("_HH",colnames(data.M))]
data.MHH.Cntl<-data.MHH[,grep("Cntl_",colnames(data.MHH))]

ddsMHH.Cntl<-DESeqDataSetFromMatrix(
countData=data.MHH.Cntl,
colData=coldata.MHH.Cntl,
design= ~Size)

dds<-ddsMHH.Cntl
dds<- dds[ rowSums(counts(dds)) > 1, ]
#redefine small as first variable, so positive fold change values will be large males, and negative will be small males
dds$Size<-relevel(dds$Size,"Small")
dds<-DESeq(dds)
res<-results(dds)
res<-res[complete.cases(res$padj),]

res.sig<-res[res$padj<=0.05,]
Gene_ID<-as.data.frame(rownames(res.sig))
colnames(Gene_ID)<-"Gene_ID"
res.sig<-as.data.frame(cbind(Gene_ID,res.sig))

write.table(res.sig,file="/Users/crisledon-rettig/Desktop/Indiana University/DSX Project 2/Contrasts/LCM vs SCM/HH_large_cntl_male_vs_small_cntl_male.txt",row.names=FALSE,sep="\t")

###Control
##Male
#TH
coldata.MTH.Cntl<-coldata[coldata$Sex=="Male"&coldata$Tissue=="TH"&coldata$Treatment=="Cntl",]

data.M<-data.merge[,grep("_Male",colnames(data.merge))]
data.MTH<-data.M[,grep("_TH",colnames(data.M))]
data.MTH.Cntl<-data.MTH[,grep("Cntl_",colnames(data.MTH))]

ddsMTH.Cntl<-DESeqDataSetFromMatrix(
countData=data.MTH.Cntl,
colData=coldata.MTH.Cntl,
design= ~Size)

dds<-ddsMTH.Cntl
dds<- dds[ rowSums(counts(dds)) > 1, ]
#redefine small as first variable, so positive fold change values will be large males, and negative will be small males
dds$Size<-relevel(dds$Size,"Small")
dds<-DESeq(dds)
res<-results(dds)
res<-res[complete.cases(res$padj),]

res.sig<-res[res$padj<=0.05,]
Gene_ID<-as.data.frame(rownames(res.sig))
colnames(Gene_ID)<-"Gene_ID"
res.sig<-as.data.frame(cbind(Gene_ID,res.sig))

write.table(res.sig,file="/Users/crisledon-rettig/Desktop/Indiana University/DSX Project 2/Contrasts/LCM vs SCM/TH_large_cntl_male_vs_small_cntl_male.txt",row.names=FALSE,sep="\t")

###Control
##Male
#Brain
coldata.MBRN.Cntl<-coldata[coldata$Sex=="Male"&coldata$Tissue=="Brain"&coldata$Treatment=="Cntl",]

data.M<-data.merge[,grep("_Male",colnames(data.merge))]
data.MBRN<-data.M[,grep("_Brain",colnames(data.M))]
data.MBRN.Cntl<-data.MBRN[,grep("Cntl_",colnames(data.MBRN))]

ddsMBRN.Cntl<-DESeqDataSetFromMatrix(
countData=data.MBRN.Cntl,
colData=coldata.MBRN.Cntl,
design= ~Size)

dds<-ddsMBRN.Cntl
dds<- dds[ rowSums(counts(dds)) > 1, ]
#redefine small as first variable, so positive fold change values will be large males, and negative will be small males
dds$Size<-relevel(dds$Size,"Small")
dds<-DESeq(dds)
res<-results(dds)
res<-res[complete.cases(res$padj),]

res.sig<-res[res$padj<=0.05,]
Gene_ID<-as.data.frame(rownames(res.sig))
colnames(Gene_ID)<-"Gene_ID"
res.sig<-as.data.frame(cbind(Gene_ID,res.sig))

write.table(res.sig,file="/Users/crisledon-rettig/Desktop/Indiana University/DSX Project 2/Contrasts/LCM vs SCM/BRN_large_cntl_male_vs_small_cntl_male.txt",row.names=FALSE,sep="\t")

###Control
##Male
#Gen
coldata.MGEN.Cntl<-coldata[coldata$Sex=="Male"&coldata$Tissue=="Gen"&coldata$Treatment=="Cntl",]

data.M<-data.merge[,grep("_Male",colnames(data.merge))]
data.MGEN<-data.M[,grep("_Gen",colnames(data.M))]
data.MGEN.Cntl<-data.MGEN[,grep("Cntl_",colnames(data.MGEN))]

ddsMGEN.Cntl<-DESeqDataSetFromMatrix(
countData=data.MGEN.Cntl,
colData=coldata.MGEN.Cntl,
design= ~Size)

dds<-ddsMGEN.Cntl
dds<- dds[ rowSums(counts(dds)) > 1, ]
#redefine small as first variable, so positive fold change values will be large males, and negative will be small males
dds$Size<-relevel(dds$Size,"Small")
dds<-DESeq(dds)
res<-results(dds)
res<-res[complete.cases(res$padj),]

res.sig<-res[res$padj<=0.05,]
Gene_ID<-as.data.frame(rownames(res.sig))
colnames(Gene_ID)<-"Gene_ID"
res.sig<-as.data.frame(cbind(Gene_ID,res.sig))

write.table(res.sig,file="/Users/crisledon-rettig/Desktop/Indiana University/DSX Project 2/Contrasts/LCM vs SCM/GEN_large_cntl_male_vs_small_cntl_male.txt",row.names=FALSE,sep="\t")






###Control_dsx
##Male
#HH
coldata.MHH.CD<-coldata[coldata$Sex=="Male"&coldata$Tissue=="HH"&coldata$Size=="Small",]

data.M<-data.merge[,grep("_Male",colnames(data.merge))]
data.MHH<-data.M[,grep("_HH",colnames(data.M))]
data.MHH.Small<-data.MHH[,grep("sm",colnames(data.MHH))]

ddsMHH.Sm<-DESeqDataSetFromMatrix(
countData=data.MHH.Small,
colData=coldata.MHH.CD,
design= ~Treatment)

dds<-ddsMHH.Sm
dds<- dds[ rowSums(counts(dds)) > 1, ]
dds$Size<-relevel(dds$Treatment,"Cntl")
dds<-DESeq(dds)
res<-results(dds)
res<-res[complete.cases(res$padj),]

res.sig<-res[res$padj<=0.05,]
#Gene_ID<-as.data.frame(rownames(res.sig))
#colnames(Gene_ID)<-"Gene_ID"
#res.sig<-as.data.frame(cbind(Gene_ID,res.sig))

#write.table(res.sig,file="/Users/crisledon-rettig/Desktop/Indiana University/DSX Project 2/Contrasts/LCM vs #SCM/HH_large_cntl_male_vs_small_cntl_male.txt",row.names=FALSE,sep="\t")

###Control_dsx
##Male
#TH
coldata.MTH.CD<-coldata[coldata$Sex=="Male"&coldata$Tissue=="TH"&coldata$Size=="Small",]

data.M<-data.merge[,grep("_Male",colnames(data.merge))]
data.MTH<-data.M[,grep("_TH",colnames(data.M))]
data.MTH.Small<-data.MTH[,grep("sm",colnames(data.MTH))]

ddsMTH.Sm<-DESeqDataSetFromMatrix(
countData=data.MTH.Small,
colData=coldata.MTH.CD,
design= ~Treatment)

dds<-ddsMTH.Sm
dds<- dds[ rowSums(counts(dds)) > 1, ]
dds$Size<-relevel(dds$Treatment,"Cntl")
dds<-DESeq(dds)
res<-results(dds)
res<-res[complete.cases(res$padj),]
res.sig<-res[res$padj<=0.05,]

Gene_ID<-as.data.frame(rownames(res.sig))
colnames(Gene_ID)<-"Gene_ID"
x<-merge(Gene_ID,names,by="Gene_ID")

###Control_dsx
##Male
#Brain
coldata.MBRN.CD<-coldata[coldata$Sex=="Male"&coldata$Tissue=="Brain"&coldata$Size=="Small",]

data.M<-data.merge[,grep("_Male",colnames(data.merge))]
data.MBRN<-data.M[,grep("_Brain",colnames(data.M))]
data.MBRN.Small<-data.MBRN[,grep("sm",colnames(data.MBRN))]

ddsMBRN.Sm<-DESeqDataSetFromMatrix(
countData=data.MBRN.Small,
colData=coldata.MBRN.CD,
design= ~Treatment)

dds<-ddsMBRN.Sm
dds<- dds[ rowSums(counts(dds)) > 1, ]
dds$Size<-relevel(dds$Treatment,"Cntl")
dds<-DESeq(dds)
res<-results(dds)
res<-res[complete.cases(res$padj),]
res.sig<-res[res$padj<=0.05,]

Gene_ID<-as.data.frame(rownames(res.sig))
colnames(Gene_ID)<-"Gene_ID"
x<-merge(Gene_ID,names,by="Gene_ID")

###Control_dsx
##Male
#Gen
coldata.MGEN.CD<-coldata[coldata$Sex=="Male"&coldata$Tissue=="Gen"&coldata$Size=="Small",]

data.M<-data.merge[,grep("_Male",colnames(data.merge))]
data.MGEN<-data.M[,grep("_Gen",colnames(data.M))]
data.MGEN.Small<-data.MGEN[,grep("sm",colnames(data.MGEN))]

ddsMGEN.Sm<-DESeqDataSetFromMatrix(
countData=data.MGEN.Small,
colData=coldata.MGEN.CD,
design= ~Treatment)

dds<-ddsMGEN.Sm
dds<- dds[ rowSums(counts(dds)) > 1, ]
dds$Size<-relevel(dds$Treatment,"Cntl")
dds<-DESeq(dds)
res<-results(dds)
res<-res[complete.cases(res$padj),]
res.sig<-res[res$padj<=0.05,]

Gene_ID<-as.data.frame(rownames(res.sig))
colnames(Gene_ID)<-"Gene_ID"
x<-merge(Gene_ID,names,by="Gene_ID")