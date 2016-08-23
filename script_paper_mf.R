#Basics
#DESeq2 contrasts between control females and control males
data.gen<-read.table(file.choose(),header=T,sep="\t",quote="",)
data.brn<-read.table(file.choose(),header=T,sep="\t",quote="",)
data.hh<-read.table(file.choose(),header=T,sep="\t",quote="",)
data.th<-read.table(file.choose(),header=T,sep="\t",quote="",)

data.gen<-data.gen[complete.cases(data.gen[,13]),]
data.brn<-data.brn[complete.cases(data.brn[,13]),]
data.th<-data.th[complete.cases(data.th[,13]),]
data.hh<-data.hh[complete.cases(data.hh[,13]),]

data.gen$Tiss<-"gen"
data.brn$Tiss<-"brn"
data.hh$Tiss<-"hh"
data.th$Tiss<-"th"

data.all<-rbind(data.gen,data.brn,data.hh,data.th)
data.all$Tiss<-as.factor(data.all$Tiss)
 
data.all$SexBias[data.all$log2FoldChange>=0]<-"F"
data.all$SexBias[data.all$log2FoldChange<=0]<-"M"
data.all$SexBias<-as.factor(data.all$SexBias)

#Sex-biased genes in control animals that are significant at an adjusted pvalue of 0.05

data.all.05<-data.all[data.all$padj<=0.05,]

data.all.05.F<-data.all.05[data.all.05$SexBias=="F",]
data.all.05.M<-data.all.05[data.all.05$SexBias=="M",]

#To arrive at the number of unique significantly female- and male-biased genes (2,720 and 1,565; MS page 9, Line 228), the following code was used.

dim(data.all.05.F[unique(data.all.05.F$Gene.ID),])
dim(data.all.05.M[unique(data.all.05.M$Gene.ID),])

#DESeq2 contrasts between dsxRNAi-injected females and dsxRNAi-injected males (DF-DM)
data.gen.D<-read.table(file.choose(),header=T,sep="\t",quote="",)
data.brn.D<-read.table(file.choose(),header=T,sep="\t",quote="",)
data.hh.D<-read.table(file.choose(),header=T,sep="\t",quote="",)
data.th.D<-read.table(file.choose(),header=T,sep="\t",quote="",)

data.gen.D<-data.gen.D[complete.cases(data.gen.D[,13]),]
data.brn.D<-data.brn.D[complete.cases(data.brn.D[,13]),]
data.th.D<-data.th.D[complete.cases(data.th.D[,13]),]
data.hh.D<-data.hh.D[complete.cases(data.hh.D[,13]),]

data.gen.D$Tiss<-"gen"
data.brn.D$Tiss<-"brn"
data.hh.D$Tiss<-"hh"
data.th.D$Tiss<-"th"

data.all.D<-rbind(data.gen.D,data.brn.D,data.hh.D,data.th.D)
data.all.D$Tiss<-as.factor(data.all.D$Tiss)
 
data.all.D$SexBias[data.all.D$log2FoldChange>=0]<-"F"
data.all.D$SexBias[data.all.D$log2FoldChange<=0]<-"M"
data.all.D$SexBias<-as.factor(data.all.D$SexBias)

#Sex-biased genes in dsx-injected animals that are significant at an adjusted pvalue of 0.05
 
data.all.D.05<-data.all.D[data.all.D$padj<=0.05,]

data.all.D.05.F<-data.all.D.05[data.all.D.05$SexBias=="F",]
data.all.D.05.M<-data.all.D.05[data.all.D.05$SexBias=="M",]

#DESeq2 contrasts between control males and dsx males (CM-DM)
data.gen.dmcm<-read.table(file.choose(),header=T,sep="\t",quote="")
data.brn.dmcm<-read.table(file.choose(),header=T,sep="\t",quote="")
data.hh.dmcm<-read.table(file.choose(),header=T,sep="\t",quote="")
data.th.dmcm<-read.table(file.choose(),header=T,sep="\t",quote="")

data.gen.dmcm<-data.gen.dmcm[complete.cases(data.gen.dmcm[,13]),]
data.brn.dmcm<-data.brn.dmcm[complete.cases(data.brn.dmcm[,13]),]
data.th.dmcm<-data.th.dmcm[complete.cases(data.th.dmcm[,13]),]
data.hh.dmcm<-data.hh.dmcm[complete.cases(data.hh.dmcm[,13]),]

data.gen.dmcm$Tiss<-"gen"
data.brn.dmcm$Tiss<-"brn"
data.hh.dmcm$Tiss<-"hh"
data.th.dmcm$Tiss<-"th"

data.all.dmcm<-rbind(data.gen.dmcm,data.brn.dmcm,data.hh.dmcm,data.th.dmcm)
data.all.dmcm$Tiss<-as.factor(data.all.dmcm$Tiss)

data.all.dmcm$Bias[data.all.dmcm$log2FoldChange>=0]<-"D"
data.all.dmcm$Bias[data.all.dmcm$log2FoldChange<=0]<-"C"
data.all.dmcm$Bias<-as.factor(data.all.dmcm$Bias)

data.all.dmcm.05<-data.all.dmcm[data.all.dmcm$padj<=0.05,]

#DESeq2 contrasts between control females and dsx females (CF-DF)
data.gen.dfcf<-read.table(file.choose(),header=T,sep="\t",quote="")
data.brn.dfcf<-read.table(file.choose(),header=T,sep="\t",quote="")
data.hh.dfcf<-read.table(file.choose(),header=T,sep="\t",quote="")
data.th.dfcf<-read.table(file.choose(),header=T,sep="\t",quote="")

data.gen.dfcf<-data.gen.dfcf[complete.cases(data.gen.dfcf[,13]),]
data.brn.dfcf<-data.brn.dfcf[complete.cases(data.brn.dfcf[,13]),]
data.th.dfcf<-data.th.dfcf[complete.cases(data.th.dfcf[,13]),]
data.hh.dfcf<-data.hh.dfcf[complete.cases(data.hh.dfcf[,13]),]

data.gen.dfcf$Tiss<-"gen"
data.brn.dfcf$Tiss<-"brn"
data.hh.dfcf$Tiss<-"hh"
data.th.dfcf$Tiss<-"th"

data.all.dfcf<-rbind(data.gen.dfcf,data.brn.dfcf,data.hh.dfcf,data.th.dfcf)
data.all.dfcf$Tiss<-as.factor(data.all.dfcf$Tiss)

data.all.dfcf$Bias[data.all.dfcf$log2FoldChange>=0]<-"D"
data.all.dfcf$Bias[data.all.dfcf$log2FoldChange<=0]<-"C"
data.all.dfcf$Bias<-as.factor(data.all.dfcf$Bias)

#dsx-mediated genes in female animals that are significant at an adjusted pvalue of 0.05

data.all.dfcf.05<-data.all.dfcf[data.all.dfcf$padj<=0.05,]

#Figure 2, Panel A, Fold Change Regression
data.all.D.F<-data.all.D[data.all.D$SexBias=="F",]
data.all.D.M<-data.all.D[data.all.D$SexBias=="M",]

data.all.F<-data.all[data.all$SexBias=="F",]
data.all.M<-data.all[data.all$SexBias=="M",]

data.same.F<-merge(data.all.F,data.all.D.F,by.x=c("Gene.ID","Tiss"),by.y=c("Gene.ID","Tiss"),all=FALSE)
data.same.M<-merge(data.all.M,data.all.D.M,by.x=c("Gene.ID","Tiss"),by.y=c("Gene.ID","Tiss"),all=FALSE)

data.same.05.F<-merge(data.all.05.F,data.all.D.F,by.x=c("Gene.ID","Tiss"),by.y=c("Gene.ID","Tiss"),all=FALSE)
data.same.05.M<-merge(data.all.05.M,data.all.D.M,by.x=c("Gene.ID","Tiss"),by.y=c("Gene.ID","Tiss"),all=FALSE)

par(oma=c(0,1,0,0))
plot(log2FoldChange.y~log2FoldChange.x,data=data.same.F, pch=20,bg="gray90",col="gray80",xlab="",ylab="",main="",ylim=c(0,5),xlim=c(0,5),type="n")
points(log2FoldChange.y~log2FoldChange.x,data=data.same.F,pch=20,bg="gray90",col="gray80",cex=1)
points(abs(log2FoldChange.y)~abs(log2FoldChange.x),data=data.same.M, pch=20, bg = "gray90",col="gray80",cex=1)
abline(a=0,b=1,col="blue")
points(abs(log2FoldChange.y)~abs(log2FoldChange.x),data=data.same.05.F, pch=20, col = rgb(red=1,green=0.25, blue=0.2,alpha=0.4),cex=1)
points(abs(log2FoldChange.y)~abs(log2FoldChange.x),data=data.same.05.M, pch=20, col = rgb(red=0,green=0.6, blue=.8,alpha=0.4),cex=1)
mtext("LogFC Controls",1,cex=1.3,line=3,font=2)
mtext("LogFC DSX Knockdown",2,cex=1.3,line=3,font=2)


#Figure 2, Panel B, sex-biased in control animals
data.gen.05<-data.all.05[data.all.05$Tiss=="gen",]
data.brn.05<-data.all.05[data.all.05$Tiss=="brn",]
data.th.05<-data.all.05[data.all.05$Tiss=="th",]
data.hh.05<-data.all.05[data.all.05$Tiss=="hh",]

data.gen.05F<-data.gen.05[data.gen.05$SexBias=="F",]
data.gen.05M<-data.gen.05[data.gen.05$SexBias=="M",]

data.brn.05F<-data.brn.05[data.brn.05$SexBias=="F",]
data.brn.05M<-data.brn.05[data.brn.05$SexBias=="M",]

data.th.05F<-data.th.05[data.th.05$SexBias=="F",]
data.th.05M<-data.th.05[data.th.05$SexBias=="M",]

data.hh.05F<-data.hh.05[data.hh.05$SexBias=="F",]
data.hh.05M<-data.hh.05[data.hh.05$SexBias=="M",]

gen.bar.05<-c(dim(data.gen.05F)[1],dim(data.gen.05M)[1])
brn.bar.05<-c(dim(data.brn.05F)[1],dim(data.brn.05M)[1])
th.bar.05<-c(dim(data.th.05F)[1],dim(data.th.05M)[1])
hh.bar.05<-c(dim(data.hh.05F)[1],dim(data.hh.05M)[1])

data.05.bar<-rbind(brn.bar.05,gen.bar.05,th.bar.05,hh.bar.05)
data.05.bar<-as.data.frame(t(data.05.bar))
colnames(data.05.bar)<-c("brain","genitalia","thor_horn","head_horn")
rownames(data.05.bar)<-c("female","male")

par(oma=c(0,1,0,0))
barplot(as.matrix(data.05.bar),beside=TRUE,main="",ylim=c(0,3500),cex.main=1.3,xaxt="n",col=c("indianred","lightskyblue"))
box()
legend("topright", inset=0.03, c("Female-biased","Male-biased"), cex=2, fill=c("indianred","lightskyblue"))
mtext("Biased Genes",2,cex=1.7,line=3.2,font=2)
mtext("brn           gen             th             hh",1,cex=1.7,line=2,font=2)
mtext("control animals",3,cex=1.7,font=2, line=2)

#Figure 2, Panel B, sex-biased as assessed in control animals, in dsxRNAi background
#Finding what genes that were originally sex-biased in control animals are still sex-biased in dsxRNAi animals
data.diff.05.F<-merge(data.all.05.F,data.all.D.05.F,by.x=c("Gene.ID","Tiss"),by.y=c("Gene.ID","Tiss"),all=FALSE)
data.diff.05.M<-merge(data.all.05.M,data.all.D.05.M,by.x=c("Gene.ID","Tiss"),by.y=c("Gene.ID","Tiss"),all=FALSE)

data.diff.05.F.gen<-data.diff.05.F[data.diff.05.F$Tiss=="gen",]
data.diff.05.F.brn<-data.diff.05.F[data.diff.05.F$Tiss=="brn",]
data.diff.05.F.th<-data.diff.05.F[data.diff.05.F$Tiss=="th",]
data.diff.05.F.hh<-data.diff.05.F[data.diff.05.F$Tiss=="hh",]
data.diff.05.M.gen<-data.diff.05.M[data.diff.05.M$Tiss=="gen",]
data.diff.05.M.brn<-data.diff.05.M[data.diff.05.M$Tiss=="brn",]
data.diff.05.M.th<-data.diff.05.M[data.diff.05.M$Tiss=="th",]
data.diff.05.M.hh<-data.diff.05.M[data.diff.05.M$Tiss=="hh",]
gen.diff.bar.05<-c(dim(data.diff.05.F.gen)[1],dim(data.diff.05.M.gen)[1])
brn.diff.bar.05<-c(dim(data.diff.05.F.brn)[1],dim(data.diff.05.M.brn)[1])
th.diff.bar.05<-c(dim(data.diff.05.F.th)[1],dim(data.diff.05.M.th)[1])
hh.diff.bar.05<-c(dim(data.diff.05.F.hh)[1],dim(data.diff.05.M.hh)[1])
data.diff.05.bar<-rbind(brn.diff.bar.05,gen.diff.bar.05,th.diff.bar.05,hh.diff.bar.05)

data.diff.05.bar<-as.data.frame(t(data.diff.05.bar))

par(oma=c(0,1,0,0))
barplot(as.matrix(data.diff.05.bar),beside=TRUE,main="",ylim=c(0,3500),cex.main=1.3,xaxt="n",col=c("indianred","lightskyblue"),density=30,angle=c(135,45))
box()
legend("topright", inset=0.03, c("Female-biased","Male-biased"), cex=2, fill=c("indianred","lightskyblue"),angle=c(135,45),density=30)
mtext("Biased Genes",2,cex=1.7,line=3.2,font=2)
mtext("brn           gen             th             hh",1,cex=1.7,line=2,font=2)
mtext("dsx-knockdown animals",3,cex=1.7,font=2,line=2)

#To arrive at the number of genes that were sex-biased in control animals and still sex-biased in dsxRNAi-injected animal, within homologous tissues, the following code was used (resulting numbers should be 1093 female-biased genes and 503 male-biased genes for a total of 1596 sex-biased genes. Dividing the number of sex-biased genes as assessed in control animals - 4285 - by 1596 reveals a 2.7-fold change; Caption Figure 2B, MS page 11)

dim(data.diff.05.F[unique(data.diff.05.F$Gene.ID),])
dim(data.diff.05.M[unique(data.diff.05.M$Gene.ID),])

#Figure 2, Panel C, dsx-mediated in females
data.gen.dfcf.05<-data.all.dfcf.05[data.all.dfcf.05$Tiss=="gen",]
data.brn.dfcf.05<-data.all.dfcf.05[data.all.dfcf.05$Tiss=="brn",]
data.th.dfcf.05<-data.all.dfcf.05[data.all.dfcf.05$Tiss=="th",]
data.hh.dfcf.05<-data.all.dfcf.05[data.all.dfcf.05$Tiss=="hh",]

data.gen.dfcf.05D<-data.gen.dfcf.05[data.gen.dfcf.05$Bias=="D",]
data.gen.dfcf.05C<-data.gen.dfcf.05[data.gen.dfcf.05$Bias=="C",]

data.brn.dfcf.05D<-data.brn.dfcf.05[data.brn.dfcf.05$Bias=="D",]
data.brn.dfcf.05C<-data.brn.dfcf.05[data.brn.dfcf.05$Bias=="C",]

data.th.dfcf.05D<-data.th.dfcf.05[data.th.dfcf.05$Bias=="D",]
data.th.dfcf.05C<-data.th.dfcf.05[data.th.dfcf.05$Bias=="C",]

data.hh.dfcf.05D<-data.hh.dfcf.05[data.hh.dfcf.05$Bias=="D",]
data.hh.dfcf.05C<-data.hh.dfcf.05[data.hh.dfcf.05$Bias=="C",]

gen.bar.dfcf.05<-c(dim(data.gen.dfcf.05D)[1],dim(data.gen.dfcf.05C)[1])
brn.bar.dfcf.05<-c(dim(data.brn.dfcf.05D)[1],dim(data.brn.dfcf.05C)[1])
th.bar.dfcf.05<-c(dim(data.th.dfcf.05D)[1],dim(data.th.dfcf.05C)[1])
hh.bar.dfcf.05<-c(dim(data.hh.dfcf.05D)[1],dim(data.hh.dfcf.05C)[1])

data.dfcf.05.bar<-rbind(brn.bar.dfcf.05,gen.bar.dfcf.05,th.bar.dfcf.05,hh.bar.dfcf.05)
data.dfcf.05.bar<-as.data.frame(t(data.dfcf.05.bar))
colnames(data.dfcf.05.bar)<-c("brain","genitalia","thor_horn","head_horn")
rownames(data.dfcf.05.bar)<-c("DSXrnai","Control")

par(oma=c(0,1,0,0))
barplot(as.matrix(data.dfcf.05.bar),beside=TRUE,main="",ylim=c(0,700),cex.main=1.3,xaxt="n",density=c(30,110),angle=c(135,0),col="indianred")
box()
legend("topright", inset=0.03, c("RNAi","Control"), cex=2, density=c(30,110),angle=c(135,0),fill="indianred")
mtext("Biased Genes",2,cex=1.7,line=3.2,font=2)
mtext("brn           gen             th             hh",1,cex=1.7,line=2,font=2)
mtext("dsx-mediated in females",3,cex=1.7,font=2, line=2)

#Figure 2, Panel C, dsx-mediated in males
data.gen.dmcm.05<-data.all.dmcm.05[data.all.dmcm.05$Tiss=="gen",]
data.brn.dmcm.05<-data.all.dmcm.05[data.all.dmcm.05$Tiss=="brn",]
data.th.dmcm.05<-data.all.dmcm.05[data.all.dmcm.05$Tiss=="th",]
data.hh.dmcm.05<-data.all.dmcm.05[data.all.dmcm.05$Tiss=="hh",]

data.gen.dmcm.05D<-data.gen.dmcm.05[data.gen.dmcm.05$Bias=="D",]
data.gen.dmcm.05C<-data.gen.dmcm.05[data.gen.dmcm.05$Bias=="C",]

data.brn.dmcm.05D<-data.brn.dmcm.05[data.brn.dmcm.05$Bias=="D",]
data.brn.dmcm.05C<-data.brn.dmcm.05[data.brn.dmcm.05$Bias=="C",]

data.th.dmcm.05D<-data.th.dmcm.05[data.th.dmcm.05$Bias=="D",]
data.th.dmcm.05C<-data.th.dmcm.05[data.th.dmcm.05$Bias=="C",]

data.hh.dmcm.05D<-data.hh.dmcm.05[data.hh.dmcm.05$Bias=="D",]
data.hh.dmcm.05C<-data.hh.dmcm.05[data.hh.dmcm.05$Bias=="C",]

gen.bar.dmcm.05<-c(dim(data.gen.dmcm.05D)[1],dim(data.gen.dmcm.05C)[1])
brn.bar.dmcm.05<-c(dim(data.brn.dmcm.05D)[1],dim(data.brn.dmcm.05C)[1])
th.bar.dmcm.05<-c(dim(data.th.dmcm.05D)[1],dim(data.th.dmcm.05C)[1])
hh.bar.dmcm.05<-c(dim(data.hh.dmcm.05D)[1],dim(data.hh.dmcm.05C)[1])

data.dmcm.05.bar<-rbind(brn.bar.dmcm.05,gen.bar.dmcm.05,th.bar.dmcm.05,hh.bar.dmcm.05)
data.dmcm.05.bar<-as.data.frame(t(data.dmcm.05.bar))
colnames(data.dmcm.05.bar)<-c("brain","genitalia","thor_horn","head_horn")
rownames(data.dmcm.05.bar)<-c("RNAi","Control")

par(oma=c(0,1,0,0))
barplot(as.matrix(data.dmcm.05.bar),beside=TRUE,main="",ylim=c(0,700),cex.main=1.3,xaxt="n",density=c(30,110),angle=c(45,0),col="lightskyblue")
box()
legend("topright", inset=0.03, c("RNAi","Control"), cex=2, density=c(30,110),angle=c(45,0),fill="lightskyblue")
mtext("Biased Genes",2,cex=1.7,line=3.2,font=2)
mtext("brn           gen             th             hh",1,cex=1.7,line=2,font=2)
mtext("dsx-mediated in males",3,cex=1.7,font=2, line=2)

#Supplementary Figure 1, distribution of sex-biased genes across tissues
#Finding how many times a given gene is duplicated (i.e., in how many tissues does it exhibit biased expression?)

data.all.05.F$numdup<-as.numeric("1")
data.all.05.M$numdup<-as.numeric("1")

spec.tab.F<-aggregate(numdup~Gene.ID,data=data.all.05.F,FUN=sum)
spec.tab.M<-aggregate(numdup~Gene.ID,data=data.all.05.M,FUN=sum)

Tiss.1.F.frame<-spec.tab.F[spec.tab.F$numdup=="1",]
Tiss.1.M.frame<-spec.tab.M[spec.tab.M$numdup=="1",]

Tiss.2.F.frame<-spec.tab.F[spec.tab.F$numdup=="2",]
Tiss.2.M.frame<-spec.tab.M[spec.tab.M$numdup=="2",]

Tiss.3.F.frame<-spec.tab.F[spec.tab.F$numdup=="3",]
Tiss.3.M.frame<-spec.tab.M[spec.tab.M$numdup=="3",]

Tiss.4.F.frame<-spec.tab.F[spec.tab.F$numdup=="4",]
Tiss.4.M.frame<-spec.tab.M[spec.tab.M$numdup=="4",]

#Recovering tissue information after counting incidences of biased-gene expression
Tiss.1.F.test<-merge(Tiss.1.F.frame,data.all.05.F,by.x="Gene.ID",by.y="Gene.ID")
Tiss.1.M.test<-merge(Tiss.1.M.frame,data.all.05.M,by.x="Gene.ID",by.y="Gene.ID")

Tiss.2.F.test<-merge(Tiss.2.F.frame,data.all.05.F,by.x="Gene.ID",by.y="Gene.ID")
Tiss.2.M.test<-merge(Tiss.2.M.frame,data.all.05.M,by.x="Gene.ID",by.y="Gene.ID")

Tiss.3.F.test<-merge(Tiss.3.F.frame,data.all.05.F,by.x="Gene.ID",by.y="Gene.ID")
Tiss.3.M.test<-merge(Tiss.3.M.frame,data.all.05.M,by.x="Gene.ID",by.y="Gene.ID")

Tiss.4.F.test<-merge(Tiss.4.F.frame,data.all.05.F,by.x="Gene.ID",by.y="Gene.ID")
Tiss.4.M.test<-merge(Tiss.4.M.frame,data.all.05.M,by.x="Gene.ID",by.y="Gene.ID")

#Separating by tissue
Tiss.2.F.H<-Tiss.2.F.test[Tiss.2.F.test$Tiss=="hh",]
Tiss.2.F.T<-Tiss.2.F.test[Tiss.2.F.test$Tiss=="th",]
Tiss.2.F.G<-Tiss.2.F.test[Tiss.2.F.test$Tiss=="gen",]
Tiss.2.F.B<-Tiss.2.F.test[Tiss.2.F.test$Tiss=="brn",]

Tiss.2.M.H<-Tiss.2.M.test[Tiss.2.M.test$Tiss=="hh",]
Tiss.2.M.T<-Tiss.2.M.test[Tiss.2.M.test$Tiss=="th",]
Tiss.2.M.G<-Tiss.2.M.test[Tiss.2.M.test$Tiss=="gen",]
Tiss.2.M.B<-Tiss.2.M.test[Tiss.2.M.test$Tiss=="brn",]

Tiss.3.F.H<-Tiss.3.F.test[Tiss.3.F.test$Tiss=="hh",]
Tiss.3.F.T<-Tiss.3.F.test[Tiss.3.F.test$Tiss=="th",]
Tiss.3.F.G<-Tiss.3.F.test[Tiss.3.F.test$Tiss=="gen",]
Tiss.3.F.B<-Tiss.3.F.test[Tiss.3.F.test$Tiss=="brn",]

Tiss.3.M.H<-Tiss.3.M.test[Tiss.3.M.test$Tiss=="hh",]
Tiss.3.M.T<-Tiss.3.M.test[Tiss.3.M.test$Tiss=="th",]
Tiss.3.M.G<-Tiss.3.M.test[Tiss.3.M.test$Tiss=="gen",]
Tiss.3.M.B<-Tiss.3.M.test[Tiss.3.M.test$Tiss=="brn",]

Tiss.4.F.H<-Tiss.4.F.test[Tiss.3.F.test$Tiss=="hh",] 
Tiss.4.F.T<-Tiss.4.F.test[Tiss.3.F.test$Tiss=="th",]
Tiss.4.F.G<-Tiss.4.F.test[Tiss.3.F.test$Tiss=="gen",]
Tiss.4.F.B<-Tiss.4.F.test[Tiss.3.F.test$Tiss=="brn",]

Tiss.4.M.H<-Tiss.4.M.test[Tiss.4.M.test$Tiss=="hh",]
Tiss.4.M.T<-Tiss.4.M.test[Tiss.4.M.test$Tiss=="th",]
Tiss.4.M.G<-Tiss.4.M.test[Tiss.4.M.test$Tiss=="gen",]
Tiss.4.M.B<-Tiss.4.M.test[Tiss.4.M.test$Tiss=="brn",]

#The number of sex-biased genes found exclusively head horn tissue
dim(Tiss.1.F.test[Tiss.1.F.test$Tiss=="hh",])[1]
dim(Tiss.1.M.test[Tiss.1.M.test$Tiss=="hh",])[1]

#The number of sex-biased genes found exclusively thoracic horn tissue
dim(Tiss.1.F.test[Tiss.1.F.test$Tiss=="th",])[1]
dim(Tiss.1.M.test[Tiss.1.M.test$Tiss=="th",])[1]

#The number of sex-biased genes found exclusively brain tissue
dim(Tiss.1.F.test[Tiss.1.F.test$Tiss=="brn",])[1]
dim(Tiss.1.M.test[Tiss.1.M.test$Tiss=="brn",])[1]

#The number of sex-biased genes found exclusively genitalia tissue
dim(Tiss.1.F.test[Tiss.1.F.test$Tiss=="gen",])[1]
dim(Tiss.1.M.test[Tiss.1.M.test$Tiss=="gen",])[1]

#The number of sex-biased genes found across head horn/thoracic horn tissue
dim(merge(Tiss.2.F.H,Tiss.2.F.T,by="Gene.ID"))[1]
dim(merge(Tiss.2.M.H,Tiss.2.M.T,by="Gene.ID"))[1]

#The number of sex-biased genes found across head horn/genitalia tissue
dim(merge(Tiss.2.F.H,Tiss.2.F.G,by="Gene.ID"))[1]
dim(merge(Tiss.2.M.H,Tiss.2.M.G,by="Gene.ID"))[1]

#The number of sex-biased genes found across thoracic horn/genitalia tissue
dim(merge(Tiss.2.F.T,Tiss.2.F.G,by="Gene.ID"))[1]
dim(merge(Tiss.2.M.T,Tiss.2.M.G,by="Gene.ID"))[1]

#The number of sex-biased genes found across thead horn/brain tissue
dim(merge(Tiss.2.F.H,Tiss.2.F.B,by="Gene.ID"))[1]
dim(merge(Tiss.2.M.H,Tiss.2.M.B,by="Gene.ID"))[1]

#The number of sex-biased genes found across thoracic horn/brain tissue
dim(merge(Tiss.2.F.T,Tiss.2.F.B,by="Gene.ID"))[1]
dim(merge(Tiss.2.M.T,Tiss.2.M.B,by="Gene.ID"))[1]

#The number of sex-biased genes found across genitalia/brain tissue
dim(merge(Tiss.2.F.G,Tiss.2.F.B,by="Gene.ID"))[1]
dim(merge(Tiss.2.M.G,Tiss.2.M.B,by="Gene.ID"))[1]

#The number of sex-biased genes found across head horn/thoracic horn/genitalia tissue
Tiss.3.F.HT<-merge(Tiss.3.F.H,Tiss.3.F.T,by="Gene.ID")
dim(merge(Tiss.3.F.HT,Tiss.3.F.G,by="Gene.ID"))[1]

Tiss.3.M.HT<-merge(Tiss.3.M.H,Tiss.3.M.T,by="Gene.ID")
dim(merge(Tiss.3.M.HT,Tiss.3.M.G,by="Gene.ID"))[1]

#The number of sex-biased genes found across head horn/thoracic horn/brain tissue
dim(merge(Tiss.3.F.HT,Tiss.3.F.B,by="Gene.ID"))[1]
dim(merge(Tiss.3.M.HT,Tiss.3.M.B,by="Gene.ID"))[1]

#The number of sex-biased genes found across head horn/genitalia/brain tissue
Tiss.3.F.HG<-merge(Tiss.3.F.H,Tiss.3.F.G,by="Gene.ID")
dim(merge(Tiss.3.F.HG,Tiss.3.F.B,by="Gene.ID"))[1]

Tiss.3.M.HG<-merge(Tiss.3.M.H,Tiss.3.M.G,by="Gene.ID")
dim(merge(Tiss.3.M.HG,Tiss.3.M.B,by="Gene.ID"))[1]

#The number of sex-biased genes found across brain/thoracic horn/genitalia tissue
Tiss.3.F.BT<-merge(Tiss.3.F.B,Tiss.3.F.T,by="Gene.ID")
dim(merge(Tiss.3.F.BT,Tiss.3.F.G,by="Gene.ID"))[1]

Tiss.3.M.BT<-merge(Tiss.3.M.B,Tiss.3.M.T,by="Gene.ID")
dim(merge(Tiss.3.M.BT,Tiss.3.M.G,by="Gene.ID"))[1]

#The number of sex-biased genes found across all tissues
dim(Tiss.4.F.frame)[1]
dim(Tiss.4.M.frame)[1]

#Figure 3 and Supplementary Figure 2, Changes in sex-biased gene expression across different treatment groups
#Use with file "OnthophagusFPKM.annot.tsv"
counts<-read.table(file.choose(),header=T,sep="\t",quote="")
counts<-counts[,c(1,8:103)]

counts$CMB<-rowMeans(counts[2:7],na.rm=TRUE)
counts$CMT<-rowMeans(counts[8:13],na.rm=TRUE)
counts$CMH<-rowMeans(counts[14:19],na.rm=TRUE)
counts$CMG<-rowMeans(counts[20:25],na.rm=TRUE)

counts$CFB<-rowMeans(counts[26:31],na.rm=TRUE)
counts$CFT<-rowMeans(counts[32:37],na.rm=TRUE)
counts$CFH<-rowMeans(counts[38:43],na.rm=TRUE)
counts$CFG<-rowMeans(counts[44:49],na.rm=TRUE)

counts$DMB<-rowMeans(counts[50:55],na.rm=TRUE)
counts$DMT<-rowMeans(counts[56:61],na.rm=TRUE)
counts$DMH<-rowMeans(counts[62:67],na.rm=TRUE)
counts$DMG<-rowMeans(counts[68:73],na.rm=TRUE)

counts$DFB<-rowMeans(counts[74:79],na.rm=TRUE)
counts$DFT<-rowMeans(counts[80:85],na.rm=TRUE)
counts$DFH<-rowMeans(counts[86:91],na.rm=TRUE)
counts$DFG<-rowMeans(counts[92:97],na.rm=TRUE)

counts.new<-counts[,c(1,98:113)]

counts.hh<-counts.new[,c(1,4,8,12,16)]
counts.th<-counts.new[,c(1,3,7,11,15)]
counts.gen<-counts.new[,c(1,5,9,13,17)]
counts.brn<-counts.new[,c(1,2,6,10,14)]

#Merging information about significant male and female-biased genes over tissues (determined from DESeq2 contrasts) with expression data (FPKM)

M.biased.gen<-merge(data.gen.05M,counts.gen,by.x="Gene.ID",by.y="Gene_ID",all=FALSE)
group<-c(rep(c("CMG","CFG","DMG","DFG"),each=dim(M.biased.gen)[1]))

M.biased.brn<-merge(data.brn.05M,counts.brn,by.x="Gene.ID",by.y="Gene_ID",all=FALSE)
group<-c(rep(c("CMG","CFG","DMG","DFG"),each=dim(M.biased.brn)[1]))

M.biased.hh<-merge(data.hh.05M,counts.hh,by.x="Gene.ID",by.y="Gene_ID",all=FALSE)
group<-c(rep(c("CMG","CFG","DMG","DFG"),each=dim(M.biased.hh)[1]))

M.biased.th<-merge(data.th.05M,counts.th,by.x="Gene.ID",by.y="Gene_ID",all=FALSE)
group<-c(rep(c("CMG","CFG","DMG","DFG"),each=dim(M.biased.th)[1]))

F.biased.gen<-merge(data.gen.05F,counts.gen,by.x="Gene.ID",by.y="Gene_ID",all=FALSE)
group<-c(rep(c("CMG","CFG","DMG","DFG"),each=dim(F.biased.gen)[1]))

F.biased.brn<-merge(data.brn.05F,counts.brn,by.x="Gene.ID",by.y="Gene_ID",all=FALSE)
group<-c(rep(c("CMG","CFG","DMG","DFG"),each=dim(F.biased.brn)[1]))

F.biased.hh<-merge(data.hh.05F,counts.hh,by.x="Gene.ID",by.y="Gene_ID",all=FALSE)
group<-c(rep(c("CMG","CFG","DMG","DFG"),each=dim(F.biased.hh)[1]))

F.biased.th<-merge(data.th.05F,counts.th,by.x="Gene.ID",by.y="Gene_ID",all=FALSE)
group<-c(rep(c("CMG","CFG","DMG","DFG"),each=dim(F.biased.th)[1]))

#Graphing male-biased genes in genitalia over treatment groups (Supplementary Figure 2, Top Panel, Left)
group<-c(rep(c("CMG","CFG","DMG","DFG"),each=dim(M.biased.gen)[1]))
group<-as.factor(group)
fpkm.gen<-c(M.biased.gen$CMG,M.biased.gen$CFG,M.biased.gen$DMG,M.biased.gen$DFG)
M.bias.gen<-data.frame(group,fpkm.gen)
M.bias.gen$group<-factor(M.bias.gen$group,levels=c("CMG","DMG","DFG","CFG"))

par(oma=c(0,1,0,0))
plot(log(fpkm.gen)~group,data=M.bias.gen,main="",ylab="",col=c("royalblue4","skyblue1","indianred1","indianred4"))
mtext("Average Log Expression (FPKM)",2,cex=1.3,line=3,font=2)
mtext("Male-biased in Genitalia, padj <= 0.05",3,cex=1.3,line=1,font=2)

#Graphing male-biased genes in head horns over treatment groups (Figure 3B, Left)
group<-c(rep(c("CMH","CFH","DMH","DFH"),each=dim(M.biased.hh)[1]))
group<-as.factor(group)
fpkm.hh<-c(M.biased.hh$CMH,M.biased.hh$CFH,M.biased.hh$DMH,M.biased.hh$DFH)
M.bias.hh<-data.frame(group,fpkm.hh)
M.bias.hh$group<-factor(M.bias.hh$group,levels=c("CMH","DMH","DFH","CFH"))

par(oma=c(0,1,0,0))
plot(log(fpkm.hh)~group,data=M.bias.hh,main="",ylab="",col=c("royalblue4","skyblue1","indianred1","indianred4"))
mtext("Average Log Expression (FPKM)",2,cex=1.3,line=3,font=2)
mtext("Male-biased in Head Horns, padj <= 0.05",3,cex=1.3,line=1,font=2)

#Graphing male-biased genes in thoracic horns over treatment groups (Supplementary Figure 2, Bottom Panel, Left)
group<-c(rep(c("CMT","CFT","DMT","DFT"),each=dim(M.biased.th)[1]))
group<-as.factor(group)
fpkm.th<-c(M.biased.th$CMT,M.biased.th$CFT,M.biased.th$DMT,M.biased.th$DFT)
M.bias.th<-data.frame(group,fpkm.th)
M.bias.th$group<-factor(M.bias.th$group,levels=c("CMT","DMT","DFT","CFT"))

par(oma=c(0,1,0,0))
plot(log(fpkm.th)~group,data=M.bias.th,main="",ylab="",col=c("royalblue4","skyblue1","indianred1","indianred4"))
mtext("Average Log Expression (FPKM)",2,cex=1.3,line=3,font=2)
mtext("Male-biased in Thoracic Horns, padj <= 0.05",3,cex=1.3,line=1,font=2)

#Graphing male-biased genes in brain over treatment groups (Figure 3A, Left)
group<-c(rep(c("CMB","CFB","DMB","DFB"),each=dim(M.biased.brn)[1]))
group<-as.factor(group)
fpkm.brn<-c(M.biased.brn$CMB,M.biased.brn$CFB,M.biased.brn$DMB,M.biased.brn$DFB)
M.bias.brn<-data.frame(group,fpkm.brn)
M.bias.brn$group<-factor(M.bias.brn$group,levels=c("CMB","DMB","DFB","CFB"))

par(oma=c(0,1,0,0))
plot(log(fpkm.brn)~group,data=M.bias.brn,main="",ylab="",col=c("royalblue4","skyblue1","indianred1","indianred4"))
mtext("Average Log Expression (FPKM)",2,cex=1.3,line=3,font=2)
mtext("Male-biased in Brains, padj <= 0.05",3,cex=1.3,line=1,font=2)

#Graphing female-biased genes in genitalia over treatment groups (Supplementary Figure 2, Top Panel, Right)
group<-c(rep(c("CMG","CFG","DMG","DFG"),each=dim(F.biased.gen)[1]))
group<-as.factor(group)
fpkm.genF<-c(F.biased.gen$CMG,F.biased.gen$CFG,F.biased.gen$DMG,F.biased.gen$DFG)
F.bias.gen<-data.frame(group,fpkm.genF)
F.bias.gen$group<-factor(F.bias.gen$group,levels=c("CMG","DMG","DFG","CFG"))

par(oma=c(0,1,0,0))
plot(log(fpkm.genF)~group,data=F.bias.gen,main="",ylab="",col=c("royalblue4","skyblue1","indianred1","indianred4"))
mtext("Average Log Expression (FPKM)",2,cex=1.3,line=3,font=2)
mtext("Female-biased in Genitalia, padj <= 0.05",3,cex=1.3,line=1,font=2)

#Graphing female-biased genes in head horns over treatment groups (Figure 3B, Right)
group<-c(rep(c("CMH","CFH","DMH","DFH"),each=dim(F.biased.hh)[1]))
group<-as.factor(group)
fpkm.hhF<-c(F.biased.hh$CMH,F.biased.hh$CFH,F.biased.hh$DMH,F.biased.hh$DFH)
F.bias.hh<-data.frame(group,fpkm.hhF)
F.bias.hh$group<-factor(F.bias.hh$group,levels=c("CMH","DMH","DFH","CFH"))

par(oma=c(0,1,0,0))
plot(log(fpkm.hhF)~group,data=F.bias.hh,main="",ylab="",col=c("royalblue4","skyblue1","indianred1","indianred4"))
mtext("Average Log Expression (FPKM)",2,cex=1.3,line=3,font=2)
mtext("Female-biased in Head Horns, padj <= 0.05",3,cex=1.3,line=1,font=2)

#Graphing female-biased genes in thoracic horns over treatment groups (Supplementary Figure 2, Bottom Panel, Right)
group<-c(rep(c("CMT","CFT","DMT","DFT"),each=dim(F.biased.th)[1]))
group<-as.factor(group)
fpkm.thF<-c(F.biased.th$CMT,F.biased.th$CFT,F.biased.th$DMT,F.biased.th$DFT)
F.bias.th<-data.frame(group,fpkm.thF)
F.bias.th$group<-factor(F.bias.th$group,levels=c("CMT","DMT","DFT","CFT"))

par(oma=c(0,1,0,0))
plot(log(fpkm.thF)~group,data=F.bias.th,main="",ylab="",col=c("royalblue4","skyblue1","indianred1","indianred4"))
mtext("Average Log Expression (FPKM)",2,cex=1.3,line=3,font=2)
mtext("Female-biased in Thoracic Horns, padj <= 0.05",3,cex=1.3,line=1,font=2)

#Graphing female-biased genes in brain over treatment groups (Figure 3A, Right)
group<-c(rep(c("CMB","CFB","DMB","DFB"),each=dim(F.biased.brn)[1]))
group<-as.factor(group)
fpkm.brnF<-c(F.biased.brn$CMB,F.biased.brn$CFB,F.biased.brn$DMB,F.biased.brn$DFB)
F.bias.brn<-data.frame(group,fpkm.brnF)
F.bias.brn$group<-factor(F.bias.brn$group,levels=c("CMB","DMB","DFB","CFB"))

par(oma=c(0,1,0,0))
plot(log(fpkm.brnF)~group,data=F.bias.brn,main="",ylab="",col=c("royalblue4","skyblue1","indianred1","indianred4"))
mtext("Average Log Expression (FPKM)",2,cex=1.3,line=3,font=2)
mtext("Female-biased in Brains, padj <= 0.05",3,cex=1.3,line=1,font=2)

#Figure 3, Statistics
#Using Wilcoxon test on distributions (which are non-normal) and also apply Benjamini Hochberg correction for multiple comparisons.

#Male brains
M.bias.brn.CMDM<-M.bias.brn[M.bias.brn$group=="CMB"|M.bias.brn$group=="DMB",]
M.bias.brn.CMDF<-M.bias.brn[M.bias.brn$group=="CMB"|M.bias.brn$group=="DFB",]
M.bias.brn.CMCF<-M.bias.brn[M.bias.brn$group=="CMB"|M.bias.brn$group=="CFB",]
M.bias.brn.DMDF<-M.bias.brn[M.bias.brn$group=="DMB"|M.bias.brn$group=="DFB",]
M.bias.brn.DMCF<-M.bias.brn[M.bias.brn$group=="DMB"|M.bias.brn$group=="CFB",]
M.bias.brn.DFCF<-M.bias.brn[M.bias.brn$group=="DFB"|M.bias.brn$group=="CFB",]

p1<-wilcox.test(log(fpkm.brn)~group,data=M.bias.brn.CMDM)[[3]]
p2<-wilcox.test(log(fpkm.brn)~group,data=M.bias.brn.CMDF)[[3]]
p3<-wilcox.test(log(fpkm.brn)~group,data=M.bias.brn.CMCF)[[3]]
p4<-wilcox.test(log(fpkm.brn)~group,data=M.bias.brn.DMDF)[[3]]
p5<-wilcox.test(log(fpkm.brn)~group,data=M.bias.brn.DMCF)[[3]]
p6<-wilcox.test(log(fpkm.brn)~group,data=M.bias.brn.DFCF)[[3]]

pvals<-c(p1,p2,p3,p4,p5,p6)
p.adjust(pvals,method="BH")

#Female brains
F.bias.brn.CMDM<-F.bias.brn[F.bias.brn$group=="CMB"|F.bias.brn$group=="DMB",]
F.bias.brn.CMDF<-F.bias.brn[F.bias.brn$group=="CMB"|F.bias.brn$group=="DFB",]
F.bias.brn.CMCF<-F.bias.brn[F.bias.brn$group=="CMB"|F.bias.brn$group=="CFB",]
F.bias.brn.DMDF<-F.bias.brn[F.bias.brn$group=="DMB"|F.bias.brn$group=="DFB",]
F.bias.brn.DMCF<-F.bias.brn[F.bias.brn$group=="DMB"|F.bias.brn$group=="CFB",]
F.bias.brn.DFCF<-F.bias.brn[F.bias.brn$group=="DFB"|F.bias.brn$group=="CFB",]

p1<-wilcox.test(log(fpkm.brn)~group,data=F.bias.brn.CMDM)[[3]]
p2<-wilcox.test(log(fpkm.brn)~group,data=F.bias.brn.CMDF)[[3]]
p3<-wilcox.test(log(fpkm.brn)~group,data=F.bias.brn.CMCF)[[3]]
p4<-wilcox.test(log(fpkm.brn)~group,data=F.bias.brn.DMDF)[[3]]
p5<-wilcox.test(log(fpkm.brn)~group,data=F.bias.brn.DMCF)[[3]]
p6<-wilcox.test(log(fpkm.brn)~group,data=F.bias.brn.DFCF)[[3]]

pvals<-c(p1,p2,p3,p4,p5,p6)
p.adjust(pvals,method="BH")

#Male head horns
M.bias.hh.CMDM<-M.bias.hh[M.bias.hh$group=="CMH"|M.bias.hh$group=="DMH",]
M.bias.hh.CMDF<-M.bias.hh[M.bias.hh$group=="CMH"|M.bias.hh$group=="DFH",]
M.bias.hh.CMCF<-M.bias.hh[M.bias.hh$group=="CMH"|M.bias.hh$group=="CFH",]
M.bias.hh.DMDF<-M.bias.hh[M.bias.hh$group=="DMH"|M.bias.hh$group=="DFH",]
M.bias.hh.DMCF<-M.bias.hh[M.bias.hh$group=="DMH"|M.bias.hh$group=="CFH",]
M.bias.hh.DFCF<-M.bias.hh[M.bias.hh$group=="DFH"|M.bias.hh$group=="CFH",]

p1<-wilcox.test(log(fpkm.hh)~group,data=M.bias.hh.CMDM)[[3]]
p2<-wilcox.test(log(fpkm.hh)~group,data=M.bias.hh.CMDF)[[3]]
p3<-wilcox.test(log(fpkm.hh)~group,data=M.bias.hh.CMCF)[[3]]
p4<-wilcox.test(log(fpkm.hh)~group,data=M.bias.hh.DMDF)[[3]]
p5<-wilcox.test(log(fpkm.hh)~group,data=M.bias.hh.DMCF)[[3]]
p6<-wilcox.test(log(fpkm.hh)~group,data=M.bias.hh.DFCF)[[3]]

pvals<-c(p1,p2,p3,p4,p5,p6)
p.adjust(pvals,method="BH")

#Female head horns
F.bias.hh.CMDM<-F.bias.hh[F.bias.hh$group=="CMH"|F.bias.hh$group=="DMH",]
F.bias.hh.CMDF<-F.bias.hh[F.bias.hh$group=="CMH"|F.bias.hh$group=="DFH",]
F.bias.hh.CMCF<-F.bias.hh[F.bias.hh$group=="CMH"|F.bias.hh$group=="CFH",]
F.bias.hh.DMDF<-F.bias.hh[F.bias.hh$group=="DMH"|F.bias.hh$group=="DFH",]
F.bias.hh.DMCF<-F.bias.hh[F.bias.hh$group=="DMH"|F.bias.hh$group=="CFH",]
F.bias.hh.DFCF<-F.bias.hh[F.bias.hh$group=="DFH"|F.bias.hh$group=="CFH",]

p1<-wilcox.test(log(fpkm.hh)~group,data=F.bias.hh.CMDM)[[3]]
p2<-wilcox.test(log(fpkm.hh)~group,data=F.bias.hh.CMDF)[[3]]
p3<-wilcox.test(log(fpkm.hh)~group,data=F.bias.hh.CMCF)[[3]]
p4<-wilcox.test(log(fpkm.hh)~group,data=F.bias.hh.DMDF)[[3]]
p5<-wilcox.test(log(fpkm.hh)~group,data=F.bias.hh.DMCF)[[3]]
p6<-wilcox.test(log(fpkm.hh)~group,data=F.bias.hh.DFCF)[[3]]

pvals<-c(p1,p2,p3,p4,p5,p6)
p.adjust(pvals,method="BH")

#Male thoracic horns
M.bias.th.CMDM<-M.bias.th[M.bias.th$group=="CMT"|M.bias.th$group=="DMT",]
M.bias.th.CMDF<-M.bias.th[M.bias.th$group=="CMT"|M.bias.th$group=="DFT",]
M.bias.th.CMCF<-M.bias.th[M.bias.th$group=="CMT"|M.bias.th$group=="CFT",]
M.bias.th.DMDF<-M.bias.th[M.bias.th$group=="DMT"|M.bias.th$group=="DFT",]
M.bias.th.DMCF<-M.bias.th[M.bias.th$group=="DMT"|M.bias.th$group=="CFT",]
M.bias.th.DFCF<-M.bias.th[M.bias.th$group=="DFT"|M.bias.th$group=="CFT",]

p1<-wilcox.test(log(fpkm.th)~group,data=M.bias.th.CMDM)[[3]]
p2<-wilcox.test(log(fpkm.th)~group,data=M.bias.th.CMDF)[[3]]
p3<-wilcox.test(log(fpkm.th)~group,data=M.bias.th.CMCF)[[3]]
p4<-wilcox.test(log(fpkm.th)~group,data=M.bias.th.DMDF)[[3]]
p5<-wilcox.test(log(fpkm.th)~group,data=M.bias.th.DMCF)[[3]]
p6<-wilcox.test(log(fpkm.th)~group,data=M.bias.th.DFCF)[[3]]

pvals<-c(p1,p2,p3,p4,p5,p6)
p.adjust(pvals,method="BH")

#Female thoracic horns
F.bias.th.CMDM<-F.bias.th[F.bias.th$group=="CMT"|F.bias.th$group=="DMT",]
F.bias.th.CMDF<-F.bias.th[F.bias.th$group=="CMT"|F.bias.th$group=="DFT",]
F.bias.th.CMCF<-F.bias.th[F.bias.th$group=="CMT"|F.bias.th$group=="CFT",]
F.bias.th.DMDF<-F.bias.th[F.bias.th$group=="DMT"|F.bias.th$group=="DFT",]
F.bias.th.DMCF<-F.bias.th[F.bias.th$group=="DMT"|F.bias.th$group=="CFT",]
F.bias.th.DFCF<-F.bias.th[F.bias.th$group=="DFT"|F.bias.th$group=="CFT",]

p1<-wilcox.test(log(fpkm.th)~group,data=F.bias.th.CMDM)[[3]]
p2<-wilcox.test(log(fpkm.th)~group,data=F.bias.th.CMDF)[[3]]
p3<-wilcox.test(log(fpkm.th)~group,data=F.bias.th.CMCF)[[3]]
p4<-wilcox.test(log(fpkm.th)~group,data=F.bias.th.DMDF)[[3]]
p5<-wilcox.test(log(fpkm.th)~group,data=F.bias.th.DMCF)[[3]]
p6<-wilcox.test(log(fpkm.th)~group,data=F.bias.th.DFCF)[[3]]

pvals<-c(p1,p2,p3,p4,p5,p6)
p.adjust(pvals,method="BH")

#Male genitalia
M.bias.gen.CMDM<-M.bias.gen[M.bias.gen$group=="CMG"|M.bias.gen$group=="DMG",]
M.bias.gen.CMDF<-M.bias.gen[M.bias.gen$group=="CMG"|M.bias.gen$group=="DFG",]
M.bias.gen.CMCF<-M.bias.gen[M.bias.gen$group=="CMG"|M.bias.gen$group=="CFG",]
M.bias.gen.DMDF<-M.bias.gen[M.bias.gen$group=="DMG"|M.bias.gen$group=="DFG",]
M.bias.gen.DMCF<-M.bias.gen[M.bias.gen$group=="DMG"|M.bias.gen$group=="CFG",]
M.bias.gen.DFCF<-M.bias.gen[M.bias.gen$group=="DFG"|M.bias.gen$group=="CFG",]

p1<-wilcox.test(log(fpkm.gen)~group,data=M.bias.gen.CMDM)[[3]]
p2<-wilcox.test(log(fpkm.gen)~group,data=M.bias.gen.CMDF)[[3]]
p3<-wilcox.test(log(fpkm.gen)~group,data=M.bias.gen.CMCF)[[3]]
p4<-wilcox.test(log(fpkm.gen)~group,data=M.bias.gen.DMDF)[[3]]
p5<-wilcox.test(log(fpkm.gen)~group,data=M.bias.gen.DMCF)[[3]]
p6<-wilcox.test(log(fpkm.gen)~group,data=M.bias.gen.DFCF)[[3]]

pvals<-c(p1,p2,p3,p4,p5,p6)
p.adjust(pvals,method="BH")

#Female genitalia
F.bias.gen.CMDM<-F.bias.gen[F.bias.gen$group=="CMG"|F.bias.gen$group=="DMG",]
F.bias.gen.CMDF<-F.bias.gen[F.bias.gen$group=="CMG"|F.bias.gen$group=="DFG",]
F.bias.gen.CMCF<-F.bias.gen[F.bias.gen$group=="CMG"|F.bias.gen$group=="CFG",]
F.bias.gen.DMDF<-F.bias.gen[F.bias.gen$group=="DMG"|F.bias.gen$group=="DFG",]
F.bias.gen.DMCF<-F.bias.gen[F.bias.gen$group=="DMG"|F.bias.gen$group=="CFG",]
F.bias.gen.DFCF<-F.bias.gen[F.bias.gen$group=="DFG"|F.bias.gen$group=="CFG",]

p1<-wilcox.test(log(fpkm.gen)~group,data=M.bias.gen.CMDM)[[3]]
p2<-wilcox.test(log(fpkm.gen)~group,data=M.bias.gen.CMDF)[[3]]
p3<-wilcox.test(log(fpkm.gen)~group,data=M.bias.gen.CMCF)[[3]]
p4<-wilcox.test(log(fpkm.gen)~group,data=M.bias.gen.DMDF)[[3]]
p5<-wilcox.test(log(fpkm.gen)~group,data=M.bias.gen.DMCF)[[3]]
p6<-wilcox.test(log(fpkm.gen)~group,data=M.bias.gen.DFCF)[[3]]

pvals<-c(p1,p2,p3,p4,p5,p6)
p.adjust(pvals,method="BH")


#Supplementary Figure 3, DSX Targets
#Using with possum data "possum_out_0.005_dsxBindingSites_and_FoldChange.tsv" 

data<-read.table(file.choose(),header=T)

Trimming data so that at least 5 binding sites are present

data<-data[complete.cases(data[,7:16]),]

#For Males
data.dmcm.hh.05<-data.all.dmcm.05[data.all.dmcm.05$Tiss=="hh",]
data.dmcm.th.05<-data.all.dmcm.05[data.all.dmcm.05$Tiss=="th",]
data.dmcm.brn.05<-data.all.dmcm.05[data.all.dmcm.05$Tiss=="brn",]
data.dmcm.gen.05<-data.all.dmcm.05[data.all.dmcm.05$Tiss=="gen",]

data.dmcm.all.05<-rbind(data.dmcm.gen.05,data.dmcm.brn.05,data.dmcm.hh.05,data.dmcm.th.05)
bind.dmcm<-merge(data,data.dmcm.all.05,by.x="GeneID",by.y="Gene_ID")

#Up in dsxRNAi male head horns, that is, inhibited by dsx
up.dsx.male.hh<-as.data.frame(bind.dmcm[bind.dmcm$log2FoldChange>0&bind.dmcm$Tiss=="hh",3])
dim(up.dsx.male.hh)[1]

#Up in control male head horns, that is, activated by dsx
up.cntl.male.hh<-as.data.frame(bind.dmcm[bind.dmcm$log2FoldChange<0&bind.dmcm$Tiss=="hh",3])
dim(up.cntl.male.hh)[1]

#Up in dsxRNAi male thoracic horns, that is, inhibited by dsx
up.dsx.male.th<-as.data.frame(bind.dmcm[bind.dmcm$log2FoldChange>0&bind.dmcm$Tiss=="th",3])
dim(up.dsx.male.th)[1]

#Up in control male thoracic horns, that is, activated by dsx
up.cntl.male.th<-as.data.frame(bind.dmcm[bind.dmcm$log2FoldChange<0&bind.dmcm$Tiss=="th",3])
dim(up.cntl.male.th)[1]

#Up in dsxRNAi male brains, that is, inhibited by dsx
up.dsx.male.brn<-as.data.frame(bind.dmcm[bind.dmcm$log2FoldChange>0&bind.dmcm$Tiss=="brn",3])
dim(up.dsx.male.brn)[1]

#Up in control male brains, that is, activated by dsx
up.cntl.male.brn<-as.data.frame(bind.dmcm[bind.dmcm$log2FoldChange<0&bind.dmcm$Tiss=="brn",3])
dim(up.cntl.male.brn)[1]

#Up in dsxRNAi male genitalia, that is, inhibited by dsx
up.dsx.male.gen<-as.data.frame(bind.dmcm[bind.dmcm$log2FoldChange>0&bind.dmcm$Tiss=="gen",3])
dim(up.dsx.male.gen)[1]

#Up in control male genitalia, that is, activated by dsx
up.cntl.male.gen<-as.data.frame(bind.dmcm[bind.dmcm$log2FoldChange<0&bind.dmcm$Tiss=="gen",3])
dim(up.cntl.male.gen)[1]


#For females
data.dfcf.hh.05<-data.all.dfcf.05[data.all.dfcf.05$Tiss=="hh",]
data.dfcf.th.05<-data.all.dfcf.05[data.all.dfcf.05$Tiss=="th",]
data.dfcf.brn.05<-data.all.dfcf.05[data.all.dfcf.05$Tiss=="brn",]
data.dfcf.gen.05<-data.all.dfcf.05[data.all.dfcf.05$Tiss=="gen",]

data.dfcf.all.05<-rbind(data.dfcf.gen.05,data.dfcf.brn.05,data.dfcf.hh.05,data.dfcf.th.05)
bind.dfcf<-merge(data,data.dfcf.all.05,by.x="GeneID",by.y="Gene.ID")

#Up in dsxRNAi female head horns, that is, inhibited by dsx
up.dsx.female.hh<-as.data.frame(bind.dfcf[bind.dfcf$log2FoldChange>0&bind.dfcf$Tiss=="hh",3])
dim(up.dsx.female.hh)[1]

#Up in control female head horns, that is, activated by dsx
up.cntl.female.hh<-as.data.frame(bind.dfcf[bind.dfcf$log2FoldChange<0&bind.dfcf$Tiss=="hh",3])
dim(up.cntl.female.hh)[1]

#Up in dsxRNAi female thoracic horns, that is, inhibited by dsx
up.dsx.female.th<-as.data.frame(bind.dfcf[bind.dfcf$log2FoldChange>0&bind.dfcf$Tiss=="th",3])
dim(up.dsx.female.th)[1]

#Up in control female thoracic horns, that is, activated by dsx
up.cntl.female.th<-as.data.frame(bind.dfcf[bind.dfcf$log2FoldChange<0&bind.dfcf$Tiss=="th",3])
dim(up.cntl.female.th)[1]

#Up in dsxRNAi female brains, that is, inhibited by dsx
up.dsx.female.brn<-as.data.frame(bind.dfcf[bind.dfcf$log2FoldChange>0&bind.dfcf$Tiss=="brn",3])
dim(up.dsx.female.brn)[1]

#Up in control female brains, that is, activated by dsx
up.cntl.female.brn<-as.data.frame(bind.dfcf[bind.dfcf$log2FoldChange<0&bind.dfcf$Tiss=="brn",3])
dim(up.cntl.female.brn)[1]

#Up in dsxRNAi female genitalia, that is, inhibited by dsx
up.dsx.female.gen<-as.data.frame(bind.dfcf[bind.dfcf$log2FoldChange>0&bind.dfcf$Tiss=="gen",3])
dim(up.dsx.female.brn)[1]

#Up in control female genitalia, that is, activated by dsx
up.cntl.female.gen<-as.data.frame(bind.dfcf[bind.dfcf$log2FoldChange<0&bind.dfcf$Tiss=="gen",3])
dim(up.cntl.female.gen)[1]


#Polarity
#Opposite directions, head horns
merge_HH<-rbind(as.matrix(up.dsx.female.hh),as.matrix(up.cntl.male.hh))
merge_HH<-merge_HH[duplicated(merge_HH),]

merge_HH2<-rbind(as.matrix(up.dsx.male.hh),as.matrix(up.cntl.female.hh))
merge_HH2<-merge_HH2[duplicated(merge_HH2),]

HH_opp<-rbind(as.matrix(merge_HH2),as.matrix(merge_HH))
write.table(HH_opp,file="/Users/crisledon-rettig/Desktop/Indiana University/DSX Project/Binding Sites/Lists/Opp_HH.txt",row.names=FALSE,sep="\t")

#Same directions, head horns
merge_HH3<-rbind(as.matrix(up.dsx.female.hh),as.matrix(up.dsx.male.hh))
merge_HH3<-merge_HH3[duplicated(merge_HH3),]

merge_HH4<-rbind(as.matrix(up.cntl.male.hh),as.matrix(up.cntl.female.hh))
merge_HH4<-merge_HH4[duplicated(merge_HH4),]

HH_same<-rbind(as.matrix(merge_HH3),as.matrix(merge_HH4))
write.table(HH_same,file="/Users/crisledon-rettig/Desktop/Indiana University/DSX Project/Binding Sites/Lists/Same_HH.txt",row.names=FALSE,sep="\t")

#Opposite directions, thoracic horns
merge_TH<-rbind(as.matrix(up.dsx.female.th),as.matrix(up.cntl.male.th))
merge_TH<-merge_TH[duplicated(merge_TH),]

merge_TH2<-rbind(as.matrix(up.dsx.male.th),as.matrix(up.cntl.female.th))
merge_TH2<-merge_TH2[duplicated(merge_TH2),]

TH_opp<-rbind(as.matrix(merge_TH2),as.matrix(merge_TH))
write.table(TH_opp,file="/Users/crisledon-rettig/Desktop/Indiana University/DSX Project/Binding Sites/Lists/Opp_TH.txt",row.names=FALSE,sep="\t")

#Same directions, thoracic horns
merge_TH3<-rbind(as.matrix(up.dsx.female.th),as.matrix(up.dsx.male.th))
merge_TH3<-merge_TH3[duplicated(merge_TH3),]

merge_TH4<-rbind(as.matrix(up.cntl.male.th),as.matrix(up.cntl.female.th))
merge_TH4<-merge_TH4[duplicated(merge_TH4),]

TH_same<-rbind(as.matrix(merge_TH3),as.matrix(merge_TH4))
write.table(TH_same,file="/Users/crisledon-rettig/Desktop/Indiana University/DSX Project/Binding Sites/Lists/Same_TH.txt",row.names=FALSE,sep="\t")

#Do for other tissues

#Overlap
#first for males
up.cntl.male.hh$Tiss<-"HH"
up.dsx.male.hh$Tiss<-"HH"
up.cntl.male.th$Tiss<-"TH"
up.dsx.male.th$Tiss<-"TH"
up.cntl.male.brn$Tiss<-"B"
up.dsx.male.brn$Tiss<-"B"
up.cntl.male.gen$Tiss<-"G"
up.dsx.male.gen$Tiss<-"G"

names(up.cntl.male.hh)<-c("x","Tiss")
names(up.dsx.male.hh)<-c("x","Tiss")
HH.M.dsx<-rbind(up.cntl.male.hh,up.dsx.male.hh)

names(up.cntl.male.th)<-c("x","Tiss")
names(up.dsx.male.th)<-c("x","Tiss")
TH.M.dsx<-rbind(up.cntl.male.th,up.dsx.male.th)

names(up.cntl.male.brn)<-c("x","Tiss")
names(up.dsx.male.brn)<-c("x","Tiss")
BRN.M.dsx<-rbind(up.cntl.male.brn,up.dsx.male.brn)

names(up.cntl.male.gen)<-c("x","Tiss")
names(up.dsx.male.gen)<-c("x","Tiss")
GEN.M.dsx<-rbind(up.cntl.male.gen,up.dsx.male.gen)

all.M<-as.data.frame(rbind(as.matrix(HH.M.dsx),as.matrix(TH.M.dsx),as.matrix(BRN.M.dsx),as.matrix(GEN.M.dsx)))

#total male targets
dim(all.M[!duplicated(all.M[,"x"]),])[1]

#adding counter and tissues
all.M$numdup<-as.numeric("1")
x<-aggregate(numdup~x,data=all.M,FUN=sum)
x.1<-x[x$numdup=="1",]
x.2<-x[x$numdup=="2",]
x.3<-x[x$numdup=="3",]

all.M.HH<-all.M[all.M$Tiss=="HH",]
all.M.TH<-all.M[all.M$Tiss=="TH",]
all.M.G<-all.M[all.M$Tiss=="G",]
all.M.B<-all.M[all.M$Tiss=="B",]

#each tissue
M.HH.1<-merge(x.1,all.M.HH,by="x")
dim(M.HH.1)[1]
M.TH.1<-merge(x.1,all.M.TH,by="x")
dim(M.TH.1)[1]
M.G.1<-merge(x.1,all.M.G,by="x")
dim(M.G.1)[1]
M.B.1<-merge(x.1,all.M.B,by="x")
dim(M.B.1)[1]

#2 tissues
all.M.2<-merge(x.2,all.M,by="x")
all.M.2.HH_TH<-all.M.2[all.M.2$Tiss!="G"&all.M.2$Tiss!="B",]
all.M.2.HH_TH<-all.M.2.HH_TH[duplicated(all.M.2.HH_TH$x),]
dim(all.M.2.HH_TH)[1]

all.M.2.HH_B<-all.M.2[all.M.2$Tiss!="G"&all.M.2$Tiss!="TH",]
all.M.2.HH_B<-all.M.2.HH_B[duplicated(all.M.2.HH_B$x),]
dim(all.M.2.HH_B)[1]

all.M.2.HH_G<-all.M.2[all.M.2$Tiss!="TH"&all.M.2$Tiss!="B",]
all.M.2.HH_G<-all.M.2.HH_G[duplicated(all.M.2.HH_G$x),]
dim(all.M.2.HH_G)[1]

all.M.2.G_TH<-all.M.2[all.M.2$Tiss!="HH"&all.M.2$Tiss!="B",]
all.M.2.G_TH<-all.M.2.G_TH[duplicated(all.M.2.G_TH$x),]
dim(all.M.2.G_TH)[1]

all.M.2.B_TH<-all.M.2[all.M.2$Tiss!="HH"&all.M.2$Tiss!="G",]
all.M.2.B_TH<-all.M.2.B_TH[duplicated(all.M.2.B_TH$x),]
dim(all.M.2.B_TH)[1]

all.M.2.G_B<-all.M.2[all.M.2$Tiss!="HH"&all.M.2$Tiss!="TH",]
all.M.2.G_B<-all.M.2.G_B[duplicated(all.M.2.G_B$x),]
dim(all.M.2.G_B)[1]

#3 tissues
all.M.3<-merge(x.3,all.M,by="x")
all.M.3.HH_TH_G<-all.M.3[all.M.3$Tiss!="B",]
all.M.3.HH_TH_G<-all.M.3.HH_TH_G[unique(all.M.3.HH_TH_G$x),]
dim(all.M.3.HH_TH_G)[1]


#now for females
up.cntl.female.hh$Tiss<-"HH"
up.dsx.female.hh$Tiss<-"HH"
up.cntl.female.th$Tiss<-"TH"
up.dsx.female.th$Tiss<-"TH"
up.cntl.female.brn$Tiss<-"B"
up.dsx.female.brn$Tiss<-"B"
up.cntl.female.gen$Tiss<-"G"
up.dsx.female.gen$Tiss<-"G"

names(up.cntl.female.hh)<-c("x","Tiss")
names(up.dsx.female.hh)<-c("x","Tiss")
HH.F.dsx<-rbind(up.cntl.female.hh,up.dsx.female.hh)

names(up.cntl.female.th)<-c("x","Tiss")
names(up.dsx.female.th)<-c("x","Tiss")
TH.F.dsx<-rbind(up.cntl.female.th,up.dsx.female.th)

names(up.cntl.female.brn)<-c("x","Tiss")
names(up.dsx.female.brn)<-c("x","Tiss")
BRN.F.dsx<-rbind(up.cntl.female.brn,up.dsx.female.brn)

names(up.cntl.female.gen)<-c("x","Tiss")
names(up.dsx.female.gen)<-c("x","Tiss")
GEN.F.dsx<-rbind(up.cntl.female.gen,up.dsx.female.gen)

all.F<-as.data.frame(rbind(as.matrix(HH.F.dsx),as.matrix(TH.F.dsx),as.matrix(BRN.F.dsx),as.matrix(GEN.F.dsx)))

#total female targets
dim(all.F[!duplicated(all.F[,"x"]),])[1]

#adding counter and tissues
all.F$numdup<-as.numeric("1")
x<-aggregate(numdup~x,data=all.F,FUN=sum)
x.1<-x[x$numdup=="1",]
x.2<-x[x$numdup=="2",]
x.3<-x[x$numdup=="3",]

all.F.HH<-all.F[all.F$Tiss=="HH",]
all.F.TH<-all.F[all.F$Tiss=="TH",]
all.F.G<-all.F[all.F$Tiss=="G",]
all.F.B<-all.F[all.F$Tiss=="B",]

#each tissue
F.HH.1<-merge(x.1,all.F.HH,by="x")
dim(F.HH.1)[1]
F.TH.1<-merge(x.1,all.F.TH,by="x")
dim(F.TH.1)[1]
F.G.1<-merge(x.1,all.F.G,by="x")
dim(F.G.1)[1]
F.B.1<-merge(x.1,all.F.B,by="x")
dim(F.B.1)[1]

#2 tissues
all.F.2<-merge(x.2,all.F,by="x")
all.F.2.HH_TH<-all.F.2[all.F.2$Tiss!="G"&all.F.2$Tiss!="B",]
all.F.2.HH_TH<-all.F.2.HH_TH[duplicated(all.F.2.HH_TH$x),]
dim(all.F.2.HH_TH)[1]

all.F.2.HH_B<-all.F.2[all.F.2$Tiss!="G"&all.F.2$Tiss!="TH",]
all.F.2.HH_B<-all.F.2.HH_B[duplicated(all.F.2.HH_B$x),]
dim(all.F.2.HH_B)[1]

all.F.2.HH_G<-all.F.2[all.F.2$Tiss!="TH"&all.F.2$Tiss!="B",]
all.F.2.HH_G<-all.F.2.HH_G[duplicated(all.F.2.HH_G$x),]
dim(all.F.2.HH_G)[1]

all.F.2.G_TH<-all.F.2[all.F.2$Tiss!="HH"&all.F.2$Tiss!="B",]
all.F.2.G_TH<-all.F.2.G_TH[duplicated(all.F.2.G_TH$x),]
dim(all.F.2.G_TH)[1]

all.F.2.B_TH<-all.F.2[all.F.2$Tiss!="HH"&all.F.2$Tiss!="G",]
all.F.2.B_TH<-all.F.2.B_TH[duplicated(all.F.2.B_TH$x),]
dim(all.F.2.B_TH)[1]

all.F.2.G_B<-all.F.2[all.F.2$Tiss!="HH"&all.F.2$Tiss!="TH",]
all.F.2.G_B<-all.F.2.G_B[duplicated(all.F.2.G_B$x),]
dim(all.F.2.G_B)[1]

#3 tissues
all.F.3<-merge(x.3,all.F,by="x")
all.F.3.HH_TH_G<-all.F.3[all.F.3$Tiss!="B",]
all.F.3.HH_TH_G<-all.F.3.HH_TH_G[unique(all.F.3.HH_TH_G$x),]
dim(all.F.3.HH_TH_G)[1]
