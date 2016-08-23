data.gen<-read.table(file.choose(),header=T,sep="\t")
data.brn<-read.table(file.choose(),header=T,sep="\t")
data.hh<-read.table(file.choose(),header=T,sep="\t")
data.th<-read.table(file.choose(),header=T,sep="\t")

data.gen$Tiss<-"gen"
data.brn$Tiss<-"brn"
data.hh$Tiss<-"hh"
data.th$Tiss<-"th"

data.all<-rbind(data.gen,data.brn,data.hh,data.th)
data.all$Tiss<-as.factor(data.all$Tiss)

data.all$SexBias[data.all$logFC>=0]<-"F"
data.all$SexBias[data.all$logFC<=0]<-"M"
data.all$SexBias<-as.factor(data.all$SexBias)

data.all.05<-data.all[data.all$FDR<=0.05,]
data.all.1<-data.all[data.all$FDR<=0.1,]
data.all.2<-data.all[data.all$FDR<=0.2,]
data.all.3<-data.all[data.all$FDR<=0.3,]

#Barcharts, ordering by brain, genetalia, th then hh

data.gen.3<-data.all.3[data.all.3$Tiss=="gen",]
data.brn.3<-data.all.3[data.all.3$Tiss=="brn",]
data.th.3<-data.all.3[data.all.3$Tiss=="th",]
data.hh.3<-data.all.3[data.all.3$Tiss=="hh",]

data.gen.3F<-data.gen.3[data.gen.3$SexBias=="F",]
data.gen.3M<-data.gen.3[data.gen.3$SexBias=="M",]

data.brn.3F<-data.brn.3[data.brn.3$SexBias=="F",]
data.brn.3M<-data.brn.3[data.brn.3$SexBias=="M",]

data.th.3F<-data.th.3[data.th.3$SexBias=="F",]
data.th.3M<-data.th.3[data.th.3$SexBias=="M",]

data.hh.3F<-data.hh.3[data.hh.3$SexBias=="F",]
data.hh.3M<-data.hh.3[data.hh.3$SexBias=="M",]

gen.bar.3<-c(dim(data.gen.3F)[1],dim(data.gen.3M)[1])
brn.bar.3<-c(dim(data.brn.3F)[1],dim(data.brn.3M)[1])
th.bar.3<-c(dim(data.th.3F)[1],dim(data.th.3M)[1])
hh.bar.3<-c(dim(data.hh.3F)[1],dim(data.hh.3M)[1])

data.3.bar<-rbind(brn.bar.3,gen.bar.3,th.bar.3,hh.bar.3)
data.3.bar<-as.data.frame(t(data.3.bar))
colnames(data.3.bar)<-c("brain","genitalia","thor_horn","head_horn")
rownames(data.3.bar)<-c("female","male")

par(oma=c(0,1,0,0))
barplot(as.matrix(data.3.bar),beside=TRUE,main="",ylim=c(0,550),cex.main=1.3,legend=rownames(data.3.bar))
box()
mtext("Biased Reads",2,cex=1.3,line=3,font=2)
mtext("FDR <= 0.3",3,cex=1.3,line=1,font=2)

data.gen.1<-data.all.1[data.all.1$Tiss=="gen",]
data.brn.1<-data.all.1[data.all.1$Tiss=="brn",]
data.th.1<-data.all.1[data.all.1$Tiss=="th",]
data.hh.1<-data.all.1[data.all.1$Tiss=="hh",]

data.gen.1F<-data.gen.1[data.gen.1$SexBias=="F",]
data.gen.1M<-data.gen.1[data.gen.1$SexBias=="M",]

data.brn.1F<-data.brn.1[data.brn.1$SexBias=="F",]
data.brn.1M<-data.brn.1[data.brn.1$SexBias=="M",]

data.th.1F<-data.th.1[data.th.1$SexBias=="F",]
data.th.1M<-data.th.1[data.th.1$SexBias=="M",]

data.hh.1F<-data.hh.1[data.hh.1$SexBias=="F",]
data.hh.1M<-data.hh.1[data.hh.1$SexBias=="M",]

gen.bar.1<-c(dim(data.gen.1F)[1],dim(data.gen.1M)[1])
brn.bar.1<-c(dim(data.brn.1F)[1],dim(data.brn.1M)[1])
th.bar.1<-c(dim(data.th.1F)[1],dim(data.th.1M)[1])
hh.bar.1<-c(dim(data.hh.1F)[1],dim(data.hh.1M)[1])

data.1.bar<-rbind(brn.bar.1,gen.bar.1,th.bar.1,hh.bar.1)
data.1.bar<-as.data.frame(t(data.1.bar))
colnames(data.1.bar)<-c("brain","genitalia","thor_horn","head_horn")
rownames(data.1.bar)<-c("female","male")

par(oma=c(0,1,0,0))
barplot(as.matrix(data.1.bar),beside=TRUE,main="",ylim=c(0,550),cex.main=1.3,legend=rownames(data.1.bar))
box()
mtext("Biased Reads",2,cex=1.3,line=3,font=2)
mtext("FDR <= 0.1",3,cex=1.3,line=1,font=2)

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
barplot(as.matrix(data.05.bar),beside=TRUE,main="",ylim=c(0,550),cex.main=1.3,legend=rownames(data.05.bar))
box()
mtext("Biased Reads",2,cex=1.3,line=3,font=2)
mtext("FDR <= 0.05",3,cex=1.3,line=1,font=2)


data.all.3.uniq<-data.all.3[!duplicated(data.all.3$GeneID),]

data.all.3.uniq.brn<-data.all.3.uniq[data.all.3.uniq$Tiss=="brn",]
data.all.3.uniq.gen<-data.all.3.uniq[data.all.3.uniq$Tiss=="gen",]
data.all.3.uniq.th<-data.all.3.uniq[data.all.3.uniq$Tiss=="th",]
data.all.3.uniq.hh<-data.all.3.uniq[data.all.3.uniq$Tiss=="hh",]

data.all.3.uniq.brn.F<-data.all.3.uniq.brn[data.all.3.uniq.brn$Sex=="F",]
data.all.3.uniq.brn.M<-data.all.3.uniq.brn[data.all.3.uniq.brn$Sex=="M",]

data.all.3.uniq.gen.F<-data.all.3.uniq.gen[data.all.3.uniq.gen$Sex=="F",]
data.all.3.uniq.gen.M<-data.all.3.uniq.gen[data.all.3.uniq.gen$Sex=="M",]

data.all.3.uniq.th.F<-data.all.3.uniq.th[data.all.3.uniq.th$Sex=="F",]
data.all.3.uniq.th.M<-data.all.3.uniq.th[data.all.3.uniq.th$Sex=="M",]

data.all.3.uniq.hh.F<-data.all.3.uniq.hh[data.all.3.uniq.hh$Sex=="F",]
data.all.3.uniq.hh.M<-data.all.3.uniq.hh[data.all.3.uniq.hh$Sex=="M",]

brn.bar.3.uniq<-c(dim(data.all.3.uniq.brn.F)[1],dim(data.all.3.uniq.brn.M)[1])
gen.bar.3.uniq<-c(dim(data.all.3.uniq.gen.F)[1],dim(data.all.3.uniq.gen.M)[1])
th.bar.3.uniq<-c(dim(data.all.3.uniq.th.F)[1],dim(data.all.3.uniq.th.M)[1])
hh.bar.3.uniq<-c(dim(data.all.3.uniq.hh.F)[1],dim(data.all.3.uniq.hh.M)[1])

data.3.bar.uniq<-rbind(brn.bar.3.uniq,gen.bar.3.uniq,th.bar.3.uniq,hh.bar.3.uniq)
data.3.bar.uniq<-as.data.frame(t(data.3.bar.uniq))
colnames(data.3.bar.uniq)<-c("brain","genitalia","thor_horn","head_horn")
rownames(data.3.bar.uniq)<-c("female","male")

par(oma=c(0,1,0,0))
barplot(as.matrix(data.3.bar.uniq),beside=TRUE,main="",ylim=c(0,550),cex.main=1.3,legend=rownames(data.3.bar.uniq))
box()
mtext("Biased Reads, Unique",2,cex=1.3,line=3,font=2)
mtext("FDR <= 0.3",3,cex=1.3,line=1,font=2)

data.all.3.F<-data.all.3[data.all.3$Sex=="F",]
data.all.3.M<-data.all.3[data.all.3$Sex=="M",]

data.all.3.F.dup<-data.all.3.F[duplicated(data.all.3.F$GeneID)|duplicated(data.all.3.F$GeneID,fromLast = TRUE),]

data.all.3.M.dup<-data.all.3.M[duplicated(data.all.3.M$GeneID)|duplicated(data.all.3.M$GeneID,fromLast = TRUE),]

data.all.3.F.dup$numdup<-as.numeric("1")
data.all.3.M.dup$numdup<-as.numeric("1")

dup.tab.F<-aggregate(numdup~GeneID,data=data.all.3.F.dup,FUN=sum)
dup.tab.M<-aggregate(numdup~GeneID,data=data.all.3.M.dup,FUN=sum)

data.all.3.F<-data.all.3[data.all.3$Sex=="F",]
data.all.3.M<-data.all.3[data.all.3$Sex=="M",]

data.all.3.F$numdup<-as.numeric("1")
data.all.3.M$numdup<-as.numeric("1")

spec.tab.F<-aggregate(numdup~GeneID,data=data.all.3.F,FUN=sum)
spec.tab.M<-aggregate(numdup~GeneID,data=data.all.3.M,FUN=sum)

Tiss.1.F<-dim(spec.tab.F[spec.tab.F$numdup=="1",])[1]
Tiss.1.M<-dim(spec.tab.M[spec.tab.M$numdup=="1",])[1]

Tiss.2.F<-dim(spec.tab.F[spec.tab.F$numdup=="2",])[1]
Tiss.2.M<-dim(spec.tab.M[spec.tab.M$numdup=="2",])[1]

Tiss.3.F<-dim(spec.tab.F[spec.tab.F$numdup=="3",])[1]
Tiss.3.M<-dim(spec.tab.M[spec.tab.M$numdup=="3",])[1]

Tiss.4.F<-dim(spec.tab.F[spec.tab.F$numdup=="4",])[1]
Tiss.4.M<-dim(spec.tab.M[spec.tab.M$numdup=="4",])[1]

spec.3.bar.1<-c(Tiss.1.F,Tiss.1.M)
spec.3.bar.2<-c(Tiss.2.F,Tiss.2.M)
spec.3.bar.3<-c(Tiss.3.F,Tiss.3.M)
spec.3.bar.4<-c(Tiss.4.F,Tiss.4.M)

spec.3.bar<-rbind(spec.3.bar.1,spec.3.bar.2,spec.3.bar.3,spec.3.bar.4)
spec.3.bar<-as.data.frame(t(spec.3.bar))
colnames(spec.3.bar)<-c("1 Tiss","2 Tiss","3 Tiss","4 Tiss")
rownames(spec.3.bar)<-c("female","male")

par(oma=c(0,1,0,0))
barplot(as.matrix(spec.3.bar),beside=TRUE,main="",ylim=c(0,650),cex.main=1.3,legend=rownames(spec.3.bar))
box()
mtext("Biased Reads",2,cex=1.3,line=3,font=2)
mtext("FDR <= 0.3",3,cex=1.3,line=1,font=2)


Tiss.1.F.frame<-spec.tab.F[spec.tab.F$numdup=="1",]
Tiss.1.M.frame<-spec.tab.M[spec.tab.M$numdup=="1",]

Tiss.2.F.frame<-spec.tab.F[spec.tab.F$numdup=="2",]
Tiss.2.M.frame<-spec.tab.M[spec.tab.M$numdup=="2",]

Tiss.3.F.frame<-spec.tab.F[spec.tab.F$numdup=="3",]
Tiss.3.M.frame<-spec.tab.M[spec.tab.M$numdup=="3",]

Tiss.4.F.frame<-spec.tab.F[spec.tab.F$numdup=="4",]
Tiss.4.M.frame<-spec.tab.M[spec.tab.M$numdup=="4",]

#Tables
Tiss.1.F.test<-merge(Tiss.1.F.frame,data.all.3.F,by.x="GeneID",by.y="GeneID")
aggregate(numdup.y~Tiss,data=Tiss.1.F.test,FUN=sum)

Tiss.1.M.test<-merge(Tiss.1.M.frame,data.all.3.M,by.x="GeneID",by.y="GeneID")
aggregate(numdup.y~Tiss,data=Tiss.1.M.test,FUN=sum)

Tiss.2.F.test<-merge(Tiss.2.F.frame,data.all.3.F,by.x="GeneID",by.y="GeneID")

Tiss.2.M.test<-merge(Tiss.2.M.frame,data.all.3.M,by.x="GeneID",by.y="GeneID")

Tiss.3.F.test<-merge(Tiss.3.F.frame,data.all.3.F,by.x="GeneID",by.y="GeneID")

Tiss.3.M.test<-merge(Tiss.3.M.frame,data.all.3.M,by.x="GeneID",by.y="GeneID")

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

#Head
Tiss.1.F.H<-Tiss.1.F[Tiss.1.F$Tiss=="hh",]


#Head/Thoracic
Tiss.2.F.HT<-dim(merge(Tiss.2.F.H,Tiss.2.F.T,by="GeneID"))[1]
Tiss.2.M.HT<-dim(merge(Tiss.2.M.H,Tiss.2.M.T,by="GeneID"))[1]

#Head/Genitalia
Tiss.2.F.HG<-dim(merge(Tiss.2.F.H,Tiss.2.F.G,by="GeneID"))[1]
Tiss.2.M.HG<-dim(merge(Tiss.2.M.H,Tiss.2.M.G,by="GeneID"))[1]

#Thoracic/Genitalia
Tiss.2.F.TG<-dim(merge(Tiss.2.F.T,Tiss.2.F.G,by="GeneID"))[1]
Tiss.2.M.TG<-dim(merge(Tiss.2.M.T,Tiss.2.M.G,by="GeneID"))[1]

#Head/Brain
Tiss.2.F.HB<-dim(merge(Tiss.2.F.H,Tiss.2.F.B,by="GeneID"))[1]
Tiss.2.M.HB<-dim(merge(Tiss.2.M.H,Tiss.2.M.B,by="GeneID"))[1]

#Thoracic/Brain
Tiss.2.F.TB<-dim(merge(Tiss.2.F.T,Tiss.2.F.B,by="GeneID"))[1]
Tiss.2.M.TB<-dim(merge(Tiss.2.M.T,Tiss.2.M.B,by="GeneID"))[1]

#Genitalia/Brain
Tiss.2.F.GB<-dim(merge(Tiss.2.F.G,Tiss.2.F.B,by="GeneID"))[1]
Tiss.2.M.GB<-dim(merge(Tiss.2.M.G,Tiss.2.M.B,by="GeneID"))[1]

#Head/Thoracic/Genitalia
Tiss.3.F.HT<-merge(Tiss.3.F.H,Tiss.3.F.T,by="GeneID")
Tiss.3.F.HTG<-dim(merge(Tiss.3.F.HT,Tiss.3.F.G,by="GeneID"))[1]

Tiss.3.M.HT<-merge(Tiss.3.M.H,Tiss.3.M.T,by="GeneID")
Tiss.3.M.HTG<-dim(merge(Tiss.3.M.HT,Tiss.3.M.G,by="GeneID"))[1]

#Head/Thoracic/Brain
Tiss.3.F.HTB<-dim(merge(Tiss.3.F.HT,Tiss.3.F.B,by="GeneID"))[1]

Tiss.3.M.HTB<-dim(merge(Tiss.3.M.HT,Tiss.3.M.B,by="GeneID"))[1]

#Head/Genitalia/Brain
Tiss.3.F.HG<-merge(Tiss.3.F.H,Tiss.3.F.G,by="GeneID")
Tiss.3.F.HGB<-dim(merge(Tiss.3.F.HG,Tiss.3.F.B,by="GeneID"))[1]

Tiss.3.M.HG<-merge(Tiss.3.M.H,Tiss.3.M.G,by="GeneID")
Tiss.3.M.HGB<-dim(merge(Tiss.3.M.HG,Tiss.3.M.B,by="GeneID"))[1]

#Brain/Thoracic/Genitalia
Tiss.3.F.BT<-merge(Tiss.3.F.B,Tiss.3.F.T,by="GeneID")
Tiss.3.F.BTG<-dim(merge(Tiss.3.F.BT,Tiss.3.F.G,by="GeneID"))[1]

Tiss.3.M.BT<-merge(Tiss.3.M.B,Tiss.3.M.T,by="GeneID")
Tiss.3.M.BTG<-dim(merge(Tiss.3.M.BT,Tiss.3.M.G,by="GeneID"))[1]


data.gen.d<-read.table(file.choose(),header=T,sep="\t")
data.brn.d<-read.table(file.choose(),header=T,sep="\t")
data.hh.d<-read.table(file.choose(),header=T,sep="\t")
data.th.d<-read.table(file.choose(),header=T,sep="\t")

data.gen.d$Tiss<-"gen"
data.brn.d$Tiss<-"brn"
data.hh.d$Tiss<-"hh"
data.th.d$Tiss<-"th"

data.all.d<-rbind(data.gen.d,data.brn.d,data.hh.d,data.th.d)
data.all.d$Tiss<-as.factor(data.all.d$Tiss)

data.all.d$Bias[data.all.d$logFC>=0]<-"F"
data.all.d$Bias[data.all.d$logFC<=0]<-"M"
data.all.d$Bias<-as.factor(data.all.d$Bias)

data.all.d.05<-data.all.d[data.all.d$FDR<=0.05,]
data.all.d.1<-data.all.d[data.all.d$FDR<=0.1,]
data.all.d.2<-data.all.d[data.all.d$FDR<=0.2,]
data.all.d.3<-data.all.d[data.all.d$FDR<=0.3,]



data.gen.d.3<-data.all.d.3[data.all.d.3$Tiss=="gen",]
data.brn.d.3<-data.all.d.3[data.all.d.3$Tiss=="brn",]
data.th.d.3<-data.all.d.3[data.all.d.3$Tiss=="th",]
data.hh.d.3<-data.all.d.3[data.all.d.3$Tiss=="hh",]

data.gen.d.3F<-data.gen.d.3[data.gen.d.3$Bias=="F",]
data.gen.d.3M<-data.gen.d.3[data.gen.d.3$Bias=="M",]

data.brn.d.3F<-data.brn.d.3[data.brn.d.3$Bias=="F",]
data.brn.d.3M<-data.brn.d.3[data.brn.d.3$Bias=="M",]

data.th.d.3F<-data.th.d.3[data.th.d.3$Bias=="F",]
data.th.d.3M<-data.th.d.3[data.th.d.3$Bias=="M",]

data.hh.d.3F<-data.hh.d.3[data.hh.d.3$Bias=="F",]
data.hh.d.3M<-data.hh.d.3[data.hh.d.3$Bias=="M",]

gen.bar.d.3<-c(dim(data.gen.d.3F)[1],dim(data.gen.d.3M)[1])
brn.bar.d.3<-c(dim(data.brn.d.3F)[1],dim(data.brn.d.3M)[1])
th.bar.d.3<-c(dim(data.th.d.3F)[1],dim(data.th.d.3M)[1])
hh.bar.d.3<-c(dim(data.hh.d.3F)[1],dim(data.hh.d.3M)[1])

data.d.3.bar<-rbind(brn.bar.d.3,gen.bar.d.3,th.bar.d.3,hh.bar.d.3)
data.d.3.bar<-as.data.frame(t(data.d.3.bar))
colnames(data.d.3.bar)<-c("brain","genitalia","thor_horn","head_horn")
rownames(data.d.3.bar)<-c("female","male")

par(oma=c(0,1,0,0))
barplot(as.matrix(data.d.3.bar),beside=TRUE,main="",ylim=c(0,550),cex.main=1.3,legend=rownames(data.d.3.bar))
box()
mtext("Biased Reads",2,cex=1.3,line=3,font=2)
mtext("FDR <= 0.3",3,cex=1.3,line=1,font=2)

data.all.d.3.F<-data.all.d.3[data.all.d.3$Bias=="F",]
data.all.d.3.M<-data.all.d.3[data.all.d.3$Bias=="M",]

data.all.d.3.F$numdup<-as.numeric("1")
data.all.d.3.M$numdup<-as.numeric("1")

spec.tab.F<-aggregate(numdup~Gene.ID,data=data.all.d.3.F,FUN=sum)
spec.tab.M<-aggregate(numdup~Gene.ID,data=data.all.d.3.M,FUN=sum)

Tiss.1.F<-dim(spec.tab.F[spec.tab.F$numdup=="1",])[1]
Tiss.1.M<-dim(spec.tab.M[spec.tab.M$numdup=="1",])[1]

Tiss.2.F<-dim(spec.tab.F[spec.tab.F$numdup=="2",])[1]
Tiss.2.M<-dim(spec.tab.M[spec.tab.M$numdup=="2",])[1]

Tiss.3.F<-dim(spec.tab.F[spec.tab.F$numdup=="3",])[1]
Tiss.3.M<-dim(spec.tab.M[spec.tab.M$numdup=="3",])[1]

Tiss.4.F<-dim(spec.tab.F[spec.tab.F$numdup=="4",])[1]
Tiss.4.M<-dim(spec.tab.M[spec.tab.M$numdup=="4",])[1]

spec.3.bar.1<-c(Tiss.1.F,Tiss.1.M)
spec.3.bar.2<-c(Tiss.2.F,Tiss.2.M)
spec.3.bar.3<-c(Tiss.3.F,Tiss.3.M)
spec.3.bar.4<-c(Tiss.4.F,Tiss.4.M)

spec.3.bar<-rbind(spec.3.bar.1,spec.3.bar.2,spec.3.bar.3,spec.3.bar.4)
spec.3.bar<-as.data.frame(t(spec.3.bar))
colnames(spec.3.bar)<-c("1 Tiss","2 Tiss","3 Tiss","4 Tiss")
rownames(spec.3.bar)<-c("female","male")

par(oma=c(0,1,0,0))
barplot(as.matrix(spec.3.bar),beside=TRUE,main="",ylim=c(0,750),cex.main=1.3,legend=rownames(spec.3.bar))
box()
mtext("Biased Reads",2,cex=1.3,line=3,font=2)
mtext("FDR <= 0.3",3,cex=1.3,line=1,font=2)

data.gen<-read.table(file.choose(),header=T,sep="\t")
data.brn<-read.table(file.choose(),header=T,sep="\t")
data.hh<-read.table(file.choose(),header=T,sep="\t")
data.th<-read.table(file.choose(),header=T,sep="\t")

data.gen$Tiss<-"gen"
data.brn$Tiss<-"brn"
data.hh$Tiss<-"hh"
data.th$Tiss<-"th"

data.all<-rbind(data.gen,data.brn,data.hh,data.th)
data.all$Tiss<-as.factor(data.all$Tiss)

data.all$SexBias[data.all$logFC>=0]<-"F"
data.all$SexBias[data.all$logFC<=0]<-"M"
data.all$SexBias<-as.factor(data.all$SexBias)
data.all.3<-data.all[data.all$FDR<=0.3,]

data.all.3.F<-data.all.3[data.all.3$Sex=="F",]
data.all.3.M<-data.all.3[data.all.3$Sex=="M",]

data.all.d.3.F<-data.all.d.3[data.all.d.3$Bias=="F",]
data.all.d.3.M<-data.all.d.3[data.all.d.3$Bias=="M",]

##Data from data.all.d.3 that was in data.all.3

data.diff.3.F<-merge(data.all.3.F,data.all.d.3.F,by.x="GeneID",by.y="Gene.ID",all=FALSE)
data.diff.3.M<-merge(data.all.3.M,data.all.d.3.M,by.x="GeneID",by.y="Gene.ID",all=FALSE)

F.bar.3<-c(dim(data.all.3.F)[1],dim(data.all.d.3.F)[1],dim(data.diff.3.F)[1])
M.bar.3<-c(dim(data.all.3.M)[1],dim(data.all.d.3.M)[1],dim(data.diff.3.M)[1])
bar.3<-rbind(F.bar.3,M.bar.3)
bar.3<-as.data.frame(t(bar.3))

colnames(bar.3)<-c("Female Biased","Male Biased")
rownames(bar.3)<-c("Control","DSX RNAi","Common")

par(oma=c(0,1,0,0))
barplot(as.matrix(bar.3),beside=TRUE,main="",ylim=c(0,1000),cex.main=1.3,legend=rownames(bar.3))
box()
mtext("Biased Reads",2,cex=1.3,line=3,font=2)
mtext("FDR <= 0.3",3,cex=1.3,line=1,font=2)


##Data from data.all.d.3 that was in data.all.3

data.diff.3.F<-merge(data.all.3.F,data.all.d.3.F,by.x=c("GeneID","Tiss"),by.y=c("Gene.ID","Tiss"),all=FALSE)
data.diff.3.M<-merge(data.all.3.M,data.all.d.3.M,by.x=c("GeneID","Tiss"),by.y=c("Gene.ID","Tiss"),all=FALSE)

F.bar.3<-c(dim(data.all.3.F)[1],dim(data.all.d.3.F)[1],dim(data.diff.3.F)[1])
M.bar.3<-c(dim(data.all.3.M)[1],dim(data.all.d.3.M)[1],dim(data.diff.3.M)[1])
bar.3<-rbind(F.bar.3,M.bar.3)
bar.3<-as.data.frame(t(bar.3))

colnames(bar.3)<-c("Female Biased","Male Biased")
rownames(bar.3)<-c("Control","DSX RNAi","Common")

par(oma=c(0,1,0,0))
barplot(as.matrix(bar.3),beside=TRUE,main="",ylim=c(0,1000),cex.main=1.3,legend=rownames(bar.3))
box()
mtext("Biased Reads",2,cex=1.3,line=3,font=2)
mtext("FDR <= 0.3",3,cex=1.3,line=1,font=2)


data.diff.gen.3F<-data.diff.3.F[data.diff.3.F$Tiss=="gen",]
data.diff.gen.3M<-data.diff.3.M[data.diff.3.M$Tiss=="gen",]

data.diff.brn.3F<-data.diff.3.F[data.diff.3.F$Tiss=="brn",]
data.diff.brn.3M<-data.diff.3.M[data.diff.3.M$Tiss=="brn",]

data.diff.th.3F<-data.diff.3.F[data.diff.3.F$Tiss=="th",]
data.diff.th.3M<-data.diff.3.M[data.diff.3.M$Tiss=="th",]

data.diff.hh.3F<-data.diff.3.F[data.diff.3.F$Tiss=="hh",]
data.diff.hh.3M<-data.diff.3.M[data.diff.3.M$Tiss=="hh",]

gen.bar.diff.3<-c(dim(data.diff.gen.3F)[1],dim(data.diff.gen.3M)[1])
brn.bar.diff.3<-c(dim(data.diff.brn.3F)[1],dim(data.diff.brn.3M)[1])
th.bar.diff.3<-c(dim(data.diff.th.3F)[1],dim(data.diff.th.3M)[1])
hh.bar.diff.3<-c(dim(data.diff.hh.3F)[1],dim(data.diff.hh.3M)[1])

data.diff.3.bar<-rbind(brn.bar.diff.3,gen.bar.diff.3,th.bar.diff.3,hh.bar.diff.3)
data.diff.3.bar<-as.data.frame(t(data.diff.3.bar))
colnames(data.diff.3.bar)<-c("brain","genitalia","thor_horn","head_horn")
rownames(data.diff.3.bar)<-c("female","male")

par(oma=c(0,1,0,0))
barplot(as.matrix(data.diff.3.bar),beside=TRUE,main="",ylim=c(0,300),cex.main=1.3,legend=rownames(data.diff.3.bar))
box()
mtext("Biased Reads",2,cex=1.3,line=3,font=2)
mtext("FDR <= 0.3",3,cex=1.3,line=1,font=2)




