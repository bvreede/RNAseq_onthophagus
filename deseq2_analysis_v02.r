#read input table
counts_all <- read.table("/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/counts_all.tsv",header=T,row.names=1)

#generate phenotype data from column names
all_names <- as.vector(colnames(counts_all))
temp_li <- strsplit(all_names, "_")
temp_mat <- matrix(unlist(temp_li), ncol=5,byrow=TRUE)
category_df <- as.data.frame(temp_mat)
colnames(category_df) <- c("size","treatment","sex","tissue","sample")
rownames(category_df) <- all_names

#add group category to phenotype data
pasted <- c('size','treatment','sex','tissue')
category_df$group <- apply(category_df[,pasted],1,paste,collapse='')
category_df$group <- as.factor(category_df$group)

#run deseq analysis
large_data_dds <- DESeqDataSetFromMatrix(countData = counts_all, colData=category_df, design=~group)
large_data_dds <- DESeq(large_data_dds)

#extract pairwise contrasts
#nutrition
cmh <- results(large_data_dds, contrast=c("group","smallctrlmalehh","largectrlmalehh"))
dmh <- results(large_data_dds, contrast=c("group","smalldsxmalehh","largedsxmalehh"))
cmt <- results(large_data_dds, contrast=c("group","smallctrlmaleth","largectrlmaleth"))
dmt <- results(large_data_dds, contrast=c("group","smalldsxmaleth","largedsxmaleth"))
cmg <- results(large_data_dds, contrast=c("group","smallctrlmalegen","largectrlmalegen"))
dmg <- results(large_data_dds, contrast=c("group","smalldsxmalegen","largedsxmalegen"))
cmb <- results(large_data_dds, contrast=c("group","smallctrlmalebrain","largectrlmalebrain"))
dmb <- results(large_data_dds, contrast=c("group","smalldsxmalebrain","largedsxmalebrain"))

#treatment
smh <- results(large_data_dds, contrast=c("group","smallctrlmalehh","smalldsxmalehh"))
smt <- results(large_data_dds, contrast=c("group","smallctrlmaleth","smalldsxmaleth"))
smg <- results(large_data_dds, contrast=c("group","smallctrlmalegen","smalldsxmalegen"))
smb <- results(large_data_dds, contrast=c("group","smallctrlmalebrain","smalldsxmalebrain"))
lmh <- results(large_data_dds, contrast=c("group","largectrlmalehh","largedsxmalehh"))
lmt <- results(large_data_dds, contrast=c("group","largectrlmaleth","largedsxmaleth"))
lmg <- results(large_data_dds, contrast=c("group","largectrlmalegen","largedsxmalegen"))
lmb <- results(large_data_dds, contrast=c("group","largectrlmalebrain","largedsxmalebrain"))
lfh <- results(large_data_dds, contrast=c("group","largectrlfemalehh","largedsxfemalehh"))
lft <- results(large_data_dds, contrast=c("group","largectrlfemaleth","largedsxfemaleth"))
lfg <- results(large_data_dds, contrast=c("group","largectrlfemalegen","largedsxfemalegen"))
lfb <- results(large_data_dds, contrast=c("group","largectrlfemalebrain","largedsxfemalebrain"))

#sex
lch <- results(large_data_dds, contrast=c("group","largectrlfemalehh","largectrlmalehh"))
lct <- results(large_data_dds, contrast=c("group","largectrlfemaleth","largectrlmaleth"))
lcg <- results(large_data_dds, contrast=c("group","largectrlfemalegen","largectrlmalegen"))
lcb <- results(large_data_dds, contrast=c("group","largectrlfemalebrain","largectrlmalebrain"))
ldh <- results(large_data_dds, contrast=c("group","largedsxfemalehh","largedsxmalehh"))
ldt <- results(large_data_dds, contrast=c("group","largedsxfemaleth","largedsxmaleth"))
ldg <- results(large_data_dds, contrast=c("group","largedsxfemalegen","largedsxmalegen"))
ldb <- results(large_data_dds, contrast=c("group","largedsxfemalebrain","largedsxmalebrain"))

#determine number and fold change of dsxtargets
treatment.df <- data.frame(size=character(),sex=character(),tissue=character(),fc=character(),count=character(),genes=character(),stringsAsFactors=FALSE) 

#male headhorn
#extract significantly differentially expressed genes
smh_sig <- smh[complete.cases(smh$padj),]
smh_sig <- smh_sig[smh_sig$padj<=0.05,]
lmh_sig <- lmh[complete.cases(lmh$padj),]
lmh_sig <- lmh_sig[lmh_sig$padj<=0.05,]

#determine the number of significantly differentially expressed genes
n_smh <- dim(smh_sig)[1]
n_lmh <- dim(lmh_sig)[1]

#extract genes with positive or negative fold changes (i.e.: those activated or repressed by dsx)
#NB: negative fold change after dsx RNAi means they are activated by dsx.
#So n_smh_neg = the number of genes activated by dsx in small male headhorns.
smh_pos <- smh_sig[smh_sig$log2FoldChange>0,]
smh_neg <- smh_sig[smh_sig$log2FoldChange<0,]
n_smh_pos <- dim(smh_pos)[1]
n_smh_neg <- dim(smh_neg)[1]

lmh_pos <- lmh_sig[lmh_sig$log2FoldChange>0,]
lmh_neg <- lmh_sig[lmh_sig$log2FoldChange<0,]
n_lmh_pos <- dim(lmh_pos)[1]
n_lmh_neg <- dim(lmh_neg)[1]

#now look at the overlap
mh <- merge(as.data.frame(smh_sig),as.data.frame(lmh_sig),by=0)
mh_pos_s <- mh[mh$log2FoldChange.x>0,] #positive only in small
mh_neg_s <- mh[mh$log2FoldChange.x<0,] #negative only in small
mh_pos_l <- mh[mh$log2FoldChange.y>0,]
mh_neg_l <- mh[mh$log2FoldChange.y<0,]
mh_pos <- mh[mh$log2FoldChange.x>0 & mh$log2FoldChange.y>0,] #positive in both
mh_neg <- mh[mh$log2FoldChange.x<0 & mh$log2FoldChange.y<0,] #negative in both
n_mh <- dim(mh)[1]
n_mh_pos <- dim(mh_pos)[1]
n_mh_neg <- dim(mh_neg)[1]
n_mh_pos_s <- dim(mh_pos_s)[1]
n_mh_neg_s <- dim(mh_neg_s)[1]
n_mh_pos_l <- dim(mh_pos_l)[1]
n_mh_neg_l <- dim(mh_neg_l)[1]

#put the data in the new dataframe
treatment.df["smh_pos",] <- c("small","male","hh","pos",n_smh_pos,"")
treatment.df["smh_neg",] <- c("small","male","hh","neg",n_smh_neg,"")
treatment.df["lmh_pos",] <- c("large","male","hh","pos",n_lmh_pos,"")
treatment.df["lmh_neg",] <- c("large","male","hh","neg",n_lmh_neg,"")
treatment.df["mh_pos",] <- c("overlap","male","hh","pos",n_mh_pos,"")
treatment.df["mh_neg",] <- c("overlap","male","hh","neg",n_mh_neg,"")
treatment.df["mh_pos_s",] <- c("overlap-small","male","hh","pos",n_mh_pos_s,"")
treatment.df["mh_neg_s",] <- c("overlap-small","male","hh","neg",n_mh_neg_s,"")
treatment.df["mh_pos_l",] <- c("overlap-large","male","hh","pos",n_mh_pos_l,"")
treatment.df["mh_neg_l",] <- c("overlap-large","male","hh","neg",n_mh_neg_l,"")

#male thoracic horn
#extract significantly differentially expressed genes
smt_sig <- smt[complete.cases(smt$padj),]
smt_sig <- smt_sig[smt_sig$padj<=0.05,]
lmt_sig <- lmt[complete.cases(lmt$padj),]
lmt_sig <- lmt_sig[lmt_sig$padj<=0.05,]

#determine the number of significantly differentially expressed genes
n_smt <- dim(smt_sig)[1]
n_lmt <- dim(lmt_sig)[1]

#extract genes with positive or negative fold changes (i.e.: those activated or repressed by dsx)
smt_pos <- smt_sig[smt_sig$log2FoldChange>0,]
smt_neg <- smt_sig[smt_sig$log2FoldChange<0,]
n_smt_pos <- dim(smt_pos)[1]
n_smt_neg <- dim(smt_neg)[1]

lmt_pos <- lmt_sig[lmt_sig$log2FoldChange>0,]
lmt_neg <- lmt_sig[lmt_sig$log2FoldChange<0,]
n_lmt_pos <- dim(lmt_pos)[1]
n_lmt_neg <- dim(lmt_neg)[1]

#now look at the overlap
mt <- merge(as.data.frame(smt_sig),as.data.frame(lmt_sig),by=0)
mt_pos_s <- mt[mt$log2FoldChange.x>0,] #positive only in small
mt_neg_s <- mt[mt$log2FoldChange.x<0,] #negative only in small
mt_pos_l <- mt[mt$log2FoldChange.y>0,]
mt_neg_l <- mt[mt$log2FoldChange.y<0,]
mt_pos <- mt[mt$log2FoldChange.x>0 & mt$log2FoldChange.y>0,] #positive in both
mt_neg <- mt[mt$log2FoldChange.x<0 & mt$log2FoldChange.y<0,] #negative in both
n_mt <- dim(mt)[1]
n_mt_pos <- dim(mt_pos)[1]
n_mt_neg <- dim(mt_neg)[1]
n_mt_pos_s <- dim(mt_pos_s)[1]
n_mt_neg_s <- dim(mt_neg_s)[1]
n_mt_pos_l <- dim(mt_pos_l)[1]
n_mt_neg_l <- dim(mt_neg_l)[1]

#put the data in the new dataframe
treatment.df["smt_pos",] <- c("small","male","th","pos",n_smt_pos,"")
treatment.df["smt_neg",] <- c("small","male","th","neg",n_smt_neg,"")
treatment.df["lmt_pos",] <- c("large","male","th","pos",n_lmt_pos,"")
treatment.df["lmt_neg",] <- c("large","male","th","neg",n_lmt_neg,"")
treatment.df["mt_pos",] <- c("overlap","male","th","pos",n_mt_pos,"")
treatment.df["mt_neg",] <- c("overlap","male","th","neg",n_mt_neg,"")
treatment.df["mt_pos_s",] <- c("overlap-small","male","th","pos",n_mt_pos_s,"")
treatment.df["mt_neg_s",] <- c("overlap-small","male","th","neg",n_mt_neg_s,"")
treatment.df["mt_pos_l",] <- c("overlap-large","male","th","pos",n_mt_pos_l,"")
treatment.df["mt_neg_l",] <- c("overlap-large","male","th","neg",n_mt_neg_l,"")


#male genitalia
#extract significantly differentially expressed genes
smg_sig <- smg[complete.cases(smg$padj),]
smg_sig <- smg_sig[smg_sig$padj<=0.05,]
lmg_sig <- lmg[complete.cases(lmg$padj),]
lmg_sig <- lmg_sig[lmg_sig$padj<=0.05,]

#determine the number of significantly differentially expressed genes
n_smg <- dim(smg_sig)[1]
n_lmg <- dim(lmg_sig)[1]

#extract genes with positive or negative fold changes (i.e.: those activated or repressed by dsx)
smg_pos <- smg_sig[smg_sig$log2FoldChange>0,]
smg_neg <- smg_sig[smg_sig$log2FoldChange<0,]
n_smg_pos <- dim(smg_pos)[1]
n_smg_neg <- dim(smg_neg)[1]

lmg_pos <- lmg_sig[lmg_sig$log2FoldChange>0,]
lmg_neg <- lmg_sig[lmg_sig$log2FoldChange<0,]
n_lmg_pos <- dim(lmg_pos)[1]
n_lmg_neg <- dim(lmg_neg)[1]

#now look at the overlap
mg <- merge(as.data.frame(smg_sig),as.data.frame(lmg_sig),by=0)
mg_pos_s <- mg[mg$log2FoldChange.x>0,] #positive only in small
mg_neg_s <- mg[mg$log2FoldChange.x<0,] #negative only in small
mg_pos_l <- mg[mg$log2FoldChange.y>0,]
mg_neg_l <- mg[mg$log2FoldChange.y<0,]
mg_pos <- mg[mg$log2FoldChange.x>0 & mg$log2FoldChange.y>0,] #positive in both
mg_neg <- mg[mg$log2FoldChange.x<0 & mg$log2FoldChange.y<0,] #negative in both
n_mg <- dim(mg)[1]
n_mg_pos <- dim(mg_pos)[1]
n_mg_neg <- dim(mg_neg)[1]
n_mg_pos_s <- dim(mg_pos_s)[1]
n_mg_neg_s <- dim(mg_neg_s)[1]
n_mg_pos_l <- dim(mg_pos_l)[1]
n_mg_neg_l <- dim(mg_neg_l)[1]

#put the data in the new dataframe
treatment.df["smg_pos",] <- c("small","male","gen","pos",n_smg_pos,"")
treatment.df["smg_neg",] <- c("small","male","gen","neg",n_smg_neg,"")
treatment.df["lmg_pos",] <- c("large","male","gen","pos",n_lmg_pos,"")
treatment.df["lmg_neg",] <- c("large","male","gen","neg",n_lmg_neg,"")
treatment.df["mg_pos",] <- c("overlap","male","gen","pos",n_mg_pos,"")
treatment.df["mg_neg",] <- c("overlap","male","gen","neg",n_mg_neg,"")
treatment.df["mg_pos_s",] <- c("overlap-small","male","gen","pos",n_mg_pos_s,"")
treatment.df["mg_neg_s",] <- c("overlap-small","male","gen","neg",n_mg_neg_s,"")
treatment.df["mg_pos_l",] <- c("overlap-large","male","gen","pos",n_mg_pos_l,"")
treatment.df["mg_neg_l",] <- c("overlap-large","male","gen","neg",n_mg_neg_l,"")


#male brain
#extract significantly differentially expressed genes
smb_sig <- smb[complete.cases(smb$padj),]
smb_sig <- smb_sig[smb_sig$padj<=0.05,]
lmb_sig <- lmb[complete.cases(lmb$padj),]
lmb_sig <- lmb_sig[lmb_sig$padj<=0.05,]

#determine the number of significantly differentially expressed genes
n_smb <- dim(smb_sig)[1]
n_lmb <- dim(lmb_sig)[1]

#extract genes with positive or negative fold changes (i.e.: those activated or repressed by dsx)
smb_pos <- smb_sig[smb_sig$log2FoldChange>0,]
smb_neg <- smb_sig[smb_sig$log2FoldChange<0,]
n_smb_pos <- dim(smb_pos)[1]
n_smb_neg <- dim(smb_neg)[1]

lmb_pos <- lmb_sig[lmb_sig$log2FoldChange>0,]
lmb_neg <- lmb_sig[lmb_sig$log2FoldChange<0,]
n_lmb_pos <- dim(lmb_pos)[1]
n_lmb_neg <- dim(lmb_neg)[1]

#now look at the overlap
mb <- merge(as.data.frame(smb_sig),as.data.frame(lmb_sig),by=0)
mb_pos_s <- mb[mb$log2FoldChange.x>0,] #positive only in small
mb_neg_s <- mb[mb$log2FoldChange.x<0,] #negative only in small
mb_pos_l <- mb[mb$log2FoldChange.y>0,]
mb_neg_l <- mb[mb$log2FoldChange.y<0,]
mb_pos <- mb[mb$log2FoldChange.x>0 & mb$log2FoldChange.y>0,] #positive in both
mb_neg <- mb[mb$log2FoldChange.x<0 & mb$log2FoldChange.y<0,] #negative in both
n_mb <- dim(mb)[1]
n_mb_pos <- dim(mb_pos)[1]
n_mb_neg <- dim(mb_neg)[1]
n_mb_pos_s <- dim(mb_pos_s)[1]
n_mb_neg_s <- dim(mb_neg_s)[1]
n_mb_pos_l <- dim(mb_pos_l)[1]
n_mb_neg_l <- dim(mb_neg_l)[1]

#put the data in the new dataframe
treatment.df["smb_pos",] <- c("small","male","brain","pos",n_smb_pos,"")
treatment.df["smb_neg",] <- c("small","male","brain","neg",n_smb_neg,"")
treatment.df["lmb_pos",] <- c("large","male","brain","pos",n_lmb_pos,"")
treatment.df["lmb_neg",] <- c("large","male","brain","neg",n_lmb_neg,"")
treatment.df["mb_pos",] <- c("overlap","male","brain","pos",n_mb_pos,"")
treatment.df["mb_neg",] <- c("overlap","male","brain","neg",n_mb_neg,"")
treatment.df["mb_pos_s",] <- c("overlap-small","male","brain","pos",n_mb_pos_s,"")
treatment.df["mb_neg_s",] <- c("overlap-small","male","brain","neg",n_mb_neg_s,"")
treatment.df["mb_pos_l",] <- c("overlap-large","male","brain","pos",n_mb_pos_l,"")
treatment.df["mb_neg_l",] <- c("overlap-large","male","brain","neg",n_mb_neg_l,"")

#large headhorn
#extract significantly differentially expressed genes
lfh_sig <- lfh[complete.cases(lfh$padj),]
lfh_sig <- lfh_sig[lfh_sig$padj<=0.05,]

#determine the number of significantly differentially expressed genes
n_lfh <- dim(lfh_sig)[1]

#extract genes with positive or negative fold changes (i.e.: those activated or repressed by dsx)
lfh_pos <- lfh_sig[lfh_sig$log2FoldChange>0,]
lfh_neg <- lfh_sig[lfh_sig$log2FoldChange<0,]
n_lfh_pos <- dim(lfh_pos)[1]
n_lfh_neg <- dim(lfh_neg)[1]

#now look at the overlap
lh <- merge(as.data.frame(lfh_sig),as.data.frame(lmh_sig),by=0)
lh_pos_f <- lh[lh$log2FoldChange.x>0,] #positive only in female
lh_neg_f <- lh[lh$log2FoldChange.x<0,] #negative only in female
lh_pos_m <- lh[lh$log2FoldChange.y>0,]
lh_neg_m <- lh[lh$log2FoldChange.y<0,]
lh_pos <- lh[lh$log2FoldChange.x>0 & lh$log2FoldChange.y>0,] #positive in both
lh_neg <- lh[lh$log2FoldChange.x<0 & lh$log2FoldChange.y<0,] #negative in both
n_lh <- dim(lh)[1]
n_lh_pos <- dim(lh_pos)[1]
n_lh_neg <- dim(lh_neg)[1]
n_lh_pos_f <- dim(lh_pos_f)[1]
n_lh_neg_f <- dim(lh_neg_f)[1]
n_lh_pos_m <- dim(lh_pos_m)[1]
n_lh_neg_m <- dim(lh_neg_m)[1]

#put the data in the new dataframe
treatment.df["lfh_pos",] <- c("large","female","hh","pos",n_lfh_pos,"")
treatment.df["lfh_neg",] <- c("large","female","hh","neg",n_lfh_neg,"")
treatment.df["lh_pos",] <- c("large","overlap","hh","pos",n_lh_pos,"")
treatment.df["lh_neg",] <- c("large","overlap","hh","neg",n_lh_neg,"")
treatment.df["lh_pos_f",] <- c("large","overlap-female","hh","pos",n_lh_pos_f,"")
treatment.df["lh_neg_f",] <- c("large","overlap-female","hh","neg",n_lh_neg_f,"")
treatment.df["lh_pos_m",] <- c("large","overlap-male","hh","pos",n_lh_pos_m,"")
treatment.df["lh_neg_m",] <- c("large","overlap-male","hh","neg",n_lh_neg_m,"")

#large thoracic horn
#extract significantly differentially expressed genes
lft_sig <- lft[complete.cases(lft$padj),]
lft_sig <- lft_sig[lft_sig$padj<=0.05,]

#determine the number of significantly differentially expressed genes
n_lft <- dim(lft_sig)[1]

#extract genes with positive or negative fold changes (i.e.: those activated or repressed by dsx)
lft_pos <- lft_sig[lft_sig$log2FoldChange>0,]
lft_neg <- lft_sig[lft_sig$log2FoldChange<0,]
n_lft_pos <- dim(lft_pos)[1]
n_lft_neg <- dim(lft_neg)[1]

#now look at the overlap
lt <- merge(as.data.frame(lft_sig),as.data.frame(lmt_sig),by=0)
lt_pos_f <- lt[lt$log2FoldChange.x>0,] #positive only in female
lt_neg_f <- lt[lt$log2FoldChange.x<0,] #negative only in female
lt_pos_m <- lt[lt$log2FoldChange.y>0,]
lt_neg_m <- lt[lt$log2FoldChange.y<0,]
lt_pos <- lt[lt$log2FoldChange.x>0 & lt$log2FoldChange.y>0,] #positive in both
lt_neg <- lt[lt$log2FoldChange.x<0 & lt$log2FoldChange.y<0,] #negative in both
n_lt <- dim(lt)[1]
n_lt_pos <- dim(lt_pos)[1]
n_lt_neg <- dim(lt_neg)[1]
n_lt_pos_f <- dim(lt_pos_f)[1]
n_lt_neg_f <- dim(lt_neg_f)[1]
n_lt_pos_m <- dim(lt_pos_m)[1]
n_lt_neg_m <- dim(lt_neg_m)[1]

#put the data in the new dataframe
treatment.df["lft_pos",] <- c("large","female","th","pos",n_lft_pos,"")
treatment.df["lft_neg",] <- c("large","female","th","neg",n_lft_neg,"")
treatment.df["lt_pos",] <- c("large","overlap","th","pos",n_lt_pos,"")
treatment.df["lt_neg",] <- c("large","overlap","th","neg",n_lt_neg,"")
treatment.df["lt_pos_f",] <- c("large","overlap-female","th","pos",n_lt_pos_f,"")
treatment.df["lt_neg_f",] <- c("large","overlap-female","th","neg",n_lt_neg_f,"")
treatment.df["lt_pos_m",] <- c("large","overlap-male","th","pos",n_lt_pos_m,"")
treatment.df["lt_neg_m",] <- c("large","overlap-male","th","neg",n_lt_neg_m,"")


#large genitalia
#extract significantly differentially expressed genes
lfg_sig <- lfg[complete.cases(lfg$padj),]
lfg_sig <- lfg_sig[lfg_sig$padj<=0.05,]

#determine the number of significantly differentially expressed genes
n_lfg <- dim(lfg_sig)[1]

#extract genes with positive or negative fold changes (i.e.: those activated or repressed by dsx)
lfg_pos <- lfg_sig[lfg_sig$log2FoldChange>0,]
lfg_neg <- lfg_sig[lfg_sig$log2FoldChange<0,]
n_lfg_pos <- dim(lfg_pos)[1]
n_lfg_neg <- dim(lfg_neg)[1]

#now look at the overlap
lg <- merge(as.data.frame(lfg_sig),as.data.frame(lmg_sig),by=0)
lg_pos_f <- lg[lg$log2FoldChange.x>0,] #positive only in female
lg_neg_f <- lg[lg$log2FoldChange.x<0,] #negative only in female
lg_pos_m <- lg[lg$log2FoldChange.y>0,]
lg_neg_m <- lg[lg$log2FoldChange.y<0,]
lg_pos <- lg[lg$log2FoldChange.x>0 & lg$log2FoldChange.y>0,] #positive in both
lg_neg <- lg[lg$log2FoldChange.x<0 & lg$log2FoldChange.y<0,] #negative in both
n_lg <- dim(lg)[1]
n_lg_pos <- dim(lg_pos)[1]
n_lg_neg <- dim(lg_neg)[1]
n_lg_pos_f <- dim(lg_pos_f)[1]
n_lg_neg_f <- dim(lg_neg_f)[1]
n_lg_pos_m <- dim(lg_pos_m)[1]
n_lg_neg_m <- dim(lg_neg_m)[1]

#put the data in the new dataframe
treatment.df["lfg_pos",] <- c("large","female","gen","pos",n_lfg_pos,"")
treatment.df["lfg_neg",] <- c("large","female","gen","neg",n_lfg_neg,"")
treatment.df["lg_pos",] <- c("large","overlap","gen","pos",n_lg_pos,"")
treatment.df["lg_neg",] <- c("large","overlap","gen","neg",n_lg_neg,"")
treatment.df["lg_pos_f",] <- c("large","overlap-female","gen","pos",n_lg_pos_f,"")
treatment.df["lg_neg_f",] <- c("large","overlap-female","gen","neg",n_lg_neg_f,"")
treatment.df["lg_pos_m",] <- c("large","overlap-male","gen","pos",n_lg_pos_m,"")
treatment.df["lg_neg_m",] <- c("large","overlap-male","gen","neg",n_lg_neg_m,"")

#large brain
#extract significantly differentially expressed genes
lfb_sig <- lfb[complete.cases(lfb$padj),]
lfb_sig <- lfb_sig[lfb_sig$padj<=0.05,]

#determine the number of significantly differentially expressed genes
n_lfb <- dim(lfb_sig)[1]

#extract genes with positive or negative fold changes (i.e.: those activated or repressed by dsx)
lfb_pos <- lfb_sig[lfb_sig$log2FoldChange>0,]
lfb_neg <- lfb_sig[lfb_sig$log2FoldChange<0,]
n_lfb_pos <- dim(lfb_pos)[1]
n_lfb_neg <- dim(lfb_neg)[1]

#now look at the overlap
lb <- merge(as.data.frame(lfb_sig),as.data.frame(lmb_sig),by=0)
lb_pos_f <- lb[lb$log2FoldChange.x>0,] #positive in female
lb_neg_f <- lb[lb$log2FoldChange.x<0,] #negative in female
lb_pos_m <- lb[lb$log2FoldChange.y>0,]
lb_neg_m <- lb[lb$log2FoldChange.y<0,]
lb_pos <- lb[lb$log2FoldChange.x>0 & lb$log2FoldChange.y>0,] #positive in both
lb_neg <- lb[lb$log2FoldChange.x<0 & lb$log2FoldChange.y<0,] #negative in both
n_lb <- dim(lb)[1]
n_lb_pos <- dim(lb_pos)[1]
n_lb_neg <- dim(lb_neg)[1]
n_lb_pos_f <- dim(lb_pos_f)[1]
n_lb_neg_f <- dim(lb_neg_f)[1]
n_lb_pos_m <- dim(lb_pos_m)[1]
n_lb_neg_m <- dim(lb_neg_m)[1]

#put the data in the new dataframe
treatment.df["lfb_pos",] <- c("large","female","brain","pos",n_lfb_pos,"")
treatment.df["lfb_neg",] <- c("large","female","brain","neg",n_lfb_neg,"")
treatment.df["lb_pos",] <- c("large","overlap","brain","pos",n_lb_pos,"")
treatment.df["lb_neg",] <- c("large","overlap","brain","neg",n_lb_neg,"")
treatment.df["lb_pos_f",] <- c("large","overlap-female","brain","pos",n_lb_pos_f,"")
treatment.df["lb_neg_f",] <- c("large","overlap-female","brain","neg",n_lb_neg_f,"")
treatment.df["lb_pos_m",] <- c("large","overlap-male","brain","pos",n_lb_pos_m,"")
treatment.df["lb_neg_m",] <- c("large","overlap-male","brain","neg",n_lb_neg_m,"")


#save the databases to files
write.table(lmh_sig, file="/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/treatment/lmh_cd.tsv", quote=FALSE, sep="\t")
write.table(lmt_sig, file="/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/treatment/lmt_cd.tsv", quote=FALSE, sep="\t")
write.table(lmg_sig, file="/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/treatment/lmg_cd.tsv", quote=FALSE, sep="\t")
write.table(lmb_sig, file="/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/treatment/lmb_cd.tsv", quote=FALSE, sep="\t")
write.table(smh_sig, file="/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/treatment/smh_cd.tsv", quote=FALSE, sep="\t")
write.table(smt_sig, file="/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/treatment/smt_cd.tsv", quote=FALSE, sep="\t")
write.table(smg_sig, file="/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/treatment/smg_cd.tsv", quote=FALSE, sep="\t")
write.table(smb_sig, file="/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/treatment/smb_cd.tsv", quote=FALSE, sep="\t")
write.table(lfh_sig, file="/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/treatment/lfh_cd.tsv", quote=FALSE, sep="\t")
write.table(lft_sig, file="/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/treatment/lft_cd.tsv", quote=FALSE, sep="\t")
write.table(lfg_sig, file="/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/treatment/lfg_cd.tsv", quote=FALSE, sep="\t")
write.table(lfb_sig, file="/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/treatment/lfb_cd.tsv", quote=FALSE, sep="\t")

#overlap in datasets for size and sex
#x.x is small male
#y.x is large male
#x.y is large female
#y.y is large male
hh_ss <- merge(mh,lh,by="Row.names")
th_ss <- merge(mt,lt,by="Row.names")
gen_ss <- merge(mg,lg,by="Row.names")
brain_ss <- merge(mb,lb,by="Row.names")

#only hh database has genes, so let's figure them out
hh_ss_pos_sm <- hh_ss[hh_ss$log2FoldChange.x.x>0,]
hh_ss_neg_sm <- hh_ss[hh_ss$log2FoldChange.x.x<0,]
hh_ss_pos_lm <- hh_ss[hh_ss$log2FoldChange.y.x>0,]
hh_ss_neg_lm <- hh_ss[hh_ss$log2FoldChange.y.x<0,]
hh_ss_pos_lf <- hh_ss[hh_ss$log2FoldChange.x.y>0,]
hh_ss_neg_lf <- hh_ss[hh_ss$log2FoldChange.x.y<0,]
n_hh_ss <- dim(hh_ss)[1]
n_hh_ss_pos_sm <- dim(hh_ss_pos_sm)[1]
n_hh_ss_neg_sm <- dim(hh_ss_neg_sm)[1]
n_hh_ss_pos_lm <- dim(hh_ss_pos_lm)[1]
n_hh_ss_neg_lm <- dim(hh_ss_neg_lm)[1]
n_hh_ss_pos_lf <- dim(hh_ss_pos_lf)[1]
n_hh_ss_neg_lf <- dim(hh_ss_neg_lf)[1]

#but because they are needed in future calculations, here go the other tissues
#thoracic horn
th_ss_pos_sm <- th_ss[th_ss$log2FoldChange.x.x>0,]
th_ss_neg_sm <- th_ss[th_ss$log2FoldChange.x.x<0,]
th_ss_pos_lm <- th_ss[th_ss$log2FoldChange.y.x>0,]
th_ss_neg_lm <- th_ss[th_ss$log2FoldChange.y.x<0,]
th_ss_pos_lf <- th_ss[th_ss$log2FoldChange.x.y>0,]
th_ss_neg_lf <- th_ss[th_ss$log2FoldChange.x.y<0,]
n_th_ss <- dim(th_ss)[1]
n_th_ss_pos_sm <- dim(th_ss_pos_sm)[1]
n_th_ss_neg_sm <- dim(th_ss_neg_sm)[1]
n_th_ss_pos_lm <- dim(th_ss_pos_lm)[1]
n_th_ss_neg_lm <- dim(th_ss_neg_lm)[1]
n_th_ss_pos_lf <- dim(th_ss_pos_lf)[1]
n_th_ss_neg_lf <- dim(th_ss_neg_lf)[1]

#genitalia
gen_ss_pos_sm <- gen_ss[gen_ss$log2FoldChange.x.x>0,]
gen_ss_neg_sm <- gen_ss[gen_ss$log2FoldChange.x.x<0,]
gen_ss_pos_lm <- gen_ss[gen_ss$log2FoldChange.y.x>0,]
gen_ss_neg_lm <- gen_ss[gen_ss$log2FoldChange.y.x<0,]
gen_ss_pos_lf <- gen_ss[gen_ss$log2FoldChange.x.y>0,]
gen_ss_neg_lf <- gen_ss[gen_ss$log2FoldChange.x.y<0,]
n_gen_ss <- dim(gen_ss)[1]
n_gen_ss_pos_sm <- dim(gen_ss_pos_sm)[1]
n_gen_ss_neg_sm <- dim(gen_ss_neg_sm)[1]
n_gen_ss_pos_lm <- dim(gen_ss_pos_lm)[1]
n_gen_ss_neg_lm <- dim(gen_ss_neg_lm)[1]
n_gen_ss_pos_lf <- dim(gen_ss_pos_lf)[1]
n_gen_ss_neg_lf <- dim(gen_ss_neg_lf)[1]

#brain
brain_ss_pos_sm <- brain_ss[brain_ss$log2FoldChange.x.x>0,]
brain_ss_neg_sm <- brain_ss[brain_ss$log2FoldChange.x.x<0,]
brain_ss_pos_lm <- brain_ss[brain_ss$log2FoldChange.y.x>0,]
brain_ss_neg_lm <- brain_ss[brain_ss$log2FoldChange.y.x<0,]
brain_ss_pos_lf <- brain_ss[brain_ss$log2FoldChange.x.y>0,]
brain_ss_neg_lf <- brain_ss[brain_ss$log2FoldChange.x.y<0,]
n_brain_ss <- dim(brain_ss)[1]
n_brain_ss_pos_sm <- dim(brain_ss_pos_sm)[1]
n_brain_ss_neg_sm <- dim(brain_ss_neg_sm)[1]
n_brain_ss_pos_lm <- dim(brain_ss_pos_lm)[1]
n_brain_ss_neg_lm <- dim(brain_ss_neg_lm)[1]
n_brain_ss_pos_lf <- dim(brain_ss_pos_lf)[1]
n_brain_ss_neg_lf <- dim(brain_ss_neg_lf)[1]

#put this data in the dataframe
treatment.df["hh_ss_pos_sm",] <- c("overlap-small","overlap-male","hh","pos",n_hh_ss_pos_sm,"")
treatment.df["hh_ss_neg_sm",] <- c("overlap-small","overlap-male","hh","neg",n_hh_ss_neg_sm,"")
treatment.df["hh_ss_pos_lm",] <- c("overlap-large","overlap-male","hh","pos",n_hh_ss_pos_lm,"")
treatment.df["hh_ss_neg_lm",] <- c("overlap-large","overlap-male","hh","neg",n_hh_ss_neg_lm,"")
treatment.df["hh_ss_pos_lf",] <- c("overlap-large","overlap-female","hh","pos",n_hh_ss_pos_lf,"")
treatment.df["hh_ss_neg_lf",] <- c("overlap-large","overlap-female","hh","neg",n_hh_ss_neg_lf,"")

#save the treatment results table
write.table(treatment.df, file="/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/treatment-df.tsv", quote=FALSE, sep="\t")


#fill in the venn diagram
#HEADHORNS
#pos is upregulated by dsx; neg is downregulated
#large males
lmh_down <- n_lmh_neg-n_lh_neg_m-n_mh_neg_l+n_hh_ss_neg_lm
lmh_up <- n_lmh_pos-n_lh_pos_m-n_mh_pos_l+n_hh_ss_pos_lm
mh_down_l <- n_mh_neg_l-n_hh_ss_neg_lm
mh_up_l <- n_mh_pos_l-n_hh_ss_pos_lm
lh_down_m <- n_lh_neg_m-n_hh_ss_neg_lm
lh_up_m <- n_lh_pos_m-n_hh_ss_pos_lm
h_down_lm <- n_hh_ss_neg_lm
h_up_lm <- n_hh_ss_pos_lm

#small males
h_down_sm <- n_hh_ss_neg_sm
h_up_sm <- n_hh_ss_pos_sm
mh_down_s <- n_mh_neg_s-h_down_sm
mh_up_s <- n_mh_pos_s-h_up_sm
smh_down <- n_smh_neg-mh_down_s-h_down_sm
smh_up <- n_smh_pos-mh_up_s-h_up_sm

#large females
h_down_lf <- n_hh_ss_neg_lf
h_up_lf <- n_hh_ss_pos_lf
lh_down_f <- n_lh_neg_f-h_down_lf
lh_up_f <- n_lh_pos_f-h_up_lf
lfh_down <- n_lfh_neg-lh_down_f-h_down_lf
lfh_up <- n_lfh_pos-lh_up_f-h_up_lf

#THORACIC HORNS
#large males
lmt_down <- n_lmt_neg-n_lt_neg_m-n_mt_neg_l+n_th_ss_neg_lm
lmt_up <- n_lmt_pos-n_lt_pos_m-n_mt_pos_l+n_th_ss_pos_lm
mt_down_l <- n_mt_neg_l-n_th_ss_neg_lm
mt_up_l <- n_mt_pos_l-n_th_ss_pos_lm
lt_down_m <- n_lt_neg_m-n_th_ss_neg_lm
lt_up_m <- n_lt_pos_m-n_th_ss_pos_lm
t_down_lm <- n_th_ss_neg_lm
t_up_lm <- n_th_ss_pos_lm

#small males
t_down_sm <- n_th_ss_neg_sm
t_up_sm <- n_th_ss_pos_sm
mt_down_s <- n_mt_neg_s-t_down_sm
mt_up_s <- n_mt_pos_s-t_up_sm
smt_down <- n_smt_neg-mt_down_s-t_down_sm
smt_up <- n_smt_pos-mt_up_s-t_up_sm

#large females
t_down_lf <- n_th_ss_neg_lf
t_up_lf <- n_th_ss_pos_lf
lt_down_f <- n_lt_neg_f-t_down_lf
lt_up_f <- n_lt_pos_f-t_up_lf
lft_down <- n_lft_neg-lt_down_f-t_down_lf
lft_up <- n_lft_pos-lt_up_f-t_up_lf

#GENITALIA
#large males
lmg_down <- n_lmg_neg-n_lg_neg_m-n_mg_neg_l+n_gen_ss_neg_lm
lmg_up <- n_lmg_pos-n_lg_pos_m-n_mg_pos_l+n_gen_ss_pos_lm
mg_down_l <- n_mg_neg_l-n_gen_ss_neg_lm
mg_up_l <- n_mg_pos_l-n_gen_ss_pos_lm
lg_down_m <- n_lg_neg_m-n_gen_ss_neg_lm
lg_up_m <- n_lg_pos_m-n_gen_ss_pos_lm
g_down_lm <- n_gen_ss_neg_lm
g_up_lm <- n_gen_ss_pos_lm

#small males
g_down_sm <- n_gen_ss_neg_sm
g_up_sm <- n_gen_ss_pos_sm
mg_down_s <- n_mg_neg_s-g_down_sm
mg_up_s <- n_mg_pos_s-g_up_sm
smg_down <- n_smg_neg-mg_down_s-g_down_sm
smg_up <- n_smg_pos-mg_up_s-g_up_sm

#large females
g_down_lf <- n_gen_ss_neg_lf
g_up_lf <- n_gen_ss_pos_lf
lg_down_f <- n_lg_neg_f-g_down_lf
lg_up_f <- n_lg_pos_f-g_up_lf
lfg_down <- n_lfg_neg-lg_down_f-g_down_lf
lfg_up <- n_lfg_pos-lg_up_f-g_up_lf

#BRAINS
#large males
lmb_down <- n_lmb_neg-n_lb_neg_m-n_mb_neg_l+n_brain_ss_neg_lm
lmb_up <- n_lmb_pos-n_lb_pos_m-n_mb_pos_l+n_brain_ss_pos_lm
mb_down_l <- n_mb_neg_l-n_brain_ss_neg_lm
mb_up_l <- n_mb_pos_l-n_brain_ss_pos_lm
lb_down_m <- n_lb_neg_m-n_brain_ss_neg_lm
lb_up_m <- n_lb_pos_m-n_brain_ss_pos_lm
b_down_lm <- n_brain_ss_neg_lm
b_up_lm <- n_brain_ss_pos_lm

#small males
b_down_sm <- n_brain_ss_neg_sm
b_up_sm <- n_brain_ss_pos_sm
mb_down_s <- n_mb_neg_s-b_down_sm
mb_up_s <- n_mb_pos_s-b_up_sm
smb_down <- n_smb_neg-mb_down_s-b_down_sm
smb_up <- n_smb_pos-mb_up_s-b_up_sm

#large females
b_down_lf <- n_brain_ss_neg_lf
b_up_lf <- n_brain_ss_pos_lf
lb_down_f <- n_lb_neg_f-b_down_lf
lb_up_f <- n_lb_pos_f-b_up_lf
lfb_down <- n_lfb_neg-lb_down_f-b_down_lf
lfb_up <- n_lfb_pos-lb_up_f-b_up_lf


##NUTRITIONAL RESPONSE##
#headhorn
#extract significantly differentially expressed genes
cmh_sig <- cmh[complete.cases(cmh$padj),]
cmh_sig <- cmh_sig[cmh_sig$padj<=0.05,]
dmh_sig <- dmh[complete.cases(dmh$padj),]
dmh_sig <- dmh_sig[dmh_sig$padj<=0.05,]

#thoracic horn
#extract significantly differentially expressed genes
cmt_sig <- cmt[complete.cases(cmt$padj),]
cmt_sig <- cmt_sig[cmt_sig$padj<=0.05,]
dmt_sig <- dmt[complete.cases(dmt$padj),]
dmt_sig <- dmt_sig[dmt_sig$padj<=0.05,]

#genitalia
#extract significantly differentially expressed genes
cmg_sig <- cmg[complete.cases(cmg$padj),]
cmg_sig <- cmg_sig[cmg_sig$padj<=0.05,]
dmg_sig <- dmg[complete.cases(dmg$padj),]
dmg_sig <- dmg_sig[dmg_sig$padj<=0.05,]

#brain
#extract significantly differentially expressed genes
cmb_sig <- cmb[complete.cases(cmb$padj),]
cmb_sig <- cmb_sig[cmb_sig$padj<=0.05,]
dmb_sig <- dmb[complete.cases(dmb$padj),]
dmb_sig <- dmb_sig[dmb_sig$padj<=0.05,]

#overlap
cmh_dmh_merge <- merge(cmh_sig,dmh_sig,by=0)
cmt_dmt_merge <- merge(cmt_sig,dmt_sig,by=0)
cmg_dmg_merge <- merge(cmg_sig,dmg_sig,by=0)
cmb_dmb_merge <- merge(cmb_sig,dmb_sig,by=0)
dim(cmh_dmh_merge)
dim(cmt_dmt_merge)
dim(cmg_dmg_merge)
dim(cmb_dmb_merge)

#plot results
#headhorn
hnu <- dim(cmh_dmh_merge)[1] #hh nutritional response unresponsive to dsx
hnd <- dim(cmh_sig)[1] - hnu #hh nutritional response dsx mediated
hni <- dim(dmh_sig)[1] - hnu #hh nutritional response dsx inhibited

#thoracic horn
tnu <- dim(cmt_dmt_merge)[1] #th nutritional response unresponsive to dsx
tnd <- dim(cmt_sig)[1] - tnu #th nutritional response dsx mediated
tni <- dim(dmt_sig)[1] - tnu #th nutritional response dsx inhibited

#genitalia
gnu <- dim(cmg_dmg_merge)[1] #gen nutritional response unresponsive to dsx
gnd <- dim(cmg_sig)[1] - gnu #gen nutritional response dsx mediated
gni <- dim(dmg_sig)[1] - gnu #gen nutritional response dsx inhibited

#brain
bnu <- dim(cmb_dmb_merge)[1] #brain nutritional response unresponsive to dsx
bnd <- dim(cmb_sig)[1] - bnu #brain nutritional response dsx mediated
bni <- dim(dmb_sig)[1] - bnu #brain nutritional response dsx inhibited

dsxmed <- c(hnd,tnd,gnd,bnd)
dsxinh <- c(hni,tni,gni,bni)
dsxunr <- c(hnu,tnu,gnu,bnu)
nutrDf <- data.frame(dsxmed,dsxinh,dsxunr)
colnames(nutrDf) <- c("mediated","inhibited","unresponsive")
rownames(nutrDf) <-c("hh","th","gen","brain")
nutrplotdata <- as.matrix(t(nutrDf))
barplot(nutrplotdata,main="nutritional response affected by dsx",xlab="tissue",legend=colnames(nutrDf),col=c("forestgreen","brown","burlywood"),beside=T)


## just for fun, do the same with sex
#extract significantly differentially expressed genes
#hh
lch_sig <- lch[complete.cases(lch$padj),]
lch_sig <- lch_sig[lch_sig$padj<=0.05,]
ldh_sig <- ldh[complete.cases(ldh$padj),]
ldh_sig <- ldh_sig[ldh_sig$padj<=0.05,]
lch_ldh_merge <- merge(lch_sig,ldh_sig,by=0)

#th
lct_sig <- lct[complete.cases(lct$padj),]
lct_sig <- lct_sig[lct_sig$padj<=0.05,]
ldt_sig <- ldt[complete.cases(ldt$padj),]
ldt_sig <- ldt_sig[ldt_sig$padj<=0.05,]
lct_ldt_merge <- merge(lct_sig,ldt_sig,by=0)

#gen
lcg_sig <- lcg[complete.cases(lcg$padj),]
lcg_sig <- lcg_sig[lcg_sig$padj<=0.05,]
ldg_sig <- ldg[complete.cases(ldg$padj),]
ldg_sig <- ldg_sig[ldg_sig$padj<=0.05,]
lcg_ldg_merge <- merge(lcg_sig,ldg_sig,by=0)

#brain
lcb_sig <- lcb[complete.cases(lcb$padj),]
lcb_sig <- lcb_sig[lcb_sig$padj<=0.05,]
ldb_sig <- ldb[complete.cases(ldb$padj),]
ldb_sig <- ldb_sig[ldb_sig$padj<=0.05,]
lcb_ldb_merge <- merge(lcb_sig,ldb_sig,by=0)

#determine the number of significantly differentially expressed genes
n_lch <- dim(lch_sig)[1]
n_ldh <- dim(ldh_sig)[1]
n_lct <- dim(lct_sig)[1]
n_ldt <- dim(ldt_sig)[1]
n_lcg <- dim(lcg_sig)[1]
n_ldg <- dim(ldg_sig)[1]
n_lcb <- dim(lcb_sig)[1]
n_ldb <- dim(ldb_sig)[1]

#plot results
#headhorn
hsu <- dim(lch_ldh_merge)[1] #hh sex bias unresponsive to dsx
hsd <- dim(lch_sig)[1] - hsu #hh sex bias dsx mediated
hsi <- dim(ldh_sig)[1] - hsu #hh sex bias dsx inhibited

#thoracic horn
tsu <- dim(lct_ldt_merge)[1] #th sex bias unresponsive to dsx
tsd <- dim(lct_sig)[1] - tsu #th sex bias dsx mediated
tsi <- dim(ldt_sig)[1] - tsu #th sex bias dsx inhibited

#genitalia
gsu <- dim(lcg_ldg_merge)[1] #gen sex bias unresponsive to dsx
gsd <- dim(lcg_sig)[1] - gsu #gen sex bias dsx mediated
gsi <- dim(ldg_sig)[1] - gsu #gen sex bias dsx inhibited

#brain
bsu <- dim(lcb_ldb_merge)[1] #brain sex bias unresponsive to dsx
bsd <- dim(lcb_sig)[1] - bsu #brain sex bias dsx mediated
bsi <- dim(ldb_sig)[1] - bsu #brain sex bias dsx inhibited

sb_dsxmed <- c(hsd,tsd,gsd,bsd)
sb_dsxinh <- c(hsi,tsi,gsi,bsi)
sb_dsxunr <- c(hsu,tsu,gsu,bsu)
sexDf <- data.frame(sb_dsxmed,sb_dsxinh,sb_dsxunr)
colnames(sexDf) <- c("dsx-mediated","dsx-inhibited","dsx-unresponsive")
rownames(sexDf) <-c("hh","th","gen","brain")
sexplotdata <- as.matrix(t(sexDf))
barplot(sexplotdata,main="sex bias and the role of dsx",xlab="tissue",legend=colnames(sexDf),col=c("forestgreen","brown","burlywood"),beside=T)


## VENN DIAGRAMS
#headhorn
cmhgenes <- rownames(cmh_sig)
dmhgenes <- rownames(dmh_sig)
lmhgenes <- rownames(lmh_sig)
smhgenes <- rownames(smh_sig)

#genes that are dsx targeted in small and large males ONLY (not in both!)
dt.mh <- setdiff(union(lmhgenes,smhgenes),intersect(lmhgenes,smhgenes))
#genes that mediate nutritional response downstream of dsx
nrdm.mh <- setdiff(cmhgenes,dmhgenes) #differentially expressed between sm-LM in control, but not in dsx
#genes that buffer the nutritional response downstream of dsx
nrdi.mh <- setdiff(dmhgenes,cmhgenes) #differentially expressed between sm-LM in dsx, but not in ctrl
#genes that are unresponsive to dsx in the nutritional response
nrdu.mh <- intersect(dmhgenes,cmhgenes)

#overlap between above categories
unresp.allmale.mh <- intersect(dt.mh,nrdu.mh)
inhib.allmale.mh <- intersect(dt.mh,nrdi.mh)
medi.allmale.mh <- intersect(dt.mh,nrdm.mh)

#what to write in venn diagrams
#mediated by dsx, no overlap
length(nrdm.mh)-length(medi.allmale.mh)
#mediated by dsx, overlap with males
length(medi.allmale.mh)
#inhibited by dsx, no overlap
length(nrdi.mh)-length(inhib.allmale.mh)
#inhibited by dsx, overlap with males
length(inhib.allmale.mh)
#males, no overlap
length(dt.mh)-length(medi.allmale.mh)-length(inhib.allmale.mh)

#thoracic horn
cmtgenes <- rownames(cmt_sig)
dmtgenes <- rownames(dmt_sig)
lmtgenes <- rownames(lmt_sig)
smtgenes <- rownames(smt_sig)

#genes that are dsx targeted in small and large males ONLY (not in both!)
dt.mt <- setdiff(union(lmtgenes,smtgenes),intersect(lmtgenes,smtgenes))
#genes that mediate nutritional response downstream of dsx
nrdm.mt <- setdiff(cmtgenes,dmtgenes) #differentially expressed between sm-LM in control, but not in dsx
#genes that buffer the nutritional response downstream of dsx
nrdi.mt <- setdiff(dmtgenes,cmtgenes) #differentially expressed between sm-LM in dsx, but not in ctrl
#genes that are unresponsive to dsx in the nutritional response
nrdu.mt <- intersect(dmtgenes,cmtgenes)

#overlap between above categories
unresp.allmale.mt <- intersect(dt.mt,nrdu.mt)
inhib.allmale.mt <- intersect(dt.mt,nrdi.mt)
medi.allmale.mt <- intersect(dt.mt,nrdm.mt)

#what to write in venn diagrams
#mediated by dsx, no overlap
length(nrdm.mt)-length(medi.allmale.mt)
#mediated by dsx, overlap with males
length(medi.allmale.mt)
#inhibited by dsx, no overlap
length(nrdi.mt)-length(inhib.allmale.mt)
#inhibited by dsx, overlap with males
length(inhib.allmale.mt)
#males, no overlap
length(dt.mt)-length(medi.allmale.mt)-length(inhib.allmale.mt)

#genitalia
cmggenes <- rownames(cmg_sig)
dmggenes <- rownames(dmg_sig)
lmggenes <- rownames(lmg_sig)
smggenes <- rownames(smg_sig)

#genes that are dsx targeted in small and large males ONLY (not in both!)
dt.mg <- setdiff(union(lmggenes,smggenes),intersect(lmggenes,smggenes))
#genes that mediate nutritional response downstream of dsx
nrdm.mg <- setdiff(cmggenes,dmggenes) #differentially expressed between sm-LM in control, but not in dsx
#genes that buffer the nutritional response downstream of dsx
nrdi.mg <- setdiff(dmggenes,cmggenes) #differentially expressed between sm-LM in dsx, but not in ctrl
#genes that are unresponsive to dsx in the nutritional response
nrdu.mg <- intersect(dmggenes,cmggenes)

#overlap between above categories
unresp.allmale.mg <- intersect(dt.mg,nrdu.mg)
inhib.allmale.mg <- intersect(dt.mg,nrdi.mg)
medi.allmale.mg <- intersect(dt.mg,nrdm.mg)

#what to write in venn diagrams
#mediated by dsx, no overlap
length(nrdm.mg)-length(medi.allmale.mg)
#mediated by dsx, overlap with males
length(medi.allmale.mg)
#inhibited by dsx, no overlap
length(nrdi.mg)-length(inhib.allmale.mg)
#inhibited by dsx, overlap with males
length(inhib.allmale.mg)
#males, no overlap
length(dt.mg)-length(medi.allmale.mg)-length(inhib.allmale.mg)

#brain
cmbgenes <- rownames(cmb_sig)
dmbgenes <- rownames(dmb_sig)
lmbgenes <- rownames(lmb_sig)
smbgenes <- rownames(smb_sig)

#genes that are dsx targeted in small and large males ONLY (not in both!)
dt.mb <- setdiff(union(lmbgenes,smbgenes),intersect(lmbgenes,smbgenes))
#genes that mediate nutritional response downstream of dsx
nrdm.mb <- setdiff(cmbgenes,dmbgenes) #differentially expressed between sm-LM in control, but not in dsx
#genes that buffer the nutritional response downstream of dsx
nrdi.mb <- setdiff(dmbgenes,cmbgenes) #differentially expressed between sm-LM in dsx, but not in ctrl
#genes that are unresponsive to dsx in the nutritional response
nrdu.mb <- intersect(dmbgenes,cmbgenes)

#overlap between above categories
unresp.allmale.mb <- intersect(dt.mb,nrdu.mb)
inhib.allmale.mb <- intersect(dt.mb,nrdi.mb)
medi.allmale.mb <- intersect(dt.mb,nrdm.mb)

#what to write in venn diagrams
#mediated by dsx, no overlap
length(nrdm.mb)-length(medi.allmale.mb)
#mediated by dsx, overlap with males
length(medi.allmale.mb)
#inhibited by dsx, no overlap
length(nrdi.mb)-length(inhib.allmale.mb)
#inhibited by dsx, overlap with males
length(inhib.allmale.mb)
#males, no overlap
length(dt.mb)-length(medi.allmale.mb)-length(inhib.allmale.mb)


#genes of interest
smoothened <- "OTAU001567"
doublesex <- "OTAU004153"

#ask if these genes are present in certain categories...
smoothened %in% lmhgenes
smoothened %in% nrdm.mh
doublesex %in% cmhgenes



#save results
save.image(file="/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/phase1_2a_deseq2.RData")
