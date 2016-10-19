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

#save results
save.image(file="/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/phase1_2a_deseq2.RData")
