#install and import bioclite and deseq2 (only needs to be done once)
source("http://bioconductor.org/biocLite.R")
biocLite("DESeq2")

#load deseq2
library(DESeq2)

#read the database
counts_all <- read.table("/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/counts_all.tsv",header=T,row.names=1)

#make a matrix of categories based on the header names
all_names <- as.vector(colnames(counts_all))
temp_li <- strsplit(all_names, "_")
temp_mat <- matrix(unlist(temp_li), ncol=5,byrow=TRUE)
category_df <- as.data.frame(temp_mat)
colnames(category_df) <- c("size","treatment","sex","tissue","sample")
rownames(category_df) <- all_names

### FUNCTION FOR DATA SELECTION ###
#use: selection(size="",treatment="",sex="",tissue=""), using "all" instead of the descriptor in the variable for which you want to assess differential expression.
#e.g. select_lmh <- selection(size="large",sex="male",tissue="hh")
#select_lmh$cols gives the category/colname database; select_lmh$counts gives the counttable.

selection <- function(size="all",treatment="all",sex="all",tissue="all"){
#select the comparison you want to make:
#select size
if(size=="all"){counts_s <- counts_all
category_s <- category_df}
else{counts_s <- counts_all[,grep(size,colnames(counts_all))]
category_s <- category_df[category_df$size==size,]}
#select treatment
if(treatment=="all"){counts_st <- counts_s
category_st <- category_s}
else{counts_st <- counts_s[,grep(treatment,colnames(counts_s))]
category_st <- category_s[category_s$treatment==treatment,]}
#select sex
if(sex=="all"){counts_sts <- counts_st
category_sts <- category_st}
else{sex.adj <- paste0("_",sex,"_")
counts_sts<-counts_st[,grep(sex.adj,colnames(counts_st))]
category_sts <- category_st[category_st$sex==sex,]}
#select tissue
counts_stst<-counts_sts[,grep(tissue,colnames(counts_sts))]
category_stst <- category_sts[category_sts$tissue==tissue,]
#return category and counts dataframes
returnlist <- list("counts" = counts_stst, "cols" = category_stst)
return(returnlist)
}

### FUNCTION FOR RUNNING DESEQ2 WITH CUSTOMIZABLE COMPARISON ###
#use: deseq_automator(comparison="variable", selection="all")
#e.g. deseq_automator(comparison="treatment",selection="lmh")

deseq_automator <- function(comparison="comparison",selection=""){
#generate the first deseq dataset, entering the comparison you want to look at
if(comparison=="treatment"){
dds_selected<-DESeqDataSetFromMatrix(
countData=data_selection$counts,
colData=data_selection$cols,
design= ~treatment)}
if(comparison=="size"){
dds_selected<-DESeqDataSetFromMatrix(
countData=data_selection$counts,
colData=data_selection$cols,
design= ~size)}
if(comparison=="sex"){
dds_selected<-DESeqDataSetFromMatrix(
countData=data_selection$counts,
colData=data_selection$cols,
design= ~sex)}
#put the now generated dataset into dds, and remove rows that have less than one read count
dds <- dds_selected
dds <- dds[rowSums(counts(dds))>1,]

#determine the appropriate order for fold change depending on the comparison (variable mentioned is first, so associated with negative fold changes)
if(comparison=="size"){
dds$size <- relevel(dds$size,"small")}
if(comparison=="sex"){
dds$sex <- relevel(dds$sex,"female")}
if(comparison=="treatment"){
dds$treatment <- relevel(dds$treatment,"ctrl")}

#run deseq and look at the results
dds <- DESeq(dds)
res <- results(dds)

#select the appropriate genes: first, only those with non-NA p values; then: those with <0.05 adjusted p value (padj)
res <-res[complete.cases(res$padj),]
res_sig <- res[res$padj<=0.05,]

#write the results to file (and add a header to the row names)
Gene_ID <- as.data.frame(rownames(res_sig))
colnames(Gene_ID) <- "Gene_ID"
res_sig <- as.data.frame(cbind(Gene_ID,res_sig))

#define the name for the result file
outfile <- paste0("/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/",comparison,selection,".tsv")

write.table(res_sig,file=outfile,row.names=FALSE,sep="\t",quote=FALSE)
}


### THIS IS WHERE THE CODE CAN BE CUSTOMIZED BASED ON EACH COMPARISON ###

#use the above functions:

##TREATMENT##
#large-male-headhorn
data_selection <- selection(size="large",treatment="all",sex="male",tissue="hh")
deseq_automator(comparison="treatment",selection="lmh")

#large-male-thoracic
data_selection <- selection(size="large",treatment="all",sex="male",tissue="th")
deseq_automator(comparison="treatment",selection="lmt")

#large-male-brain
data_selection <- selection(size="large",treatment="all",sex="male",tissue="brain")
deseq_automator(comparison="treatment",selection="lmb")

#large-male-genitalia
data_selection <- selection(size="large",treatment="all",sex="male",tissue="gen")
deseq_automator(comparison="treatment",selection="lmg")

#large-female-headhorn
data_selection <- selection(size="large",treatment="all",sex="female",tissue="hh")
deseq_automator(comparison="treatment",selection="lfh")

#large-female-thoracic
data_selection <- selection(size="large",treatment="all",sex="female",tissue="th")
deseq_automator(comparison="treatment",selection="lft")

#large-female-brain
data_selection <- selection(size="large",treatment="all",sex="female",tissue="brain")
deseq_automator(comparison="treatment",selection="lfb")

#large-female-genitalia
data_selection <- selection(size="large",treatment="all",sex="female",tissue="gen")
deseq_automator(comparison="treatment",selection="lfg")

#small-male-headhorn
data_selection <- selection(size="small",treatment="all",sex="male",tissue="hh")
deseq_automator(comparison="treatment",selection="smh")

#small-male-thoracic
data_selection <- selection(size="small",treatment="all",sex="male",tissue="th")
deseq_automator(comparison="treatment",selection="smt")

#small-male-brain
data_selection <- selection(size="small",treatment="all",sex="male",tissue="brain")
deseq_automator(comparison="treatment",selection="smb")

#small-male-genitalia
data_selection <- selection(size="small",treatment="all",sex="male",tissue="gen")
deseq_automator(comparison="treatment",selection="smg")

##CONDITION##
#control-male-headhorn
data_selection <- selection(size="all",treatment="ctrl",sex="male",tissue="hh")
deseq_automator(comparison="size",selection="cmh")

#control-male-thoracic
data_selection <- selection(size="all",treatment="ctrl",sex="male",tissue="th")
deseq_automator(comparison="size",selection="cmt")

#control-male-brain
data_selection <- selection(size="all",treatment="ctrl",sex="male",tissue="brain")
deseq_automator(comparison="size",selection="cmb")

#control-male-genitalia
data_selection <- selection(size="all",treatment="ctrl",sex="male",tissue="gen")
deseq_automator(comparison="size",selection="cmg")

#dsxRNAI-male-headhorn
data_selection <- selection(size="all",treatment="dsx",sex="male",tissue="hh")
deseq_automator(comparison="size",selection="dmh")

#dsxRNAI-male-thoracic
data_selection <- selection(size="all",treatment="dsx",sex="male",tissue="th")
deseq_automator(comparison="size",selection="dmt")

#dsxRNAI-male-brain
data_selection <- selection(size="all",treatment="dsx",sex="male",tissue="brain")
deseq_automator(comparison="size",selection="dmb")

#dsxRNAI-male-genitalia
data_selection <- selection(size="all",treatment="dsx",sex="male",tissue="gen")
deseq_automator(comparison="size",selection="dmg")

##SEX##
#large-control-headhorn
data_selection <- selection(size="large",treatment="ctrl",sex="all",tissue="hh")
deseq_automator(comparison="sex",selection="lch")

#large-control-thoracic
data_selection <- selection(size="large",treatment="ctrl",sex="all",tissue="th")
deseq_automator(comparison="sex",selection="lct")

#large-control-brain
data_selection <- selection(size="large",treatment="ctrl",sex="all",tissue="brain")
deseq_automator(comparison="sex",selection="lcb")

#large-control-genitalia
data_selection <- selection(size="large",treatment="ctrl",sex="all",tissue="gen")
deseq_automator(comparison="sex",selection="lcg")

#large-dsxRNAi-headhorn
data_selection <- selection(size="large",treatment="dsx",sex="all",tissue="hh")
deseq_automator(comparison="sex",selection="ldh")

#large-dsxRNAi-thoracic
data_selection <- selection(size="large",treatment="dsx",sex="all",tissue="th")
deseq_automator(comparison="sex",selection="ldt")

#large-dsxRNAi-brain
data_selection <- selection(size="large",treatment="dsx",sex="all",tissue="brain")
deseq_automator(comparison="sex",selection="ldb")

#large-dsxRNAi-genitalia
data_selection <- selection(size="large",treatment="dsx",sex="all",tissue="gen")
deseq_automator(comparison="sex",selection="ldg")