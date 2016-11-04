
import os
#import pandas as pd
sourcedir = "/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/Onthophagus_FoldChange"

"""
for f in os.listdir(sourcedir):
	if ".annot.tsv" not in f:
		continue
	
	infile = open("%s/%s" %(sourcedir,f))
	outfile = open("%s/%s" %(sourcedir,f.replace(".tsv",".csv")), "w")
	count=0
	for line in infile:
		lili = line.split('\t')
		for l in lili[:-1]:
			if "thioadenosine phosphorylase" in l:
				l = "S-methyl-5-thioadenosine phosphorylase" #correction due to weird character that does not read well in R!
			outfile.write("%s;"%l.strip())
		outfile.write("%s\n" %lili[-1].strip())
		count +=1
	outfile.close()
	print "number of lines in file %s:" %f, count
"""
"""	
infile = open("/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/OnthophagusFPKM.annot.tsv")
outfile = open("/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/gene_annotation_table.csv", "w")


count=0
for line in infile:
	lili = line.split('\t')
	for l in lili[:6]:
		l=l.replace(' ','_')
		l=l.replace("'","")
		l=l.replace(",","")
		outfile.write("%s;"%l.strip())
	outfile.write("%s\n" %lili[6].strip())
	count +=1
outfile.close()
print "total no. rows:", count
"""

"""	
outfile = open("/Users/BarbaraMaria/Box Sync/PROJECT doublesex/Scripts/github-barbara/onthophagus/to_R.txt","w")
for f in os.listdir(sourcedir):
	if ".annot.csv" not in f:
		continue
	if "DESeq" not in f:
		continue
	outfile.write('res <- read.csv("%s/%s", sep=";")\n' %(sourcedir,f))
	outfile.write('res <- res[complete.cases(res$padj),]\n')
	outfile.write('res_sig <- res[res$padj <= 0.05,]\n')
	outfile.write('dim(res_sig)\n')
"""

#gofile = pd.DataFrame.from_csv("/Users/BarbaraMaria/Box Sync/Onthophagus/OtauGenomicData/otau_blastP_uniprot_GO.txt", sep="\t")
#gofile.to_csv("/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/Otau-genedata/godata2.tsv", sep='\t')
#print gofile["Gene.ontology..molecular.function."][371:377]

def problemfixer(l):
	l = l.replace("5'","5pr")
	l = l.replace("3'","3pr")
	l = l.replace("2'","2pr")
	#l = l.replace("-PA","")
	#l = l.replace("; ","+")
	#l = l.replace(", ","+")
	l=l.replace("'","") #Jerdon's jumping ant, Anna's something or other
	return l
"""
goin= open("/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/Otau-genedata/godata2.tsv")
goout = open("/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/Otau-genedata/godata3.tsv","w")
count = 0
dbdict = {}
for line in goin:
	count +=1
	lili = line.split('\t')
	dbdict[count] = lili
	for l in lili[:-1]:
		l = problemfixer(l)
		goout.write("%s\t"%l.strip())
	lastl = problemfixer(lili[-1].strip())
	goout.write("%s\n" %lastl)

goout.close()
print dbdict[10214]
"""

"""
genein = open("/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/Otau-genedata/Otau_genedata.tsv")
geneout = open("/Users/BarbaraMaria/Dropbox/projects/2016_dsx-phase2/data/Otau-genedata/Otau_genedata2.tsv","w")

count = 0
dbdict = {}
for line in genein:
	count +=1
	lili = line.split('\t')
	dbdict[count] = lili
	for l in lili[:-1]:
		l = problemfixer(l)
		geneout.write("%s\t"%l.strip())
	lastl = problemfixer(lili[-1].strip())
	geneout.write("%s\n" %lastl)

geneout.close()
"""


