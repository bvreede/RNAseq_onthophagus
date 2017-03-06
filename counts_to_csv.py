#!/bin/python

import os, sys, datetime

#input and output
htseqdir = "/Volumes/HD3v2/barbara/htseq_counts"
outfile = "/Volumes/HD3v2/barbara/countstable.csv"

#in case a counts table already exists, make sure a copy is saved.
if os.path.exists(outfile): 
	rightnow = str(datetime.datetime.now()).replace(':','-').replace(' ','_').split('.')[0]
	os.system("cp %s %s_copy_%s.csv" %(outfile,outfile[:-4],rightnow))

out = open(outfile,"w")

#make a list of gene IDs, from OTAU000001 to OTAU017483
allgenes = []
for n in range(1,17484):
	gene = "OTAU" + str("%06d"%n)
	allgenes.append(gene)

#make a dataframe that will contain the reads
alldata = [['']] #empty list, containing only the first element, which is also a list.
for g in allgenes:
	alldata.append([g]) #for each gene, add the gene name as first element inside a list to the larger list

#fill the dataframe with contents of each file in htseqdir	
for file in os.listdir(htseqdir):
	filesplit = file.split('_trimmed')
	if len(filesplit) != 2: #ensure that only proper htseq outputfiles are used
		continue
	sname = filesplit[0]
	#open the file and read the contents
	sample = open("%s/%s" %(htseqdir,file))
	#put the contents in a dictionary
	sdict = {}
	for line in sample:
		gene,read = line.split()
		sdict[gene] = read
	#write the dictionary to the 'alldata' table
	for d in alldata:
		gene = d[0]
		if gene == "":
			d.append(sname) #adds the sample name to the first element in the list
		else:
			d.append(sdict[gene])


#write the dataframe to the outputfile
for line in alldata:
	outstr = ""
	for e in line:
		outstr += e
		outstr += ','
	out.write(outstr[:-1])
	out.write('\n')
