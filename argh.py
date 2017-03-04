#!/bin/python

#Data merging and Trimmomatic trimming
#Adapted from Ed Zattara's bash scripts
#Written by Barbara Vreede

import os,sys

# Establish folders necessary
rawreads = "/Volumes/HD1v2/dsxRNAseq/barbara/rawreads"
jointreads = "/Volumes/HD3v2/barbara/jointreads"
trimmed = "/Volumes/HD3v2/barbara/trimmedreads"
meta = "/Volumes/HD1v2/dsxRNAseq/barbara/meta"


samplenos = []
filenames = {}

#Go through files and unzip and change names if necessary.
#This is also used to collect the sample numbers
for i in os.listdir(rawreads):
	if i[0] == '.':
		continue
	gsf = i.split('-')[0]
	#remove A/B/C from GSF number
	gsf = gsf[:-1]
	#determine the sample number
	fileno = i.split('_S')[1].split('_')[0] #sample number
	if gsf == 'GSF1120': #GSF1120 is different, because there are two numbers in the filename, so this needs to be fixed.
		continue
	samplenos.append(fileno)
	filenames[fileno] = i

#make the sample list unique
samplenos = list(set(samplenos))
print "Total raw samples:", len(samplenos)

#Merge sequence files from the same sample.
for n in samplenos:
	# First, a very elaborate sequence of finding the correct name for the cat file
	# determine sample name as used in filename
	sname = "_S%s_" %n
	# find filename
	fname = filenames[n].lower()
	# determine sex and treatment from filename
	if 'dsxm' in fname:
		treatment = "dsx"
		sex = "male"
	elif 'dsxf' in fname:
		treatment = "dsx"
		sex = "female"
	elif 'ctrlm' in fname:
		treatment = "ctrl"
		sex = "male"
	elif 'ctrlf' in fname:
		treatment = "ctrl"
		sex = "female"
	else:
		print "Could not determine sex/treatment for file %s (sample S%s)." %(fname,n)
		continue
	# determine tissue from filename
	if '-che' in fname:
		tissue = 'HH'
	elif '-the' in fname:
		tissue = 'TH'
	elif '-br' in fname:
		tissue = 'BR'
	elif '-gen' in fname:
		tissue = 'GEN'
	else:
		print "Could not determine tissue for file %s (sample S%s)." %(fname,n)
		continue
	# determine size from filename !! FROM PROJECT ID!
	if "863" in fname:
		size = 'large'
	elif "1120" in fname:
		size = 'small'
	elif "1287" in fname:
		size = 'small'
	else:
		print "Could not determine size for file %s (sample S%s)." %(fname,n)
		continue
	
	# concatenate all files from this sample to the new file.
	jointname = "%s/%s_%s_%s_%s_S%s.fastq" %(jointreads,size,sex,treatment,tissue,n)
	
	if os.path.exists(jointname):
		print "The file for sample S%s already exists." %n
	else:
		os.system("touch %s" %jointname)