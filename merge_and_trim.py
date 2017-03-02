#!/bin/python

#Data merging and Trimmomatic trimming
#Adapted from Ed Zattara's bash scripts
#Written by Barbara Vreede

import os,sys

# Establish folders necessary
rawreads = "/Volumes/HD1V2/dsxRNAseq/barbara/rawreads"
jointreads = "/Volumes/HD1V2/dsxRNAseq/barbara/jointreads"
trimmed = "/Volumes/HD1V2/dsxRNAseq/barbara/trimmedreads"
meta = "/Volumes/HD1V2/dsxRNAseq/barbara/meta"

if not os.path.exists(rawreads):
	sys.exit("Input folder 'rawreads' (%s) is not available. Check the path and try again." %rawreads)

if not os.path.exists(jointreads):
	os.system("mkdir %s" %jointreads)

if not os.path.exists(trimmed):
	os.system("mkdir %s" %trimmed)
	
if not os.path.exists(meta):
	os.system("mkdir %s" %meta)

#Translation meta file
gsftf_fn = "%s/GSF1120_sampletranslation.csv"
if os.path.exists("%s/%s" %(meta,gsftf_fn)):
	gsftf = open("%s/GSF1120_sampletranslation.csv" %meta,"a")
else:
	gsftf = open("%s/GSF1120_sampletranslation.csv" %meta,"w")
	gsftf.write("filename,new_filename,sampleID,new_sampleID\n")

samplenos = []
filenames = {}

#Go through files and unzip and change names if necessary.
#This is also used to collect the sample numbers
for i in os.listdir(rawreads):
	if i[0] == '.':
		continue
	#unzip if necessary
	if i[-3:] == ".gz":
		print "Unzipping", i
		os.system("gunzip %s/%s" %(rawreads,i))
		i = i.replace('.gz','') #update filename
	gsf = i.split('-')[0]
	if gsf == 'GSF':
		gsf = 'GSF'+i.split('-')[1]
	#remove A/B/C from GSF number
	gsf = gsf[:-1]
	#determine the sample number
	fileno = i.split('_S')[1].split('_')[0] #sample number
	if gsf == 'GSF1120': #GSF1120 is different, because there are two numbers in the filename
		subfileno = i.split('_S')[0].split('-')[-1] 
		# create a filename so that only one sample ID is written in the filename
		newi = i.replace('-' + subfileno,'')
		subfileno = subfileno.replace('L','')
		try:
			int(subfileno)
		except ValueError: # this means the subfilenumber is not a number, which means the filename has already been changed.
			continue
		newi = newi.replace('_S%s_' %fileno,'_S%s_' %subfileno)
		# put the old and new name in the meta database
		gsftf.write("%s,%s,%s,%s\n" %(i,newi,fileno,subfileno))
		print "Changing filename for", i
		os.system("mv %s/%s %s/%s" %(rawreads,i,rawreads,newi)) #change the name of the file
		# apply the new sample number
		fileno = subfileno
	samplenos.append(fileno)
	filenames[fileno] = i
gsftf.close()

#make the sample list unique
samplenos = list(set(samplenos))

#Merge sequence files from the same sample.
for n in samplenos:
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
	print fname, sname, treatment, sex