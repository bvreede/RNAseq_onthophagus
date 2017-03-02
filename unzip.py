import os



for file in os.listdir("/Volumes/HD1v2/dsxRNAseq/barbara/rawreads/"):
	if file[-3:] != ".gz":
		continue
	os.system("gunzip /Volumes/HD1v2/dsxRNAseq/barbara/rawreads/%s" %file)