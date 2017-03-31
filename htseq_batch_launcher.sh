#!/bin/bash
#batch htseq alignment script


#add bowtie and samtools to path
PATH=$PATH:/Users/bioinformatics/bin/bowtie2-2.2.5:/Users/bioinformatics/bin/samtools-1.3.1
export PATH

tophatdir=/Users/bioinformatics/bin/tophat-2.0.14.OSX_x86_64
gff=/Volumes/HD1v2/otau/Onthophagus_Tophat_Index/OTAU.Models.gff3
outputdir=/Volumes/HD3v2/barbara/htseq_counts
inputdir=/Volumes/HD3v2/barbara/tophat_alignments

cd $inputdir

for file in small_female*; do

echo "Starting count analysis for $file..."

htseq-count -f bam -r pos -s no -t gene -i ID  $file/accepted_hits.bam $gff > $outputdir/$file.counts

echo "Count analysis of $file complete."

done