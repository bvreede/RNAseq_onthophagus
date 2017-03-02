#!/bin/bash
#batch htseq alignment script

PATH=$PATH:/Users/bioinformatics/bin/bowtie2-2.2.5
export PATH

tophatdir=/Users/bioinformatics/bin/tophat-2.0.14.OSX_x86_64
gffdir=/Volumes/HD1/GenomicResources/OtauGenome/i5kdata/scaffold/analyses/BCM_version_0.5.3/consensus_gene_set
outputdir=/Volumes/HD1/otau/i5kalignments/htseq_counts
inputdir=/Volumes/HD1/otau/i5kalignments

cd $inputdir

for file in WA*; do

echo "Starting alignment of $file"

htseq-count -f bam -r pos -s no -t gene -i ID  $file/accepted_hits_sorted.bam $gffdir/OTAU.Models.gff3 > $outputdir/$file.counts

echo "Finished aligning $file!"

done