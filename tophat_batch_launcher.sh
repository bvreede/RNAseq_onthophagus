#!/bin/bash
#batch tophat alignment script

#add bowtie and samtools to path
PATH=$PATH:/Users/bioinformatics/bin/bowtie2-2.2.5:/Users/bioinformatics/bin/samtools-1.3.1
export PATH

tophatdir=/Users/bioinformatics/bin/tophat-2.0.14.OSX_x86_64
indexdir=/Volumes/HD1v2/otau/Onthophagus_Tophat_Index
outputdir=/Volumes/HD3v2/barbara/tophat_alignments_smallF
inputdir=/Volumes/HD1v2/dsxRNAseq/barbara/trimmedreads_smallF

cd $inputdir

for file in *; do
echo "Starting alignment of $file..."
$tophatdir/tophat2 -o $outputdir/$file --transcriptome-index=$indexdir/transcriptome-index -p8 --b2-very-sensitive --read-edit-dist 2 $indexdir/Otaur.scaffolds $file
echo "Alignment of $file complete."

done