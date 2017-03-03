#!/bin/bash
#batch tophat alignment script

PATH=$PATH:/Users/bioinformatics/bin/bowtie2-2.2.5
export PATH

tophatdir=/Users/bioinformatics/bin/tophat-2.0.14.OSX_x86_64
indexdir=/Volumes/HD1v2/otau/Onthophagus_Tophat_Index
outputdir=/Volumes/HD3v2/barbara/tophat_alignments
inputdir=/Volumes/HD3v2/barbara/trimmedreads

cd $inputdir

for file in *; do
echo "Starting alignment of $file"
$tophatdir/tophat2 -o $outputdir/$file --transcriptome-index=$indexdir/transcriptome-index -p8 --b2-very-sensitive --read-edit-dist 2 $indexdir/Otaur.scaffolds $file
echo "Finished aligning $file!"

done