#!/bin/bash
#Data merging and Trimmomatic trimming
#Should process 

rawreads=/Volumes/HD2/dsx_RNAseq_phase2/rawreads
jointreads=/Volumes/HD2/dsx_RNAseq_phase2/jointreads
trimmed=/Volumes/HD2/dsx_RNAseq_phase2/trimmedreads

# Samples 97 to 102 - dsxM BRN
gzcat -d $rawreads/GSF1120*97* > $jointreads/GSF1120-097-dsxM-BRN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-097-dsxM-BRN_R1_001.fastq $trimmed/GSF1120-097-dsxM-BRN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*98* > $jointreads/GSF1120-098-dsxM-BRN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-098-dsxM-BRN_R1_001.fastq $trimmed/GSF1120-098-dsxM-BRN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*99* > $jointreads/GSF1120-099-dsxM-BRN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-099-dsxM-BRN_R1_001.fastq $trimmed/GSF1120-099-dsxM-BRN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*100* > $jointreads/GSF1120-100-dsxM-BRN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-100-dsxM-BRN_R1_001.fastq $trimmed/GSF1120-100-dsxM-BRN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*101* > $jointreads/GSF1120-101-dsxM-BRN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-101-dsxM-BRN_R1_001.fastq $trimmed/GSF1120-101-dsxM-BRN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*102* > $jointreads/GSF1120-102-dsxM-BRN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-102-dsxM-BRN_R1_001.fastq $trimmed/GSF1120-102-dsxM-BRN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Samples 97 to 102 - dsxM BRN done"

# Samples 103 to 108 - ctrlM BRN
gzcat -d $rawreads/GSF1120*103* > $jointreads/GSF1120-103-ctrlM-BRN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-103-ctrlM-BRN_R1_001.fastq $trimmed/GSF1120-103-ctrlM-BRN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*104* > $jointreads/GSF1120-104-ctrlM-BRN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-104-ctrlM-BRN_R1_001.fastq $trimmed/GSF1120-104-ctrlM-BRN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*105* > $jointreads/GSF1120-105-ctrlM-BRN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-105-ctrlM-BRN_R1_001.fastq $trimmed/GSF1120-105-ctrlM-BRN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*106* > $jointreads/GSF1120-106-ctrlM-BRN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-106-ctrlM-BRN_R1_001.fastq $trimmed/GSF1120-106-ctrlM-BRN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*107* > $jointreads/GSF1120-107-ctrlM-BRN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-107-ctrlM-BRN_R1_001.fastq $trimmed/GSF1120-107-ctrlM-BRN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*108* > $jointreads/GSF1120-108-ctrlM-BRN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-108-ctrlM-BRN_R1_001.fastq $trimmed/GSF1120-108-ctrlM-BRN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Samples 103 to 108 - ctrlM BRN"

# Samples 109 to 114 - dsxM THE
gzcat -d $rawreads/GSF1120*109* > $jointreads/GSF1120-109-dsxM-THE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-109-dsxM-THE_R1_001.fastq $trimmed/GSF1120-109-dsxM-THE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*110* > $jointreads/GSF1120-110-dsxM-THE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-110-dsxM-THE_R1_001.fastq $trimmed/GSF1120-110-dsxM-THE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*111* > $jointreads/GSF1120-111-dsxM-THE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-111-dsxM-THE_R1_001.fastq $trimmed/GSF1120-111-dsxM-THE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*112* > $jointreads/GSF1120-112-dsxM-THE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-112-dsxM-THE_R1_001.fastq $trimmed/GSF1120-112-dsxM-THE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*113* > $jointreads/GSF1120-113-dsxM-THE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-113-dsxM-THE_R1_001.fastq $trimmed/GSF1120-113-dsxM-THE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*114* > $jointreads/GSF1120-114-dsxM-THE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-114-dsxM-THE_R1_001.fastq $trimmed/GSF1120-114-dsxM-THE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Samples 109 to 114 - dsxM THE"

# Samples 115 to 120 - ctrlM THE
gzcat -d $rawreads/GSF1120*115* > $jointreads/GSF1120-115-ctrlM-THE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-115-ctrlM-THE_R1_001.fastq $trimmed/GSF1120-115-ctrlM-THE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*116* > $jointreads/GSF1120-116-ctrlM-THE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-116-ctrlM-THE_R1_001.fastq $trimmed/GSF1120-116-ctrlM-THE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*117* > $jointreads/GSF1120-117-ctrlM-THE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-117-ctrlM-THE_R1_001.fastq $trimmed/GSF1120-117-ctrlM-THE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*118* > $jointreads/GSF1120-118-ctrlM-THE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-118-ctrlM-THE_R1_001.fastq $trimmed/GSF1120-118-ctrlM-THE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*119* > $jointreads/GSF1120-119-ctrlM-THE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-119-ctrlM-THE_R1_001.fastq $trimmed/GSF1120-119-ctrlM-THE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*120* > $jointreads/GSF1120-120-ctrlM-THE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-120-ctrlM-THE_R1_001.fastq $trimmed/GSF1120-120-ctrlM-THE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Samples 115 to 120 - ctrlM THE"

# Samples 121 to 126 - dsxM CHE
gzcat -d $rawreads/GSF1120*121* > $jointreads/GSF1120-121-dsxM-CHE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-121-dsxM-CHE_R1_001.fastq $trimmed/GSF1120-121-dsxM-CHE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*122* > $jointreads/GSF1120-122-dsxM-CHE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-122-dsxM-CHE_R1_001.fastq $trimmed/GSF1120-122-dsxM-CHE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*123* > $jointreads/GSF1120-123-dsxM-CHE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-123-dsxM-CHE_R1_001.fastq $trimmed/GSF1120-123-dsxM-CHE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*124* > $jointreads/GSF1120-124-dsxM-CHE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-124-dsxM-CHE_R1_001.fastq $trimmed/GSF1120-124-dsxM-CHE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*125* > $jointreads/GSF1120-125-dsxM-CHE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-125-dsxM-CHE_R1_001.fastq $trimmed/GSF1120-125-dsxM-CHE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*126* > $jointreads/GSF1120-126-dsxM-CHE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-126-dsxM-CHE_R1_001.fastq $trimmed/GSF1120-126-dsxM-CHE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Samples 121 to 126 - dsxM CHE"

# Samples 127 to 132 - ctrlM CHE
gzcat -d $rawreads/GSF1120*127* > $jointreads/GSF1120-127-ctrlM-CHE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-127-ctrlM-CHE_R1_001.fastq $trimmed/GSF1120-127-ctrlM-CHE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*128* > $jointreads/GSF1120-128-ctrlM-CHE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-128-ctrlM-CHE_R1_001.fastq $trimmed/GSF1120-128-ctrlM-CHE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*129* > $jointreads/GSF1120-129-ctrlM-CHE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-129-ctrlM-CHE_R1_001.fastq $trimmed/GSF1120-129-ctrlM-CHE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*130* > $jointreads/GSF1120-130-ctrlM-CHE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-130-ctrlM-CHE_R1_001.fastq $trimmed/GSF1120-130-ctrlM-CHE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*131* > $jointreads/GSF1120-131-ctrlM-CHE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-131-ctrlM-CHE_R1_001.fastq $trimmed/GSF1120-131-ctrlM-CHE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*132* > $jointreads/GSF1120-132-ctrlM-CHE_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-132-ctrlM-CHE_R1_001.fastq $trimmed/GSF1120-132-ctrlM-CHE_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Samples 127 to 132 - ctrlM CHE"

# Samples 133 to 138 - dsxM GEN
gzcat -d $rawreads/GSF1120*133* > $jointreads/GSF1120-133-dsxM-GEN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-133-dsxM-GEN_R1_001.fastq $trimmed/GSF1120-133-dsxM-GEN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*134* > $jointreads/GSF1120-134-dsxM-GEN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-134-dsxM-GEN_R1_001.fastq $trimmed/GSF1120-134-dsxM-GEN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*135* > $jointreads/GSF1120-135-dsxM-GEN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-135-dsxM-GEN_R1_001.fastq $trimmed/GSF1120-135-dsxM-GEN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*136* > $jointreads/GSF1120-136-dsxM-GEN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-136-dsxM-GEN_R1_001.fastq $trimmed/GSF1120-136-dsxM-GEN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*137* > $jointreads/GSF1120-137-dsxM-GEN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-137-dsxM-GEN_R1_001.fastq $trimmed/GSF1120-137-dsxM-GEN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*138* > $jointreads/GSF1120-138-dsxM-GEN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-138-dsxM-GEN_R1_001.fastq $trimmed/GSF1120-138-dsxM-GEN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Samples 133 to 138 - dsxM GEN"

# Samples 139 to 142 - ctrlM GEN
gzcat -d $rawreads/GSF1120*139* > $jointreads/GSF1120-139-ctrlM-GEN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-139-ctrlM-GEN_R1_001.fastq $trimmed/GSF1120-139-ctrlM-GEN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*140* > $jointreads/GSF1120-140-ctrlM-GEN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-140-ctrlM-GEN_R1_001.fastq $trimmed/GSF1120-140-ctrlM-GEN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*141* > $jointreads/GSF1120-141-ctrlM-GEN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-141-ctrlM-GEN_R1_001.fastq $trimmed/GSF1120-141-ctrlM-GEN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*142* > $jointreads/GSF1120-142-ctrlM-GEN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-142-ctrlM-GEN_R1_001.fastq $trimmed/GSF1120-142-ctrlM-GEN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*143* > $jointreads/GSF1120-143-ctrlM-GEN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-143-ctrlM-GEN_R1_001.fastq $trimmed/GSF1120-143-ctrlM-GEN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

gzcat -d $rawreads/GSF1120*144* > $jointreads/GSF1120-144-ctrlM-GEN_R1_001.fastq
java -jar ~/bin/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads 8 -phred33 $jointreads/GSF1120-144-ctrlM-GEN_R1_001.fastq $trimmed/GSF1120-144-ctrlM-GEN_R1_001.trimmed.fastq ILLUMINACLIP:TruSeq3-SE:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

echo "Samples 139 to 142 - ctrlM GEN"
