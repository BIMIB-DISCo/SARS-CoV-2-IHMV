#!/bin/bash
#  REQUIRED EXTERNAL TOOLS:
#    trimmomatic 0.39
#    bwa 0.7.17
#    samtools 1.6
#    picard 2.22.1
#    varscan 2.4.4

# SETTINGS
filename="runs.txt"

while read line; do

# performing variant calling (paired-end)
echo $line

basedir="/data/BimiB/share/SARS-CoV-2/raw_data/PRJNA647529/"
genomeFa="/data/BimiB/share/SARS-CoV-2/reference/SARS-CoV-2-ANC.fasta"
resultsdir="/data/BimiB/share/SARS-CoV-2/results/PRJNA647529/"
fastqTrimDir=${resultsdir}trimmed/${line}/
bwaSamDir=${resultsdir}bwasam/${line}/
bamDir=${resultsdir}bamDir/${line}/
vcfDir=${resultsdir}vcfDir/${line}/
coverageDir=${resultsdir}coverage/${line}/
jobs=16
sampleName=$line

#===B1===

echo "trimmomatic -- trimming fastq"

mkdir -p $fastqTrimDir
if [ ! -d "$fastqTrimDir" ]; then
    echo "Error mkdir"
    exit 1
fi

trimmomatic PE -threads $jobs -phred33 -summary ${fastqTrimDir}${sampleName}_trim.summary -quiet -validatePairs ${basedir}${sampleName}_1.fastq.gz ${basedir}${sampleName}_2.fastq.gz ${fastqTrimDir}${sampleName}_1_paired.trim.fastq.gz ${fastqTrimDir}${sampleName}_1_unpaired.trim.fastq.gz ${fastqTrimDir}${sampleName}_2_paired.trim.fastq.gz ${fastqTrimDir}${sampleName}_2_unpaired.trim.fastq.gz LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40

#===E1===

#===B2===

echo "bwa mem -- mapping reads to a reference SARS-CoV-2"

mkdir -p $bwaSamDir
if [ ! -d "$bwaSamDir" ]; then
    echo "Error mkdir"
    exit 1
fi

bwa mem -t $jobs $genomeFa ${fastqTrimDir}${sampleName}_1_paired.trim.fastq.gz ${fastqTrimDir}${sampleName}_2_paired.trim.fastq.gz > ${bwaSamDir}${sampleName}_aln.sam

#===E2===

#===B3===

echo "samtools -- building sorted bam"

mkdir -p $bamDir
if [ ! -d "$bamDir" ]; then
    echo "Error mkdir"
    exit 1
fi

samtools view -bT $genomeFa ${bwaSamDir}${sampleName}_aln.sam > ${bamDir}${sampleName}_aln.bam
samtools sort ${bamDir}${sampleName}_aln.bam > ${bamDir}${sampleName}_aln.sorted.bam

#===E3===

#===B4===

echo "picard -- removing duplicates"

picard MarkDuplicates I=${bamDir}${sampleName}_aln.sorted.bam O=${bamDir}${sampleName}_aln.sorted_no_duplicates.bam M=${bamDir}${sampleName}_aln.sorted_no_duplicates_metrics.txt REMOVE_DUPLICATES=true
samtools index -b ${bamDir}${sampleName}_aln.sorted_no_duplicates.bam

#===E4===

#===B5===

echo "samtools mpileup -- building mpileup"

samtools mpileup -f $genomeFa ${bamDir}${sampleName}_aln.sorted_no_duplicates.bam --output ${bamDir}${sampleName}_aln.sorted_no_duplicates.mpileup

#===E5===

#===B6===

echo "varscan -- calling SNPs"

mkdir -p $vcfDir
if [ ! -d "$vcfDir" ]; then
    echo "Error mkdir"
    exit 1
fi

varscan pileup2snp ${bamDir}${sampleName}_aln.sorted_no_duplicates.mpileup --min-var-freq 0.01 --p-value 1 > ${vcfDir}${sampleName}_aln.sorted_no_duplicates.vcf

#===E6===

#===B7===

mkdir -p $coverageDir
if [ ! -d "$coverageDir" ]; then
    echo "Error mkdir"
    exit 1
fi

echo "samtools depth -- extracting coverage information"

samtools depth -a ${bamDir}${sampleName}_aln.sorted_no_duplicates.bam > ${coverageDir}${sampleName}_aln.sorted_no_duplicates.txt

#===E7===

#===B8===

echo "removing temporary files"

rm -r $fastqTrimDir
rm -r $bwaSamDir
rm ${bamDir}${sampleName}_aln.bam
rm ${bamDir}${sampleName}_aln.sorted.bam
rm ${bamDir}${sampleName}_aln.sorted_no_duplicates_metrics.txt
rm ${bamDir}${sampleName}_aln.sorted_no_duplicates.mpileup

#===E8===

done < $filename
