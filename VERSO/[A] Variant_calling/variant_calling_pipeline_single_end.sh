#

### variant_calling_pipeline_single_end.sh
### Name of file says it all
### Of course, this script assumes an underlying UN*X/Linux platform.

### See the file LICENSE in the main directory for licensing
### information.

###  REQUIRED EXTERNAL TOOLS:
###    trimmomatic
###    bwa
###    samtools
###    picard
###    varscan


### SETTINGS

## Directories.

basedir="/analysis_directory/raw_data/"
resultsdir="/analysis_directory/results/"


## Derived directories.

fastqTrimDir=${resultsdir}trimmed/
bwaSamDir=${resultsdir}bwasam/
bamDir=${resultsdir}bamDir/
vcfDir=${resultsdir}vcfDir/


# Name of the single-end sample to be analyzed.

sampleName="ABCDE"


# Reference genome.
# Reference to be used in the analysis.

genomeFa="analysis_directory/reference/reference.fasta"


# Other setup.

jobs=4


#===B1===

echo "trimmomatic -- trimming fastq"

mkdir -p $fastqTrimDir
if [ ! -d "$fastqTrimDir" ]; then
    echo "Error mkdir $fastqTrimDir"
    exit 1
fi

# We use custom parameters as described in Bastola, Anup, et al. "The
# first 2019 novel coronavirus case in Nepal." The Lancet Infectious
# Diseases 20.3 (2020): 279-280.

# LEADING:20 TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:40

trimmomatic SE -threads $jobs -phred33 \
	    -summary ${fastqTrimDir}${sampleName}_trim.summary \
	    -quiet \
	    ${basedir}${sampleName}.fastq \
	    ${fastqTrimDir}${sampleName}.trim.fastq.gz \
	    LEADING:20 \
	    TRAILING:20 \
	    SLIDINGWINDOW:4:20 \
	    MINLEN:40


#===E1===

#===B2===

echo "bwa mem -- mapping reads to a reference SARS-CoV-2"

mkdir -p $bwaSamDir
if [ ! -d "$bwaSamDir" ]; then
    echo "Error mkdir $bwaSamDir"
    exit 1
fi

bwa mem -t $jobs $genomeFa \
    ${fastqTrimDir}${sampleName}.trim.fastq.gz \
    > ${bwaSamDir}${sampleName}_aln.sam

#===E2===

#===B3===

echo "samtools -- building sorted bam"

mkdir -p $bamDir
if [ ! -d "$bamDir" ]; then
    echo "Error mkdir $bamDir"
    exit 1
fi

samtools view -bT $genomeFa ${bwaSamDir}${sampleName}_aln.sam \
	 > ${bamDir}${sampleName}_aln.bam

samtools sort ${bamDir}${sampleName}_aln.bam \
	 > ${bamDir}${sampleName}_aln.sorted.bam

#===E3===


#===B4===

echo "picard -- removing duplicates"

picard MarkDuplicates \
       I=${bamDir}${sampleName}_aln.sorted.bam \
       O=${bamDir}${sampleName}_aln.sorted_no_duplicates.bam \
       M=${bamDir}${sampleName}_aln.sorted_no_duplicates_metrics.txt \
       REMOVE_DUPLICATES=true

samtools index -b ${bamDir}${sampleName}_aln.sorted_no_duplicates.bam

#===E4===


#===B5===

echo "samtools mpileup -- building mpileup"

samtools mpileup -f $genomeFa \
	 ${bamDir}${sampleName}_aln.sorted_no_duplicates.bam \
	 --output ${bamDir}${sampleName}_aln.sorted_no_duplicates.mpileup

#===E5===

#===B6==

echo "varscan -- calling SNPs"

mkdir -p $vcfDir
if [ ! -d "$vcfDir" ]; then
    echo "Error mkdir $vcfDir"
    exit 1
fi

varscan pileup2snp ${bamDir}${sampleName}_aln.sorted_no_duplicates.mpileup \
	--min-var-freq 0.01 --p-value 1 \
	> ${vcfDir}${sampleName}_aln.sorted_no_duplicates.vcf

#===E6===


### end of file -- variant_calling_pipeline_single_end.sh
