# SPAR: Small RNA-seq Portal for Analysis of sequencing expeRiments

This goal of this section is to provide guidelines on how to prepare raw sequencing files (n FASTQ format) for SPAR analysis.
Raw sequencing reads are first trimmed to remove adapters used in library preparation/sequencing (smrna_adapter_cut.sh), then mapped to the genome (e.g., hg19) (run_star_smrna.sh) to obtain aligned reads in BAM format. BAM file is then converted to bigWig raw signal track files.

For adapter trimming, use
```
bash smrna_adapter_cut.sh <raw.fastq.gz>
```
to obtain trimmed.fastq.gz 

To prepare bigWig/BAM files for use with [SPAR webserver](https://www.lisanwanglab.org/SPAR)
```
prepare_BAM_and_bigWigs_from_fastq.sh <trimmed.fastq.gz> <output-directory> <config-file>
```
to obtain bigWig and BAM files.

Detailed instructions on

a. how to prepare raw FASTQ file
b. how to map reads to the genome and generate a BAM file
c. how to generate BigWig files

are provided below.

## Preparing reference genome (FASTA):

```
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit

wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/twoBitToFa

chmod a+x twoBitToFa

./twoBitToFa hg19.2bit hg19.fa
```

## Preparing STAR index for the reference genome

```
mkdir -p hg19/star

STAR --runMode genomeGenerate --genomeDir hg19/star --genomeFastaFiles hg19.fa --runThreadN 4
```

## Setting up configuration file

Example configuration file is given in provided *config.hg19.sh* file:
```
#SPAR config file

export HOMEDIR="${HOME}"

#absolute path to the bin directory
export BINDIR="${HOMEDIR}/bin"

#reference genome
export GENOMEBUILD=hg19

#absolute path to the STAR genome index
export STAR="${BINDIR}/STAR_2.4.0j/bin/Linux_x86_64/STAR" # STAR
export genomeDir="${HOMEDIR}/datasets/${GENOMEBUILD}/star/"  # STAR genome index
export GENOMEFA="${HOMEDIR}/datasets/${GENOMEBUILD}/${GENOMEBUILD}.fa"

#absolute path to pre-installed STAR, samtools, AWK, etc
export SAMTOOLS="${BINDIR}/samtools-1.2/samtools"
export BEDTOOLS="${BINDIR}/bedtools-2.26/bedtools2/bin/bedtools"

# UCSC tools
export BGTOBIGWIG="${BINDIR}/bedGraphToBigWig"
export BEDTOBIGBED="${BINDIR}/bedToBigBed"

# Adapter trimming
export CUTADAPT="${BINDIR}/cutadapt-1.8.1/bin/cutadapt"

#mapping parameters for STAR
export maxMismatchCnt=0 # maximum number of genomic mismatches
export maxMapCnt=100 # maximum number of places a read can map to
export minMappedLength=14 # minimum *mapped* length
export maxReadLength=44 # maximum read length
export max5pClip=1 # maximum allowed 5' clip; all read with 5' clipping > max will be discarded
export keep5pClipped=0 # by default all reads clipped at 5' are excluded from analysis
```
Please set locations of the programs/tools in the config file to the locations in your system.

### Download/install if necessary any missing programs/tools:

STAR
```
wget https://github.com/alexdobin/STAR/archive/2.5.4b.tar.gz
```

UCSC tools
```
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedGraphToBigWig
```

Samtools
```
apt-get install samtools
```


## Recompiling bam2bedgraph binary (if necessary)

### Linux:
Install htslib
```
wget https://github.com/samtools/htslib/releases/download/1.7/htslib-1.7.tar.bz2
tar xjvf htslib-1.7.tar.bz2
cd htslib-1.7
make
make install
```

Compile bam2bedgraph:
```
g++ -O2 bam2bedgraph.cpp -o bam2bedgraph_interval -lhts
```

### MAC OS:

```
brew install htslib

g++ -O2 bam2bedgraph.cpp -o bam2bedgraph_interval -lhts
```

## Trimming reads in FASTQ.gz

To trim small RNA-seq adapters (small RNA adapters v1.0, v1.5, TruSeq):
```
bash smrna_adapter_cut.sh raw.fastq.gz
```
This will produce trimmed files in
```
raw.fastq.trimmed.gz.
```

## Preparing bigWig and BAM files from trimmed FASTQ

To prepare bigWig files and BAM from the example trimmed and gzipped FASTQ file:
```
bash prepare_BAM_and_bigWigs_from_fastq.sh example-data.fastq.gz test_out config.hg19.sh
```
Output files will be saved into *test_out* directory:
```
Mar 28 16:52:18 ..... Mapping reads
Mar 28 16:52:18 ..... Started STAR run
Mar 28 16:52:18 ..... Started mapping

...

bigWig files:
Positive strand: test_out/raw.pos.bigWig
Negative strand: test_out/raw.neg.bigWig

BAM file:
test_out/Aligned.out.filtered.hardClipped.sorted.bam
```
Expected outputs are included in the *example_out* directory.
