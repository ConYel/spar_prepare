# SPAR config file

export HOMEDIR="${HOME}"

#absolute path to the bin directory
export BINDIR="${HOMEDIR}/bin"

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
export minReadLength=14 # minimum read length
export max5pClip=1 # maximum allowed 5' clip; all read with 5' clipping > max will be discarded 
export keep5pClipped=0
