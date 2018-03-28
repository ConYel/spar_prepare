set -e

# user provided input reads in FASTQ or FASTQ.gz format
FASTQ=$1 

# user provided output directory
OUTDIR=$2

# configuration file
configFile=$3

if [ $# -lt 3 ]; then
  echo "USAGE: $0 <in.fastq> <outputDIR> <configFile>"
  exit 1
fi

# load SPAR configuration file
if [ -s ${configFile} ]; then
  source "${configFile}"
else
  echo "ERROR: ${configFile} not found or empty"
  exit 1
fi

# output bigWig and BAM files
BIGWIGPLUS=${OUTDIR}/raw.pos.bigWig
BIGWIGMINUS=${OUTDIR}/raw.neg.bigWig
BAM=${OUTDIR}/Aligned.out.filtered.hardClipped.sorted.bam
chromInfoFile=${BAM}.chromInfo.txt
SPARPATH=`dirname $0`

# map reads using STAR, filter, sort and output hardclipped sorted BAM
bash ${SPARPATH}/run_star_smrna.sh ${FASTQ} ${maxMismatchCnt} ${maxMapCnt} ${OUTDIR}

# convert BAM to bedgraphs
bash ${SPARPATH}/bam_to_bedgraph2_interval_cut_cpp.sh ${BAM} ${OUTDIR} ${minReadLength} ${maxReadLength} ${keep5pClipped}

# convert bedgraphs to bigWig
bash ${SPARPATH}/bedgraph_to_bigwig.sh ${BAM}.pos.bedgraph ${chromInfoFile} ${BIGWIGPLUS}
bash ${SPARPATH}/bedgraph_to_bigwig.sh ${BAM}.neg.bedgraph ${chromInfoFile} ${BIGWIGMINUS}

echo "bigWig files:"
echo "Positive strand: ${BIGWIGPLUS}"
echo "Negative strand: ${BIGWIGMINUS}"
echo "BAM file:"
echo "${BAM}"
