# RNA PCR Primer, Index 1 (RPI1)
# 5’ CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA
#
# grep -e "TGGAATTCTCGGG" test.pos.sam --color
# grep -e "ACCCGAGAATTCCA" test.neg.sam --color

# RNA 3' Adapter (RA3)
# 5’ TGGAATTCTCGGGTGCCAAGG

FASTQ=$1 # compressed fastq.gz
THREADS=${2:-2} # change threads

if [ $# -lt 1 ]; then
  echo "USAGE: $0 <in.fastq.gz>"
  exit 1
fi

if [[ "${FASTQ}" != *.fastq.gz ]]; then
  echo "ERROR: compressed *.fastq.gz is required."
  exit 1;
fi


if [ ! -s "${FASTQ}" ]; then
  echo "ERROR: ${FASTQ} not found";
  exit 1
fi

#>illumina_smrna_1.0_3p
smrna_1_0=TCGTATGCCGTCTTCTGCTTG
#>illumina_smrna_1.5_3p(Illumina_Small_RNA_3p_Adapter_1)
smrna_1_5=ATCTCGTATGCCGTCTTCTGCTTG
#>illumina_truseq
smrna_truseq=TGGAATTCTCGGGTGCCAAGG
smrna_ribo=AGATCGGAAGAGCACACGTCT


#adapterOptionA="-a CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTC -a AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC -a ACACTCTTTCCCTACACGACGCTCTTCCGATCT -a TCTCGGCATTCCTGCTGAACCGCTCTTCCGATCTTTTTTTTTTTTVN -a AAAAAAAAAAAA -a TTTTTTTTTTTT"
#adapterOptionB="-g CAAGCAGAAGACGGCATACGAGATCGGTCTCGGCATTCCTGCTGAACCGCTCTTC -g AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATC -g ACACTCTTTCCCTACACGACGCTCTTCCGATCT -g TCTCGGCATTCCTGCTGAACCGCTCTTCCGATCTTTTTTTTTTTTVN -g AAAAAAAAAAAA -g TTTTTTTTTTTT"

#RA3="TGGAATTCTCGGGTGCCAAGG"
#PCR="CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
#PCR="CAAGCAGAAGACGGCATACGAGATCGTGATGTGACTGGAGTTCCTTGGCACCCGAGAATTCCA"
#PCRrc=$( echo "${PCR}" | rev | tr ATCG TAGC )

CUTADAPT=${CUTADAPT:-$(command -v cutadapt)}
if [ ! -x "${CUTADAPT}" ]; then
  echo "ERROR: cutadapt is not found!"
  exit 1
fi

log=${FASTQ}.cutadapt_log
>${log}
trimmed=${FASTQ/.fastq.gz/.trimmed.fastq.gz} # compressed output
untrimmed=${FASTQ/.fastq/.untrimmed.fastq}
untrimmed=/dev/null # discard untrimmed reads
adapterOptionA="-a ${smrna_1_0} -a ${smrna_1_5} -a ${smrna_truseq} -a ${smrna_ribo}"

# trimming parameters
MIN_TRIMMED_READ_LEN=14
MAX_ADAPTER_ERROR=0.06
MIN_ADAPTER_MATCH=6

#${CUTADAPT} ${adapterOptionA} -o ${trimmed}  --untrimmed-output ${untrimmed} -O 6 -m 14 -n 8 "${FASTQ}"  | tee -a ${log}
# removed --untrimmed-output option because cutadapt can't run on multiple threads
#Multi-core Cutadapt can only write to output files given by -o and -p.
#This implies that the following command-line arguments are not compatible with multi-core: --untrimmed-output 

${CUTADAPT} ${adapterOptionA} -o ${trimmed} -O ${MIN_ADAPTER_MATCH} -m ${MIN_TRIMMED_READ_LEN} -e ${MAX_ADAPTER_ERROR} --cores=${THREADS} -n 8 "${FASTQ}"  | tee -a ${log}

echo "Trimmed reads: ${trimmed}"
echo "Untrimmed reads: ${untrimmed}"
echo "Log: ${log}"
