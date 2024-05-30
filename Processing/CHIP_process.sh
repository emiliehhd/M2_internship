#!/bin/sh

# USAGE : perfom ChIP-seq processing with an input.


# ---- VARIABLES ----
READSPATH=/Users/emiliehhd/Projects/CHIPseq/11032024_JG_CHIPseq_NDD/data/fastq/
RESULTSPATH=/Users/emiliehhd/Projects/CHIPseq/11032024_JG_CHIPseq_NDD/results_mg/

PREFIX_IP=JG71
PREFIX_input=JG72

STRAIN=mg1655
GENOME=/Users/emiliehhd/Projects/CHIPseq/11032024_JG_CHIPseq_NDD/data/ref/mg1655/ecoli_k12_mg1655
GENOME_SIZE=4641652
#mg1655: 4641652
#lf82 (genome+plasmid) : 4881487
THREADS=8


# ---- IP ALIGNMENT ----
echo IP Alignment
mkdir -p ${RESULTSPATH}fastq/${PREFIX_IP}
mkdir -p ${RESULTSPATH}bam/${PREFIX_IP}

TMP_FQ_PATH_IP=${RESULTSPATH}fastq/${PREFIX_IP}
BAM_PATH_IP=${RESULTSPATH}bam/${PREFIX_IP}

bowtie2 --maxins 1000 --threads ${THREADS} \
    -x ${GENOME} -1 ${READSPATH}${PREFIX_IP}_nxq_R1.fq.gz -2 ${READSPATH}${PREFIX_IP}_nxq_R2.fq.gz \
    --maxins 1000 \
    --un-conc-gz ${TMP_FQ_PATH_IP}/${PREFIX_IP}^unmapped_${STRAIN}.gz \
    > ${BAM_PATH_IP}/${PREFIX_IP}^mapped_${STRAIN}.sam


# ---- INPUT ALIGNMENT ----
echo Input Alignment
mkdir -p ${RESULTSPATH}fastq/${PREFIX_input}
mkdir -p ${RESULTSPATH}bam/${PREFIX_input}

TMP_FQ_PATH_input=${RESULTSPATH}fastq/${PREFIX_input}
BAM_PATH_input=${RESULTSPATH}bam/${PREFIX_input}

bowtie2 --maxins 1000 --threads ${THREADS} \
    -x ${GENOME} -1 ${READSPATH}${PREFIX_input}_nxq_R1.fq.gz -2 ${READSPATH}${PREFIX_input}_nxq_R2.fq.gz \
    --maxins 1000 \
    --un-conc-gz ${TMP_FQ_PATH_input}/${PREFIX_input}^unmapped_${STRAIN}.gz \
    > ${BAM_PATH_input}/${PREFIX_input}^mapped_${STRAIN}.sam


# ---- INPUT SAMTOOLS ----
echo IP samtools
echo ${BAM_PATH_IP}/${PREFIX_IP}^mapped_${STRAIN}.sam 
samtools fixmate -@ ${THREADS} \
    --output-fmt bam \
    -m ${BAM_PATH_IP}/${PREFIX_IP}^mapped_${STRAIN}.sam - | \
    samtools sort -@ 8 --output-fmt bam -T ${BAM_PATH_IP}/${PREFIX_IP}^mapped_${STRAIN}.sam_sorting - | \
    samtools markdup -@ 8 --output-fmt bam -r - - | \
    samtools view -@ 8 --output-fmt bam -f 0x001 -f 0x002 -F 0x004 -F 0x008 -q 10 -1 -b - | \
    samtools sort -@ 8 --output-fmt bam -l 9 -T ${BAM_PATH_IP}/${PREFIX_IP}^mapped_${STRAIN}.sam_sorting2 \
    -o ${BAM_PATH_IP}/${PREFIX_IP}^mapped_${STRAIN}.bam

samtools index -@ ${THREADS} ${BAM_PATH_IP}/${PREFIX_IP}^mapped_${STRAIN}.bam


# ---- INPUT SAMTOOLS ----
echo Input samtools
samtools fixmate -@ ${THREADS} \
    --output-fmt bam \
    -m ${BAM_PATH_input}/${PREFIX_input}^mapped_${STRAIN}.sam - | \
    samtools sort -@ 8 --output-fmt bam -T ${BAM_PATH_input}/${PREFIX_input}^mapped_${STRAIN}.sam_sorting - | \
    samtools markdup -@ 8 --output-fmt bam -r - - | \
    samtools view -@ 8 --output-fmt bam -f 0x001 -f 0x002 -F 0x004 -F 0x008 -q 10 -1 -b - | \
    samtools sort -@ 8 --output-fmt bam -l 9 -T ${BAM_PATH_input}/${PREFIX_input}^mapped_${STRAIN}.sam_sorting2 \
    -o ${BAM_PATH_input}/${PREFIX_input}^mapped_${STRAIN}.bam

samtools index -@ ${THREADS} ${BAM_PATH_input}/${PREFIX_input}^mapped_${STRAIN}.bam


# ---- CREATING IP BIGWIG ----
echo IP bigwig
mkdir -p ${RESULTSPATH}tracks/${PREFIX_IP}
TRACKS_PATH_IP=${RESULTSPATH}tracks/${PREFIX_IP}

bamCoverage --bam ${BAM_PATH_IP}/${PREFIX_IP}^mapped_${STRAIN}.bam \
    --outFileName  ${TRACKS_PATH_IP}/${PREFIX_IP}^mapped_${STRAIN}.CPM.bw \
    --binSize 1 --numberOfProcessors ${THREADS} --normalizeUsing CPM --skipNonCoveredRegions --extendReads


# ---- CREATING INPUT BIGWIG ----
echo input bigwig
mkdir -p ${RESULTSPATH}tracks/${PREFIX_input}
TRACKS_PATH_INPUT=${RESULTSPATH}tracks/${PREFIX_input}

bamCoverage --bam ${BAM_PATH_input}/${PREFIX_input}^mapped_${STRAIN}.bam \
    --outFileName  ${TRACKS_PATH_INPUT}/${PREFIX_input}^mapped_${STRAIN}.CPM.bw \
    --binSize 1 --numberOfProcessors ${THREADS} --normalizeUsing CPM --skipNonCoveredRegions --extendReads


# ---- COMPARING IP/INPUT  ----
echo Comparing IP and input
bamCompare -b1 ${BAM_PATH_IP}/${PREFIX_IP}^mapped_${STRAIN}.bam \
    -b2 ${BAM_PATH_input}/${PREFIX_input}^mapped_${STRAIN}.bam \
    --outFileName ${TRACKS_PATH_IP}/${PREFIX_IP}^mapped_${STRAIN}.vs-${PREFIX_INPUT}.bw \
    --scaleFactorsMethod readCount --operation ratio --skipZeroOverZero --skipNAs --numberOfProcessors ${THREADS} \
    --binSize 5 --skipNonCoveredRegions --extendReads


# ---- CALLPEAK  ----
echo Callpeak
mkdir -p ${RESULTSPATH}peaks/${PREFIX_IP}
macs2 callpeak -t ${BAM_PATH_IP}/${PREFIX_IP}^mapped_${STRAIN}.bam \
    -c ${BAM_PATH_input}/${PREFIX_input}^mapped_${STRAIN}.bam \
    --format BAMPE --gsize ${GENOME_SIZE} \
    --outdir ${RESULTSPATH}peaks/${PREFIX_IP} \
    --name ${PREFIX_IP}_vs-${PREFIX_input}_genome-${STRAIN}


rm -rf ${RESULTSPATH}fastq/${PREFIX_IP}
rm -rf ${RESULTSPATH}fastq/${PREFIX_input}