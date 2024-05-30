#!/bin/sh

# USAGE : perfom HiC processing with hicstuff.


# ---- VARIABLES ----
WORKPATH=/Users/emiliehhd/Projects
READSPATH=./HiC/data/reads/
RESULTSPATH=./HiC/results_4k/

PREFIX=JG38
ENZYME=HpaII,HinfI
GENOME=./0_ref_genome/ecoli_mg1655/ecoli_k12_mg1655.fa
#4000
BINSIZE=4000
THREADS=8


# ---- PROCESSING ----
cd ${WORKPATH}
mkdir ${RESULTSPATH}${PREFIX}

hicstuff pipeline --aligner bowtie2 \
    --threads ${THREADS} --enzyme ${ENZYME} \
    --binning ${BINSIZE} \
    --mapping iterative --duplicates --circular --filter --plot --force --no-cleanup \
    --genome ${GENOME} ${READSPATH}${PREFIX}_nxq_R1.fq.gz ${READSPATH}${PREFIX}_nxq_R2.fq.gz \
    --outdir ${RESULTSPATH}${PREFIX} --prefix ${PREFIX}

