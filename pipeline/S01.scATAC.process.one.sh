#!/bin/sh

cell=$1
ref=$2
rawdata_dir=$3

ref_dir=/gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/database/${ref}/hisat2

mkdir -p ${cell}
cd ${cell}

# trim ME sequence
mkdir -p trim_fq

cutadapt \
    -m 25 \
    -j 4 \
    -u -1 \
    -U -1 \
    -a CTGTCTCTTATACACATCT \
    -A CTGTCTCTTATACACATCT \
    -o trim_fq/${cell}_1_trimmed.fq.gz \
    -p trim_fq/${cell}_2_trimmed.fq.gz \
    ${rawdata_dir}/${cell}/${cell}_1.fq.gz \
    ${rawdata_dir}/${cell}/${cell}_2.fq.gz


# hisat2 mapping
mkdir mapping

hisat2 \
    -X 2000 \
    -p 4 \
    --no-spliced-alignment \
    -x ${ref_dir}/genome \
    -1 trim_fq/${cell}_1_trimmed.fq.gz \
    -2 trim_fq/${cell}_2_trimmed.fq.gz \
    --summary-file mapping/${cell}_aln_sum.txt | \
    samtools view -ShuF 4 -f 2 -q 30 - | \
    samtools sort - -o mapping/${cell}_f2q30.bam

# picard mark duplication
java -jar -Xmx4g \
    /gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/picard_2.19.0/picard.jar \
    MarkDuplicates \
    INPUT=mapping/${cell}_f2q30.bam \
    OUTPUT=${cell}_f2q30_pmd.bam \
    REMOVE_DUPLICATES=true \
    ASSUME_SORTED=true \
    METRICS_FILE=${cell}_f2q30_pmd.out

cd ..
mv scATAC*${cell}*err scATAC*${cell}*out log
