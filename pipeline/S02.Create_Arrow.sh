#!/bin/sh

##########################################################
###
### Data processing for scATAC-seq #2
### Input: sorted bam file
### Output: fragments file and arrow file (for ArchR)
### Author: Zhenyu Liu
### Email:  liuzhenyu@pku.edu.cn
### 
##########################################################

ref=$1
sample=$2

mkdir 02.arrow
cd 02.arrow

#==========================================================================
#      Step1 Merge single cell bam files
#==========================================================================

ls ../01.process/*/*.bam > file.list
samtools merge -@ 5 -b file.list -r ${sample}.merge.bam

samtools index ${sample}.merge.bam

#==========================================================================
#      Step2 convert to fragment files, index
#==========================================================================
sinto fragments -b ${sample}.merge.bam \
    -f ${sample}.fragments \
    -t RG \
    --use_chrom "(?i)^([1-9]|X)" \
    -p 5

sort -k 1,1 -k2,2n ${sample}.fragments > ${sample}.sort.fragments

rm -f ${sample}.fragments

bgzip ${sample}.sort.fragments

tabix -p bed ${sample}.sort.fragments.gz

#==========================================================================
#      Step3 Create arrow files for ArchR
#==========================================================================

Rscript /gpfs1/tangfuchou_pkuhpc/tangfuchou_coe/liuzhenyu/software/bin/Create_Arrow.R $ref ${sample}.sort.fragments.gz
