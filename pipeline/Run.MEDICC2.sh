#!/bin/sh
##########################################################
###
### Constuct phylogenetic tree of tumor subclones
### Input: subclone CNV profile
### Output: phylogenetic relationships
### Author: Zhenyu Liu
### 
##########################################################

medicc2=/mnt/d/WSL2/miniconda/miniconda3/bin/medicc2

$medicc2 CNV.subclone.arm.txt \
   ./MEDICC2 \
   --total-copy-numbers \
   --input-allele-columns cn_total