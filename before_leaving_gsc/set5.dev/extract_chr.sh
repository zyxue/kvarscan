#!/bin/bash

TOP_DIR=.
CELL_TYPES="normal"
# CHRS="chr1  chr2  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chr10 \
#       chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
#       chr21 chr22 chrX chrY chrM"

CHRS=chr1

function extract {
    top_dir=$1
    cell_type=$2
    chr=$3
    chr_dir=${top_dir}/reg/${chr}

    mkdir -p ${chr_dir}/fafq
    mkdir -p ${chr_dir}/bam

    paren_bam=${top_dir}/reg/all/bam/${cell_type}.bam
    # must exist
    ls ${paren_bam}
    child_bam=${chr_dir}/bam/${cell_type}.bam
    
    echo "${cell_type}, extract ${chr} bam"
    samtools view -hb -f 2 -F 2048 -o ${child_bam} ${paren_bam} ${chr}

    echo "${cell_type}, indexing"
    samtools index ${chr_dir}/bam/${cell_type}.bam

    echo "${cell_type}, sort by name"
    samtools sort -n ${chr_dir}/bam/${cell_type}.bam ${chr_dir}/bam/${cell_type}.sorted_by_name

    echo "${cell_type}, bedtools bamtofastq"
    bedtools bamtofastq \
	-i   ${chr_dir}/bam/${cell_type}.sorted_by_name.bam \
	-fq  ${chr_dir}/fafq/${cell_type}_1.fq \
	-fq2 ${chr_dir}/fafq/${cell_type}_2.fq

    echo "${cell_type}, gziping fq"
    gzip ${chr_dir}/fafq/${cell_type}_1.fq
    gzip ${chr_dir}/fafq/${cell_type}_2.fq

    rm -v ${chr_dir}/bam/${cell_type}.sorted_by_name.bam
}

for cell_type in ${CELL_TYPES}; do
    for chr in ${CHRS}; do
	cmd="extract $TOP_DIR $cell_type $chr"
	echo $cmd
	$cmd
    done
done
