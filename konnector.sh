#!/bin/bash

# Function to merge read pairs with abyss-mergepairs
function merge_pairs {
    cell_type=$1
    kon_dir=$2
    read1=$3
    read2=$4

    prefix=${kon_dir}/fafq/${cell_type}_mergepairs

    # expected output file names of abyss-mergepairs
    merge_fifo=${kon_dir}/fafq/${cell_type}_mergepairs_merged.fastq
    read1_fifo=${kon_dir}/fafq/${cell_type}_mergepairs_reads_1.fastq
    read2_fifo=${kon_dir}/fafq/${cell_type}_mergepairs_reads_2.fastq

    rm -f ${merge_fifo} ${read1_fifo} ${read2_fifo}

    mkfifo ${merge_fifo} ${read1_fifo} ${read2_fifo}

    abyss-tofastq --fasta ${merge_fifo} \
	| gzip -c > ${kon_dir}/fafq/${cell_type}_mergepairs_merged.fa.gz &

    gzip_merged_pid=$!
    cat ${read1_fifo} | gzip -c > ${kon_dir}/fafq/${cell_type}_mergepairs_reads_1.fq.gz &
    gzip_reads_1_pid=$!
    cat ${read2_fifo} | gzip -c > ${kon_dir}/fafq/${cell_type}_mergepairs_reads_2.fq.gz &
    gzip_reads_2_pid=$!

    time abyss-mergepairs \
	 -o ${kon_dir}/fafq/${cell_type}_mergepairs \
	 -q 3 $read1 $read2

    wait $gzip_merged_pid $gzip_reads_1_pid $gzip_reads_2_pid

    rm ${merge_fifo} ${read1_fifo} ${read2_fifo}	
}


# Function to join read pairs with konnector
function konnect_pairs {
    cell_type=$1
    kon_dir=$2
    read1=$3
    read2=$4
    k_mer_size=$5
    
    # path of bloom filter
    bloom=${kon_dir}/bf/${cell_type}_k${k_mer_size}.bf.gz
    ls $bloom

    merge_fifo=${kon_dir}/fafq/fifo_${cell_type}_k${k_mer_size}_merged.fa
    read1_fifo=${kon_dir}/fafq/fifo_${cell_type}_k${k_mer_size}_reads_1.fq
    read2_fifo=${kon_dir}/fafq/fifo_${cell_type}_k${k_mer_size}_reads_2.fq

    rm -f ${merge_fifo} ${read1_fifo} ${read2_fifo}

    mkfifo ${merge_fifo} ${read1_fifo} ${read2_fifo}

    cat ${merge_fifo} | gzip -c > ${kon_dir}/fafq/${cell_type}_k${k_mer_size}_merged.fa.gz &
    gzip_merged_pid=$!
    cat ${read1_fifo} | gzip -c > ${kon_dir}/fafq/${cell_type}_k${k_mer_size}_reads_1.fq.gz &
    gzip_reads_1_pid=$!
    cat ${read2_fifo} | gzip -c > ${kon_dir}/fafq/${cell_type}_k${k_mer_size}_reads_2.fq.gz &
    gzip_reads_2_pid=$!

    # for some reason -k and -i must come first --2015-03-24
    time konnector -k ${k_mer_size} -i <(zcat $bloom) \
	 -v -j 16 -P 4 -M nolimit -X 98.0 -B 100 -F 525 -p \
	 >(gzip -c >${kon_dir}/fafq/${cell_type}_k${k_mer_size}_altpaths.fa.gz) \
	 -o "${kon_dir}/fafq/fifo_${cell_type}_k${k_mer_size}" \
	 $read1 \
	 $read2

    wait $gzip_merged_pid $gzip_reads_1_pid $gzip_reads_2_pid
    rm ${merge_fifo} ${read1_fifo} ${read2_fifo}
}

function main {
    cell_type=$1
    kon_dir=$2
    readfile1=$3
    readfile2=$4

    # We run konnector from k90 to k30 for step size of 10, each time
    # propagating unkonnected pairs to the next
    first_k=90 step=-10 last_k=30
    for k_mer_size in $(seq $first_k $step $last_k); do
	stamp=${kon_dir}/fafq/${cell_type}_k${k_mer_size}.COMPLETE
	if [ -e $stamp ]
	then
	    echo "konnector k${k_mer_size} is done previously. Will not re-run ..."
	else
	    if [ ${k_mer_size} -eq $first_k ]
	    then
		konnect_pairs \
		    ${cell_type} ${kon_dir} \
		    ${readfile1} ${readfile2} ${k_mer_size}
	    else
		konnect_pairs \
		    ${cell_type} ${kon_dir} \
		    ${kon_dir}/fafq/${cell_type}_k${prev_k}_reads_1.fq.gz \
		    ${kon_dir}/fafq/${cell_type}_k${prev_k}_reads_2.fq.gz \
		    ${k_mer_size}

		# The unkonnected read files from the previous k are not needed
		# anymore.
		rm -f \
		   ${kon_dir}/fafq/${cell_type}_k${prev_k}_reads_1.fq.gz \
		   ${kon_dir}/fafq/${cell_type}_k${prev_k}_reads_2.fq.gz
	    fi
	    touch $stamp
	fi
	prev_k=${k_mer_size}
    done

    # We run abyss-mergepairs on the unkonnected pairs from k30
    stamp=${kon_dir}/fafq/${cell_type}_mergepairs.COMPLETE
    if [ -e $stamp ]
    then
	echo "abyss-mergepairs is done previously. Will not re-run ..."
    else
	merge_pairs ${cell_type} ${kon_dir} \
	    ${kon_dir}/fafq/${cell_type}_k${last_k}_reads_1.fq.gz \
	    ${kon_dir}/fafq/${cell_type}_k${last_k}_reads_2.fq.gz
	touch $stamp
    fi
}

# TOP_DIR=/projects/btl2/zxue/indel_callers/for_kon_paper/sim_richu
# CELL_TYPES="tumour"

TOP_DIR=/projects/btl2/zxue/indel_callers/for_kon_paper/NA19238
CELL_TYPES="normal"

# CHRS="chr1  chr2  chr3  chr4  chr5  chr6  chr7  chr8  chr9  chr10 \
#       chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 \
#       chr21 chr22 chrX chrY chrM"
CHRS="all"

for cell_type in ${CELL_TYPES}; do
    for chr in $CHRS; do
	chr_dir=${TOP_DIR}/reg/${chr}
	readfile1=${chr_dir}/fafq/${cell_type}_1.fq.gz
	readfile2=${chr_dir}/fafq/${cell_type}_2.fq.gz
	kon_dir=${TOP_DIR}/kon/${chr}
	mkdir -p ${kon_dir}/fafq
	cmd="main ${cell_type} ${kon_dir} ${readfile1} ${readfile2}"
	$cmd &> ${kon_dir}/fafq/${cell_type}_konnector.log
    done
done
