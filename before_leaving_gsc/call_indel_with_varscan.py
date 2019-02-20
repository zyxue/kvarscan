#!/usr/bin/env python

import os
import time
import subprocess
import multiprocessing
import argparse
from functools import update_wrapper


def decorator(d):
    "Make function d a decorator: d wraps a function fn."
    def _d(fn):
        return update_wrapper(d(fn), fn)
    update_wrapper(_d, d)
    return _d


@decorator
def timeit(f):
    def new_f(*args, **kwargs):
        bt = time.time()
        r = f(*args, **kwargs)
        et = time.time()        
        print "time spent on {0}: {1:.2f}s".format(f.func_name, et - bt)
        return r
    return new_f


def get_basename_pileup(bam):
    """generate the basename of the pileup file based on the name of the input bam
    file"""
    bam_bn = os.path.basename(bam)
    if bam_bn.endswith('.bam'):
        root, _ = os.path.splitext(bam_bn)
        output_pileup = '{0}.pileup'.format(root)
    else:
        output_pileup = '{0}.pileup'.format(bam_bn)
    return output_pileup


@timeit
def mpileup(input_bam, output_dir):
    output_pileup = os.path.join(output_dir, get_basename_pileup(input_bam))
    output_log = '{0}.log'.format(output_pileup)
    if os.path.exists(output_pileup):
        print '{0} already existed, do nothing'.format(output_pileup)
        return
    cmd = "samtools mpileup -B -f {ref_genome} {input_bam} -o {output_pileup} 2>{output_log}".format(
            ref_genome=REF_GENOME, **locals())
    print cmd
    subprocess.call(cmd, shell=True)


def mpileup_wrap(args):
    mpileup(*args)


def get_somatic_cmd(normal_pileup, tumor_pileup, output_dir, vcf=False):
    # rstrip because varscan will add .indel and .snp to the end of output
    # filenames
    commonprefix = os.path.commonprefix([normal_pileup, tumor_pileup]).rstrip('.')
    output_basename = os.path.join(output_dir, os.path.basename(commonprefix))
    # cmd = "java -jar {varscan_jar} somatic {normal_pileup} {tumor_pileup} {output_basename} --output-vcf".format(**locals())
    cmd = "java -Xmx12g -jar {varscan_jar} somatic {normal_pileup} {tumor_pileup} {output_basename}".format(
        varscan_jar=VARSCAN_JAR, **locals())
    if vcf:
        cmd += ' --output-vcf'
    print cmd
    return cmd
    

@timeit
def call_somatic_mutations_with_varscan(normal_pileup, tumor_pileup, output_dir, vcf=False):
    cmd = get_somatic_cmd(normal_pileup, tumor_pileup, output_dir, vcf)
    subprocess.call(cmd, shell=True)


def call_somatic_mutations_with_varscan_wrapped(args):
    call_somatic_mutations_with_varscan(*args)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref-genome', required=True)
    parser.add_argument('-n', '--normal-bam', required=True)
    parser.add_argument('-t', '--tumor-bam', required=True)
    parser.add_argument('-o', '--output-dir', required=True)
    options = parser.parse_args()

    REF_GENOME = os.path.expanduser(options.ref_genome)
    VARSCAN_JAR = os.path.expanduser('~/Downloads/VarScan.v2.3.7.jar')

    VARSCAN_DIR = os.path.join(options.output_dir, 'varscan')    
    PILEUP_DIR=os.path.join(VARSCAN_DIR, 'pileup')
    SOMATIC_DIR=os.path.join(VARSCAN_DIR, 'somatic')
    for __ in [VARSCAN_DIR, PILEUP_DIR, SOMATIC_DIR]:
        if not os.path.exists(__):
            os.mkdir(__)

    pool = multiprocessing.Pool(2)
    pool.map(mpileup_wrap, [[options.normal_bam, PILEUP_DIR],
                            [options.tumor_bam, PILEUP_DIR]])
    
    normal_pileup = os.path.join(
        PILEUP_DIR, get_basename_pileup(options.normal_bam))
    tumor_pileup = os.path.join(
        PILEUP_DIR, get_basename_pileup(options.tumor_bam))

    pool.map(call_somatic_mutations_with_varscan_wrapped,
             [[normal_pileup, tumor_pileup, SOMATIC_DIR],
              [normal_pileup, tumor_pileup, SOMATIC_DIR, True]])
    # test call:
    # call_indel_with_varscan.py -r ~/ref_data/hg19.fa -n /projects/btl2/zxue/indel_callers/analysis_2/lele.CPCG0100.normal.merged.chr1.175200931.bam -t /projects/btl2/zxue/indel_callers/analysis_2/lele.CPCG0100.tumor.merged.chr1.175200931.bam -o lele_output/
