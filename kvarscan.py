#!/usr/bin/env python
# -*- coding: utf-8 -*


"""
Rule when naming success file:
if there is only one major input file, e.g. input.bam,
then use input.bam.<task_name>.success
"""


from __future__ import unicode_literals

import os
import re
import itertools
import logging
logging.basicConfig(level=logging.INFO)

import ruffus as R

from utils import execute
from argparsers.kvarscan import OPTIONS


REF_GENOME = os.path.expanduser('~/ref_data/hg19.fa')
VARSCAN_JAR = os.path.expanduser('~/jars/VarScan.v2.3.7.jar')
PATH_RE = r'(.*)/(?P<chr>chr[0-9]{1,2}|chr[XYM]|all)/(?P<filetype>\S+)'


input_bams = OPTIONS.bams
print input_bams

@R.transform(input_bams, R.formatter(PATH_RE),
             ['{path[0]}/{basename[0]}.bam.bai',
              '{path[0]}/{basename[0]}.bam.index_bam.SUCCESS'])
def index_bam(input_bam, outputs):
    _, flag = outputs
    cmd = 'samtools index {0}'.format(input_bam)
    execute(cmd, flag=flag)


@R.follows(index_bam)
@R.mkdir(input_bams, R.formatter(PATH_RE), ['{subpath[0][1]}/pileup'])
@R.transform(input_bams, R.formatter(PATH_RE),
             ['{subpath[0][1]}/pileup/{basename[0]}.pileup.gz',
              '{subpath[0][1]}/pileup/{basename[0]}.pileup.SUCCESS',
              '{subpath[0][1]}/pileup/{basename[0]}.pileup_gzip.SUCCESS'])
def pileup_gzip(input_bam, outputs):
    i_bam = input_bam
    o_pileup_gz, flag1, flag2 = outputs
    o_pileup = o_pileup_gz.replace('.gz', '')
    chr_ = re.search(PATH_RE, o_pileup).group('chr')
    if chr_ == 'all':
        cmd1 = ('samtools mpileup -B -f {REF_GENOME} '
                '-o {o_pileup} {i_bam}'.format(REF_GENOME=REF_GENOME, **locals()))
    else:
        cmd1 = ('samtools mpileup -B -f {REF_GENOME} -r {chr_} '
                '-o {o_pileup} {i_bam}'.format(REF_GENOME=REF_GENOME, **locals()))
    execute(cmd1, flag=flag1)
    cmd2 = 'gzip -f {0}'.format(o_pileup)
    execute(cmd2, flag=flag2)


@R.mkdir(pileup_gzip, R.formatter(PATH_RE), ['{subpath[0][1]}/varscan_mpileup2indel'])
@R.transform(pileup_gzip, R.formatter(PATH_RE),
             ['{subpath[0][1]}/varscan_mpileup2indel/{basename[0]}.indel.vcf',
              '{subpath[0][1]}/varscan_mpileup2indel/varscanmpileup2indel.SUCCESS'])
def varscan_mpileup2indel(inputs, outputs):
    """by default: append --output-vcf"""
    pileup_gz, _, _ = inputs
    output, flag1 = outputs
    cmd1 = ('java -Xmx12g -jar {VARSCAN_JAR} mpileup2indel '
            '<(zcat {pileup_gz}) --output-vcf '
            '> {output}'.format(VARSCAN_JAR=VARSCAN_JAR, **locals()))
    execute(cmd1, flag=flag1)


@R.mkdir(pileup_gzip, R.formatter(PATH_RE), ['{subpath[0][1]}/varscan_mpileup2indel'])
@R.transform(pileup_gzip, R.formatter(PATH_RE),
             ['{subpath[0][1]}/varscan_mpileup2indel/{basename[0]}.indel.csv',
              '{subpath[0][1]}/varscan_mpileup2indel/varscanmpileup2indel.csv.SUCCESS'])
def varscan_mpileup2indel_csv(inputs, outputs):
    pileup_gz, _, _ = inputs
    output, flag1 = outputs
    cmd1 = ('java -Xmx12g -jar {VARSCAN_JAR} mpileup2indel '
            '<(zcat {pileup_gz}) > {output}'.format(VARSCAN_JAR=VARSCAN_JAR, **locals()))
    execute(cmd1, flag=flag1)


@R.transform(pileup_gzip, R.formatter(PATH_RE),
             ['{subpath[0][1]}/varscan_mpileup2indel/{basename[0]}.pset1.indel.csv',
              '{subpath[0][1]}/varscan_mpileup2indel/varscanmpileup2indel.pset1.csv.SUCCESS'])
def varscan_mpileup2indel_pset1_csv(inputs, outputs):
    pileup_gz, _, _ = inputs
    output, flag1 = outputs
    cmd1 = ('java -Xmx12g -jar {VARSCAN_JAR} mpileup2indel '
            '<(zcat {pileup_gz}) '
            '--min-var-freq 0 '
            '>{output}'.format(VARSCAN_JAR=VARSCAN_JAR, **locals()))
    execute(cmd1, flag=flag1)


@R.transform(pileup_gzip, R.formatter(PATH_RE),
             ['{subpath[0][1]}/varscan_mpileup2indel/{basename[0]}.pset2.indel.csv',
              '{subpath[0][1]}/varscan_mpileup2indel/varscanmpileup2indel.pset2.csv.SUCCESS'])
def varscan_mpileup2indel_pset2_csv(inputs, outputs):
    pileup_gz, _, _ = inputs
    output, flag1 = outputs
    cmd1 = ('java -Xmx12g -jar {VARSCAN_JAR} mpileup2indel '
            '<(zcat {pileup_gz}) '
            '--min-var-freq 0 '
            '--p-value 1 '
            '>{output}'.format(VARSCAN_JAR=VARSCAN_JAR, **locals()))
    execute(cmd1, flag=flag1)


@R.transform(pileup_gzip, R.formatter(PATH_RE),
             ['{subpath[0][1]}/varscan_mpileup2indel/{basename[0]}.pset3.indel.csv',
              '{subpath[0][1]}/varscan_mpileup2indel/varscanmpileup2indel.pset3.csv.SUCCESS'])
def varscan_mpileup2indel_pset3_csv(inputs, outputs):
    pileup_gz, _, _ = inputs
    output, flag1 = outputs
    cmd1 = ('java -Xmx12g -jar {VARSCAN_JAR} mpileup2indel '
            '<(zcat {pileup_gz}) '
            '--min-var-freq 0 '
            '--p-value 1 '
            '--strand-filter 0 '
            '>{output}'.format(VARSCAN_JAR=VARSCAN_JAR, **locals()))
    execute(cmd1, flag=flag1)


@R.transform(pileup_gzip, R.formatter(PATH_RE),
             ['{subpath[0][1]}/varscan_mpileup2indel/{basename[0]}.pset4.indel.csv',
              '{subpath[0][1]}/varscan_mpileup2indel/varscanmpileup2indel.pset4.csv.SUCCESS'])
def varscan_mpileup2indel_pset4_csv(inputs, outputs):
    pileup_gz, _, _ = inputs
    output, flag1 = outputs
    cmd1 = ('java -Xmx12g -jar {VARSCAN_JAR} mpileup2indel '
            '<(zcat {pileup_gz}) '
            '--min-var-freq 0 '
            '--p-value 1 '
            '--strand-filter 0 ' # default 1
            '--min-reads2 1 '    # default 2
            '>{output}'.format(VARSCAN_JAR=VARSCAN_JAR, **locals()))
    execute(cmd1, flag=flag1)



@R.mkdir(pileup_gzip, R.formatter(PATH_RE), ['{subpath[0][1]}/varscan_somatic'])
@R.collate(pileup_gzip, R.formatter(PATH_RE),
           ['{subpath[0][1]}/varscan_somatic/pset4.snp.vcf.gz',
            '{subpath[0][1]}/varscan_somatic/pset4.snp.vcf.gz.tbi',
            '{subpath[0][1]}/varscan_somatic/pset4.indel.vcf.gz',
            '{subpath[0][1]}/varscan_somatic/pset4.indel.vcf.gz.tbi',
            '{subpath[0][1]}/varscan_somatic/varscan_somatic_pset4.SUCCESS',
            '{subpath[0][1]}/varscan_somatic/varscan_somatic_pset4_bgzip_tabix_snp.SUCCESS',
            '{subpath[0][1]}/varscan_somatic/varscan_somatic_pset4_bgzip_tabix_indel.SUCCESS'])
def varscan_somatic_pset4(inputs, outputs):
    inputs = [_ for _ in itertools.chain(*inputs)
              if not _.endswith('.SUCCESS')]
    snp_vcf, snp_tbi, indel_vcf, indel_tbi, flag1, flag2, flag3 = outputs
    output_prefix = snp_vcf.replace('.snp.vcf.gz', '')
    # how to preserve the order of normal, tumour?, assert is just a dirty
    # quick hack
    assert len(inputs) == 2
    assert 'normal' in inputs[0]
    assert 'tumour' in inputs[1]
    cmd1 = ('java -Xmx12g -jar {VARSCAN_JAR} somatic '
            '<(zcat {normal_pileup_gz}) '
            '<(zcat {tumour_pileup_gz}) '
            '{output_prefix} '
            '--min-var-freq 0 '
            '--p-value 1 '
            '--strand-filter 0 ' # default 1
            '--min-reads2 1 '    # default 2
            '--somatic-p-value 1 '
            '--output-vcf'.format(VARSCAN_JAR=VARSCAN_JAR,
                                  normal_pileup_gz=inputs[0],
                                  tumour_pileup_gz=inputs[1],
                                  output_prefix=output_prefix))
    execute(cmd1, flag=flag1)
    cmd2 = 'bgzip -f {0} && tabix -p vcf {1}'.format(
        snp_vcf.replace('.gz', ''), snp_vcf)
    execute(cmd2, flag=flag2)
    cmd3 = 'bgzip -f {0} && tabix -p vcf {1}'.format(
        indel_vcf.replace('.gz', ''), indel_vcf)
    execute(cmd3, flag=flag3)



@R.mkdir(pileup_gzip, R.formatter(PATH_RE), ['{subpath[0][1]}/varscan_somatic'])
@R.collate(pileup_gzip, R.formatter(PATH_RE),
           ['{subpath[0][1]}/varscan_somatic/base6.snp.vcf.gz',
            '{subpath[0][1]}/varscan_somatic/base6.snp.vcf.gz.tbi',
            '{subpath[0][1]}/varscan_somatic/base6.indel.vcf.gz',
            '{subpath[0][1]}/varscan_somatic/base6.indel.vcf.gz.tbi',
            '{subpath[0][1]}/varscan_somatic/varscan_somatic_base6.SUCCESS',
            '{subpath[0][1]}/varscan_somatic/varscan_somatic_base6_bgzip_tabix_snp.SUCCESS',
            '{subpath[0][1]}/varscan_somatic/varscan_somatic_base6_bgzip_tabix_indel.SUCCESS'])
def varscan_somatic_base6(inputs, outputs):
    inputs = [_ for _ in itertools.chain(*inputs)
              if not _.endswith('.SUCCESS')]
    snp_vcf, snp_tbi, indel_vcf, indel_tbi, flag1, flag2, flag3 = outputs
    output_prefix = snp_vcf.replace('.snp.vcf.gz', '')
    # how to preserve the order of normal, tumour?, assert is just a dirty
    # quick hack
    assert len(inputs) == 2
    assert 'normal' in inputs[0]
    assert 'tumour' in inputs[1]
    cmd1 = ('java -Xmx12g -jar {VARSCAN_JAR} somatic '
            '<(zcat {normal_pileup_gz}) '
            '<(zcat {tumour_pileup_gz}) '
            '{output_prefix} '
            '--min-var-freq 0 '
            '--min-coverage 1 '
            '--min-coverage-normal 0 '
            '--min-coverage-tumor 1 '
            '--somatic-p-value 1 '
            '--min-reads2 1 '
            '--min-avg-qual 0 '
            '--output-vcf'.format(VARSCAN_JAR=VARSCAN_JAR,
                                  normal_pileup_gz=inputs[0],
                                  tumour_pileup_gz=inputs[1],
                                  output_prefix=output_prefix))
    execute(cmd1, flag=flag1)
    cmd2 = 'bgzip -f {0} && tabix -p vcf {1}'.format(
        snp_vcf.replace('.gz', ''), snp_vcf)
    execute(cmd2, flag=flag2)
    cmd3 = 'bgzip -f {0} && tabix -p vcf {1}'.format(
        indel_vcf.replace('.gz', ''), indel_vcf)
    execute(cmd3, flag=flag3)


if __name__ == "__main__":
    for _ in R.pipeline_get_task_names ():
        print _
    print os.environ['PWD']
    print '=' * 79

    logger, logger_mutex = R.cmdline.setup_logging(
        __name__, OPTIONS.log_file, OPTIONS.verbose)

    logger.info(OPTIONS.verbose)
    # with logger_mutex:
    #     logger.info("Look Ma. No hands")
    R.pipeline_printout_graph('lele.svg', draw_vertically=True)
    R.cmdline.run(OPTIONS)
