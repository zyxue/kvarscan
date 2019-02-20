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

import yaml
# expected to be run in the directory where there is a kvarscan.yaml
with open('./kvarscan.yaml') as inf:
    CONFIG = yaml.load(inf.read())

import ruffus as R


from utils import execute
from argparsers.konnector import OPTIONS

REF_GENOME = os.path.expanduser('~/ref_data/hg19.fa')
VARSCAN_JAR = os.path.expanduser('~/jars/VarScan.v2.3.7.jar')
PATH_RE = (r'(?P<prefix>.*)/'
            '(?P<protocal>reg|kon)/'
            '(?P<chr>chr[0-9]{1,2}|chr[XYM]|all)/'
            '(?P<filetype>\S+)/'
            '(?P<celltype>\S+)_[12]\.fq\.gz')


INPUT_FQS = [OPTIONS.fq1, OPTIONS.fq2]


def gen_vars(input_fqs):
    """generate variables for the task of abyss_bloom, konnector"""
    k_sizes = range(*CONFIG['abyss_bloom']['k_mer_sizes'])
    sr = re.search(PATH_RE, input_fqs[0])
    sr2 = re.search(PATH_RE, input_fqs[1])
    # should be of conventional directory hierarchy
    try:
        assert sr.groups() == sr2.groups()
    except AssertionError:
        print '{0} != {1}'.format(sr.groups(), sr2.groups())
        raise

    bfs, bf_flags, fas, fa_flags = [], [], [], []
    for k_size in k_sizes:
        # for abyss_bloom
        # bn: basename
        bf_bn = '{0}_k{1}.bf.gz'.format(sr.group('celltype'), k_size)
        bf_flag_bn = '{0}.SUCCESS'.format(bf_bn)
        bf_dir = os.path.join(sr.group('prefix'), 'kon', sr.group('chr'), 'bf')
        bf = os.path.join(bf_dir, bf_bn)
        bf_flag = os.path.join(bf_dir, bf_flag_bn)
        bfs.append(bf)
        bf_flags.append(bf_flag)

        # for konnector
        fa_all_bn = '{0}_k{1}_allpaths.fa.gz'.format(sr.group('celltype'), k_size)
        fa_mer_bn = '{0}_k{1}_merged.fa.gz'.format(sr.group('celltype'), k_size)
        fa_flag_bn = '{0}_k{1}.SUCCESS'.format(sr.group('celltype'), k_size)
        fa_dir = os.path.join(sr.group('prefix'), 'kon', sr.group('chr'), 'fafq')
        fa_all = os.path.join(fa_dir, fa_all_bn)
        fa_mer = os.path.join(fa_dir, fa_mer_bn)
        fa_flag = os.path.join(fa_dir, fa_flag_bn)
        fas.extend([fa_all, fa_mer])
        fa_flags.append(fa_flag)

    return k_sizes, bfs, bf_flags, fas, fa_flags

K_MER_SIZES, BFS, BF_FLAGS, FAS, FA_FLAGS = gen_vars(INPUT_FQS)

for __ in FAS:
    print __

for __ in FA_FLAGS:
    print __


@R.mkdir(INPUT_FQS, R.formatter(PATH_RE), ['{prefix[0]}/kon/{chr[0]}/bf'])
@R.collate(INPUT_FQS, R.formatter(), BFS + BF_FLAGS)
def abyss_bloom(input_fqs, outputs):
    fq1, fq2 = input_fqs
    for k_mer_size, bf, bf_flag in zip(K_MER_SIZES, BFS, BF_FLAGS):
        cmd = CONFIG['abyss_bloom']['cmd'].format(**locals())
        # cmd = ('abyss-bloom build -v -k {k_mer_size} -j 8 -b 3G -l 2 -q 15 - '
        #        '{fq1} {fq2} '
        #        '| gzip -c > {bf}'.format(**locals()))
        execute(cmd, flag=bf_flag)


@R.follow(abyss_bloom)
@R.mkdir(abyss_bloom, R.formatter(PATH_RE), ['{subpath[0][1]}/fafq'])
@R.collate(abyss_bloom, R.formatter(), FAS + FA_FLAGS)
def konnector(input_fqs, outputs):
    for k_mer_size, fa, fa_flag in zip(K_MER_SIZES, FAS, FA_FLAGS):


#              ['{subpath[0][1]}/pileup/{basename[0]}.pileup.gz',
#               '{subpath[0][1]}/pileup/{basename[0]}.pileup.SUCCESS',
#               '{subpath[0][1]}/pileup/{basename[0]}.pileup_gzip.SUCCESS'])
# def pileup_gzip(input_bam, outputs):
#     i_bam = input_bam
#     o_pileup_gz, flag1, flag2 = outputs
#     o_pileup = o_pileup_gz.replace('.gz', '')
#     chr_ = re.search(PATH_RE, o_pileup).group('chr')
#     if chr_ == 'all':
#         cmd1 = ('samtools mpileup -B -f {REF_GENOME} '
#                 '-o {o_pileup} {i_bam}'.format(REF_GENOME=REF_GENOME, **locals()))
#     else:
#         cmd1 = ('samtools mpileup -B -f {REF_GENOME} -r {chr_} '
#                 '-o {o_pileup} {i_bam}'.format(REF_GENOME=REF_GENOME, **locals()))
#     execute(cmd1, flag=flag1)
#     cmd2 = 'gzip -f {0}'.format(o_pileup)
#     execute(cmd2, flag=flag2)


# @R.mkdir(pileup_gzip, R.formatter(PATH_RE), ['{subpath[0][1]}/varscan_mpileup2indel'])
# @R.transform(pileup_gzip, R.formatter(PATH_RE),
#              ['{subpath[0][1]}/varscan_mpileup2indel/{basename[0]}.indel.csv',
#               '{subpath[0][1]}/varscan_mpileup2indel/varscanmpileup2indel.SUCCESS'])
# def varscan_mpileup2indel_csv(inputs, outputs):
#     pileup_gz, _, _ = inputs
#     output, flag1 = outputs
#     cmd1 = ('java -Xmx12g -jar {VARSCAN_JAR} mpileup2indel '
#             '<(zcat {pileup_gz}) > {output}'.format(VARSCAN_JAR=VARSCAN_JAR, **locals()))
#     execute(cmd1, flag=flag1)




# @R.mkdir(pileup_gzip, R.formatter(PATH_RE), ['{subpath[0][1]}/varscan_somatic'])
# @R.collate(pileup_gzip, R.formatter(PATH_RE),
#            ['{subpath[0][1]}/varscan_somatic/base6.snp.vcf.gz',
#             '{subpath[0][1]}/varscan_somatic/base6.snp.vcf.gz.tbi',
#             '{subpath[0][1]}/varscan_somatic/base6.indel.vcf.gz',
#             '{subpath[0][1]}/varscan_somatic/base6.indel.vcf.gz.tbi',
#             '{subpath[0][1]}/varscan_somatic/varscan_somatic_base6.SUCCESS',
#             '{subpath[0][1]}/varscan_somatic/varscan_somatic_base6_bgzip_tabix_snp.SUCCESS',
#             '{subpath[0][1]}/varscan_somatic/varscan_somatic_base6_bgzip_tabix_indel.SUCCESS'])
# def varscan_somatic_base6(inputs, outputs):
#     inputs = [_ for _ in itertools.chain(*inputs)
#               if not _.endswith('.SUCCESS')]
#     snp_vcf, snp_tbi, indel_vcf, indel_tbi, flag1, flag2, flag3 = outputs
#     output_prefix = snp_vcf.replace('.snp.vcf.gz', '')
#     # how to preserve the order of normal, tumour?, assert is just a dirty
#     # quick hack
#     assert len(inputs) == 2
#     assert 'normal' in inputs[0]
#     assert 'tumour' in inputs[1]
#     cmd1 = ('java -Xmx12g -jar {VARSCAN_JAR} somatic '
#             '<(zcat {normal_pileup_gz}) '
#             '<(zcat {tumour_pileup_gz}) '
#             '{output_prefix} '
#             '--min-var-freq 0 '
#             '--min-coverage 1 '
#             '--min-coverage-normal 0 '
#             '--min-coverage-tumor 1 '
#             '--somatic-p-value 1 '
#             '--min-reads2 1 '
#             '--min-avg-qual 0 '
#             '--output-vcf'.format(VARSCAN_JAR=VARSCAN_JAR,
#                                   normal_pileup_gz=inputs[0],
#                                   tumour_pileup_gz=inputs[1],
#                                   output_prefix=output_prefix))
#     execute(cmd1, flag=flag1)
#     cmd2 = 'bgzip -f {0} && tabix -p vcf {1}'.format(
#         snp_vcf.replace('.gz', ''), snp_vcf)
#     execute(cmd2, flag=flag2)
#     cmd3 = 'bgzip -f {0} && tabix -p vcf {1}'.format(
#         indel_vcf.replace('.gz', ''), indel_vcf)
#     execute(cmd3, flag=flag3)


if __name__ == "__main__":
    for _ in R.pipeline_get_task_names ():
        print _
    print os.environ['PWD']
    print '=' * 79

    # parser = R.cmdline.get_argparse(
    #     description="krvarscan",
    #     usage='require python-2.7.x',
    #     version='0.1')
    # options = parser.parse_args()
    logger, logger_mutex = R.cmdline.setup_logging(
        __name__, OPTIONS.log_file, OPTIONS.verbose)

    logger.info(OPTIONS.verbose)
    # with logger_mutex:
    #     logger.info("Look Ma. No hands")
    R.pipeline_printout_graph('lele.svg', draw_vertically=True)
    R.cmdline.run(OPTIONS)
