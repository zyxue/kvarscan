#!/usr/bin/env python

import sys

import vcf

from filter_vcf import print_progress, bgzip_and_tabix


def vcf2dict(vcf_reader):
    print 'generate dict from {0}'.format(vcf_reader.filename)
    return {_.POS: _ for _ in vcf_reader}


def vcfs2dicts(*args):
    """convert multiple vcf_readers to dicts"""
    return [vcf2dict(_) for _ in args]


def union(vcf1, vcf2, output_vcf='out.vcf'):
    """vcf1 takes precedence"""

    vcf1_reader = vcf.Reader(filename=vcf1)
    vcf2_reader = vcf.Reader(filename=vcf2)
    d1, d2 = vcfs2dicts(vcf1_reader, vcf2_reader)

    print 'unioning the two'
    d2.update(d1)

    res = d2

    print 'writing output'
    with open(output_vcf, 'wb') as opf:
        vcf_writer = vcf.Writer(opf, vcf1_reader)
        for (k, pos) in enumerate(sorted(d2.keys())):
             vcf_writer.write_record(res[pos])
        print_progress(k, vcf1_reader)

    bgzip_and_tabix(output_vcf)
    print 'DONE'

if __name__ == "__main__":
    union(sys.argv[1], sys.argv[2])
