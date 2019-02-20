#!/usr/bin/env python

import sys
import csv

import vcf

header = [
    'chrom',
    'position',
    'ref',
    'var',
    'normal_reads1',
    'normal_reads2',
    'normal_var_freq',
    'normal_gt',
    'tumor_reads1',
    'tumor_reads2',
    'tumor_var_freq',
    'tumor_gt',
    'somatic_status',
    'variant_p_value',
    'somatic_p_value',
    'tumor_reads1_plus',
    'tumor_reads1_minus',
    'tumor_reads2_plus',
    'tumor_reads2_minus',
    'normal_reads1_plus',
    'normal_reads1_minus',
    'normal_reads2_plus',
    'normal_reads2_minus']


if len(sys.argv) != 2:
    print 'usage: convert_vcf_to_tsv.py some.vcf.gz'
    sys.exit(1)

input_vcf = sys.argv[1]
vcf_fh = vcf.Reader(filename=input_vcf)
output = input_vcf.rstrip('.vcf.gz') + '.tsv'
with open(output, 'wb') as opf:
    writer = csv.DictWriter(opf, fieldnames=header, delimiter='\t')
    writer.writeheader()
    for k, rec in enumerate(vcf_fh):
        if k > 0 and (k + 1) % 1e4 == 0:
            print '{0} rows processed'.format(k + 1)
        try:
            row = [rec.CHROM, rec.POS, rec.REF, ','.join(map(str, rec.ALT)),
                   rec.samples[0]['RD'],
                   rec.samples[0]['AD'],
                   rec.samples[0]['FREQ'],
                   rec.samples[0]['GT'],
                   rec.samples[1]['RD'],
                   rec.samples[1]['AD'],
                   rec.samples[1]['FREQ'],
                   rec.samples[1]['GT'],
                   ['Germline', 'Somatic', 'LOH', 'Unkonwn'][int(rec.INFO['SS']) - 1],
                   rec.INFO['GPV'],
                   rec.INFO['SPV']]
            row += rec.samples[1]['DP4'].split(',')
            row += rec.samples[0]['DP4'].split(',')
            writer.writerow(dict(zip(header, row)))
        except AttributeError as err:
            print '{0} at Row {1}'.format(err, k+1)
            sys.exit(1)
    if (k + 1) % 1e4 == 0:
        print '{0} rows processed'.format(k + 1)
