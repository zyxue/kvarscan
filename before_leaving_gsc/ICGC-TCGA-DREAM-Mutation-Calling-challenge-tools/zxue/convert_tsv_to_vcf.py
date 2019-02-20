#!/usr/bin/env python

import argparse

import pandas as pd
import vcf

from filter_vcf import print_progress, bgzip_and_tabix

def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--tsv', required=True)
    parser.add_argument('-t', '--template-vcf', required=True,
                        help=('the template vcf, must contain the corresponding '
                              'records in tsv, not how to generate Record in PyVCF'))
    parser.add_argument('-o', '--output-vcf', required=True)
    return parser




if __name__ == "__main__":
    parser = get_parser()
    options = parser.parse_args()

    # true positive tuples
    tp_tuples = pd.read_csv(options.tsv, sep='\t',
                         index_col=['position', 'ref', 'var']
                     ).index.to_native_types()
        
    vcf_reader = vcf.Reader(filename=options.template_vcf)
    output_vcf = options.output_vcf
    output_records_count = 0
    with open(output_vcf, 'wb') as opf:
        vcf_writer = vcf.Writer(opf, vcf_reader)
        for k, rec in enumerate(vcf_reader):
            sig = (rec.POS, rec.REF, ','.join(map(str, rec.ALT))) # like a signature
            if sig in tp_tuples:
                vcf_writer.write_record(rec)
                output_records_count += 1
            print_progress(k, vcf_reader)
        print '{0} records processed'.format(k + 1, vcf_reader.filename)
        print '{0} records to write'.format(output_records_count)
    bgzip_and_tabix(output_vcf)
