#!/usr/bin/env python

import os
import sys
import argparse
import subprocess

import vcf


def get_dummy_parser():
    dp = argparse.ArgumentParser(add_help=False) # dp: dummy_parser
    dp.add_argument('-f', '--input-vcf', required=True)
    dp.add_argument('-o', '--output-vcf', default=None)
    return dp


def get_parser():
    dummy_parser = get_dummy_parser()
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(metavar='')
    sp_gen_somatic = subparsers.add_parser(
        'gen-somatic-vcf', parents=[dummy_parser],
        help='Generate a vcf file with only somatic variants.')
    sp_gen_somatic.set_defaults(func=main_gen_somatic_vcf)

    sp_filter_vcf = subparsers.add_parser(
        'filter-vcf', parents=[dummy_parser],
        help='Generate a vcf file with only somatic variants.')
    sp_filter_vcf.set_defaults(func=main_filter_vcf)

    sp_filter_vcf.add_argument(
        '-s', '--spv', dest='spv', type=float,
        default=1, required=True,
        help="cutoff for somatic p-value")
    return parser


def get_somatic_vcf_output_name(input_vcf):
    return input_vcf.\
        replace('.indel.vcf', '.somatic.indel.vcf').\
        strip('.gz')


def is_somatic(record):
    return 'SOMATIC' in record.INFO


def print_progress(k, reader=None):
    if (k > 0) and ((k + 1) % 1e4 == 0):
        if reader:
            print '{0} records processed from {1}'.format(k + 1, reader.filename)
        else:
            print '{0} records processed'.format(k + 1)


def generate_somatic_vcf(output_vcf, vcf_reader):
    with open(output_vcf, 'wb') as opf:
        vcf_writer = vcf.Writer(opf, vcf_reader)
        num_total, num_selected = 0, 0        
        k = 0 # if no records exist in vcf_reader, there will be an exception
              # at the end print, so better to initialize it
        for k, record in enumerate(vcf_reader):                
            if is_somatic(record):
                vcf_writer.write_record(record)
                num_selected += 1
            num_total += 1
            print_progress(k, vcf_reader)
    bgzip_and_tabix(output_vcf)
    print '{0} records processed'.format(k + 1, vcf_reader.filename)
    print "total: {0}, # selected, {1}, that's {2:.2%}".format(
        num_total, num_selected, num_selected / float(num_total))


def main_gen_somatic_vcf(options):
    vcf_reader = vcf.Reader(filename=options.input_vcf)
    if options.output_vcf is None:
        output_vcf = get_somatic_vcf_output_name(options.input_vcf)
    else:
        output_vcf = options.output_vcf

    if os.path.exists(output_vcf):
        print '{0} already exists'.format(output_vcf)
        sys.exit(1)

    generate_somatic_vcf(output_vcf, vcf_reader)
    

def transform_spv(num):
    return 'spv_{0}'.format('{0:f}'.format(num).replace('.', '_').rstrip('0'))


def get_filter_vcf_output_name(input_vcf, options):
    spv = transform_spv(options.spv)
    return input_vcf.\
        replace('.indel.vcf', '.{spv}.indel.vcf'.format(spv=spv)).\
        strip('.gz')


def main_filter_vcf(options):
    vcf_reader = vcf.Reader(filename=options.input_vcf)
    if options.output_vcf is None:
        output_vcf = get_filter_vcf_output_name(options.input_vcf, options)
    else:
        output_vcf = options.output_vcf
    filter_vcf(output_vcf, vcf_reader, options.spv)


def filter_vcf(output_vcf, vcf_reader, spv):
    with open(output_vcf, 'wb') as opf:
        vcf_writer = vcf.Writer(opf, vcf_reader)
        num_total, num_selected = 0, 0
        for k, record in enumerate(vcf_reader):
            if is_somatic(record) and record.INFO['SPV'] < spv:
                vcf_writer.write_record(record)
                num_selected += 1
            num_total += 1
            print_progress(k, vcf_reader)
    bgzip_and_tabix(output_vcf)
    print '{0} records processed from {1}'.format(k + 1, vcf_reader.filename)
    print "total: {0}, # selected, {1}, that's {2:.2%}".format(
        num_total, num_selected, num_selected / float(num_total))

def bgzip_and_tabix(vcf_file):
    subprocess.call('bgzip -f {0} && tabix -f -p vcf {1}'.format(
        vcf_file, vcf_file + '.gz'), shell=True)
    
if __name__ == "__main__":
    options = get_parser().parse_args()
    # _, ext = os.path.splitext(options.input_vcf)
    if not options.input_vcf.endswith('.indel.vcf.gz'):
        raise ValueError(
            "input vcf file doesn't end with .indel.vcf.gz".format(options.input_vcf))
    options.func(options)
