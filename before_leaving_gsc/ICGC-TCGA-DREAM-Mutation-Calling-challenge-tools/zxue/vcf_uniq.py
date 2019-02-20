import sys
import os

from vcf_union import vcfs2dicts

import vcf

from filter_vcf import print_progress, bgzip_and_tabix


def gen_output_filename(input_vcf, output_vcf):
    if output_vcf is None:
        dirname, basename = os.path.split(input_vcf)
        bn = basename.replace('.vcf.gz', '.uniq.vcf')
        return os.path.join(dirname, bn)
    return output_vcf


def write(rec_dict, uniq_pos, output, vcf_reader):
    print 'writing {0}'.format(output)
    with open(output, 'wb') as opf:
        vcf_writer = vcf.Writer(opf, vcf_reader)
        for k, pos in enumerate(sorted(uniq_pos)):
            vcf_writer.write_record(rec_dict[pos])
        print_progress(k, vcf_reader)
    bgzip_and_tabix(output)


def uniq(vcf1, vcf2, output_vcf1=None, output_vcf2=None):
    """
    compare vcf1 and vcf2 and generates two vcfs that only contain uniq variants
    specific to vcf1 and vcf2, respectively
    """
    vcf1_reader = vcf.Reader(filename=vcf1)
    vcf2_reader = vcf.Reader(filename=vcf2)
    d1, d2 = vcfs2dicts(vcf1_reader, vcf2_reader)

    uniq_to_d1 = sorted(set(d1.keys()) - set(d2.keys()))
    uniq_to_d2 = sorted(set(d2.keys()) - set(d1.keys()))

    output1 = gen_output_filename(vcf1, output_vcf1)
    output2 = gen_output_filename(vcf2, output_vcf2)

    write(d1, uniq_to_d1, output1, vcf1_reader)
    write(d2, uniq_to_d2, output2, vcf2_reader)

if __name__ == "__main__":
    uniq(sys.argv[1], sys.argv[2])

