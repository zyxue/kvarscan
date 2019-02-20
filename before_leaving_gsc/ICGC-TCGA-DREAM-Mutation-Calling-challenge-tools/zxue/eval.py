#!/usr/bin/env python

"""An more informative and modualized version of evaluator.py from
https://github.com/Sage-Bionetworks/ICGC-TCGA-DREAM-Mutation-Calling-challenge-tools/blob/master/evaluator.py"""


import sys
import os
import subprocess

import vcf

'''
Submission evaluation code for TCGA/ICGC/DREAM SMC
Adam Ewing, ewingad@soe.ucsc.edu
Requires PyVCF (https://github.com/jamescasbon/PyVCF)
'''

def snv_or_indel_match(r1, r2):
    for __ in ['POS', 'REF', 'ALT']:
        if getattr(r1, __) != getattr(r2, __):
            return False
    return True


def is_true_positive(subrec, trurec, vtype='SNV'):
    if vtype == 'SNV' and subrec.is_snp and trurec.is_snp:
        if snv_or_indel_match(subrec, trurec):
            return True

    if vtype == 'INDEL' and subrec.is_indel and trurec.is_indel:
        if snv_or_indel_match(subrec, trurec):
            return True

    if vtype == 'SV' and subrec.is_sv and trurec.is_sv:
        trustart, truend = expand_sv_ends(trurec)
        substart, subend = expand_sv_ends(subrec)

        # check for overlap
        if min(truend, subend) - max(trustart, substart) > 0:
            return True

    return False


def expand_sv_ends(rec):
    ''' assign start and end positions to SV calls using conf. intervals if present '''
    startpos, endpos = rec.start, rec.end
    assert rec.is_sv

    try:
        if rec.INFO.get('END'): # sometimes this is a list, sometimes it's an int
            if isinstance(rec.INFO.get('END'), list):
                endpos = int(rec.INFO.get('END')[0])
            if isinstance(rec.INFO.get('END'), int):
                endpos = int(rec.INFO.get('END'))

        if rec.INFO.get('CIPOS'):
            ci = map(int, rec.INFO.get('CIPOS'))
            if ci[0] < 0:
                startpos += ci[0]

        if rec.INFO.get('CIEND'):
            ci = map(int, rec.INFO.get('CIEND')) 
            if ci[0] > 0:
                endpos += ci[0]

    except TypeError as e:
        sys.stderr.write("error expanding sv interval: " + str(e) + " for record: " + str(rec) + "\n")

    if startpos > endpos:
        endpos, startpos = startpos, endpos

    return startpos, endpos


def countrecs(submission, truth, vtype='SNV', ignorechroms=None, truthmask=True):
    ''' return number of records in submission '''
    subvcfh = vcf.Reader(filename=submission)
    truvcfh = vcf.Reader(filename=truth)

    truchroms = dict([(trurec.CHROM, True) for trurec in truvcfh])
    subrecs = 0

    for subrec in subvcfh:
        if pass_filter(subrec):
            if (ignorechroms is None or subrec.CHROM not in ignorechroms):
                if not mask(subrec, truvcfh, truchroms, active=truthmask):
                    if subrec.is_snp and vtype == 'SNV':
                        subrecs += 1
                    if subrec.is_sv and vtype == 'SV':
                        subrecs += 1
                    if subrec.is_indel and vtype == 'INDEL':
                        subrecs += 1
    return subrecs


def is_relevant(rec):
    # 'ignore' types are always excluded
    if rec.INFO.get('SVTYPE'):
        if rec.INFO.get('SVTYPE') in ('IGN', 'MSK'):
            return False 
    return True


def is_of_right_vtype(rec, vtype):
    if ((vtype == 'SNV' and rec.is_snp)
        or (vtype == 'SV' and rec.is_sv)
        or (vtype == 'INDEL' and rec.is_indel)):
        return True
    return False


def pass_filter(rec):
    ''' Return true if a record is unfiltered or has 'PASS' in the filter field
    (pyvcf sets FILTER to None) '''
    if (rec.FILTER is None) or (rec.FILTER == '.') or (not rec.FILTER):
        return True
    return False


def evaluate(submission, truth, vtype='SNV'):
    ''' return stats on sensitivity, specificity, balanced accuracy '''
    sub_vcf_fh = vcf.Reader(filename=submission)
    true_vcf_fh = vcf.Reader(filename=truth)

    true_vcf_list = []
    for k, rec in enumerate(true_vcf_fh):
        if is_relevant(rec) and is_of_right_vtype(rec, vtype):
            true_vcf_list.append(rec)
    n_true_recs = len(true_vcf_list)
    print ("from {file}:\n\t{num_true}/{total} are of type {vtype} "
           "and relevant, that's {percentage:.2%}".format(
               file=true_vcf_fh.filename, num_true=n_true_recs,
               total=k + 1, vtype=vtype, 
               percentage=float(n_true_recs) / (k + 1)))

    tp, fp, fn, unrelevant = [], [], true_vcf_list, []
    used_truth = []
    for k, rec in enumerate(sub_vcf_fh):
        flag_tp = False
        if (is_relevant(rec)
            and is_of_right_vtype(rec, vtype)
            and pass_filter(rec)):

            try:
                for ref_rec in true_vcf_fh.fetch(rec.CHROM, rec.start, rec.end):
                    if is_true_positive(rec, ref_rec, vtype=vtype):
                        used_truth.append(str(ref_rec))
                        # convenient for convert_vcf_to_tsv.py since it has all
                        # the info
                        tp.append(rec)
                        fn.remove(ref_rec)
                        flag_tp = True
                        break
            except ValueError:
                # e.g. could not create iterator for region 'hs37d5:34917771-34917772'
                pass

            if not flag_tp:
                fp.append(rec)
        else:
            unrelevant.append(rec)

        if k > 0 and (k + 1) % 5e3 == 0:
            print '{0} records processed in {1}'.format(k + 1, sub_vcf_fh.filename)
    n_sub_recs  = k + 1
    print '{0} records processed in {1}'.format(k + 1, sub_vcf_fh.filename)

    tp_count, fp_count = len(tp), len(fp)
    sensitivity = float(tp_count) / float(n_true_recs)
    precision   = float(tp_count) / float(tp_count + fp_count)
    specificity = 1.0 - float(fp_count) / float(n_sub_recs)

    print precision
    print specificity

    print "sensitivity: {0}".format(sensitivity)
    print "specificity: {0}".format(specificity)
    balanced_accuracy = (sensitivity + specificity) / 2.0
    print "balanced accuracy: {0}".format(balanced_accuracy)

    print 'writing vcf files for debugging...'
    write_files(tp, fp, fn, submission, submission.rstrip('indel.vcf.gz'))


def write_files(tp, fp, fn, template, outprefix):
    write_vcf(tp, template, '{0}.true_positive.indel.vcf'.format(outprefix))
    write_vcf(fp, template, '{0}.false_positive.indel.vcf'.format(outprefix))
    write_vcf(fn, template, '{0}.false_negative.indel.vcf'.format(outprefix))


def write_vcf(records, template, outfile):
    reader = vcf.Reader(filename=template)
    with open(outfile, 'wb') as opf:
        print 'writing {0}'.format(outfile)
        writer = vcf.Writer(opf, reader)
        for __ in records:
            writer.write_record(__)
    bgzip_and_tabix(outfile)


def bgzip_and_tabix(vcf_file):
    subprocess.call('bgzip -f {0} && tabix -f -p vcf {1}'.format(
        vcf_file, vcf_file + '.gz'), shell=True)



if __name__ == '__main__':
    if len(sys.argv) != 4 and len(sys.argv) != 5:
        print ("standalone usage for testing:",
               sys.argv[0],
               ("<submission VCF> <truth VCF (tabix-indexed)> "
                "<SV, SNV, or INDEL> "
                "[ignore chrom list (comma-delimited, optional)]"))
        sys.exit(1)

    subvcf, truvcf, evtype = sys.argv[1:4]

    if not subvcf.endswith('.vcf') and not subvcf.endswith('.vcf.gz'):
        sys.stderr.write("submission VCF filename does not enc in .vcf or .vcf.gz\n")
        sys.exit(1)

    if not os.path.exists(truvcf + '.tbi'):
        sys.stderr.write("truth VCF does not appear to be indexed. bgzip + tabix index required.\n")
        sys.exit(1)

    if evtype not in ('SV', 'SNV', 'INDEL'):
        sys.stderr.write("last arg must be either SV, SNV, or INDEL\n")
        sys.exit(1)

    print "unmasked:"
    evaluate(subvcf, truvcf, vtype=evtype)
    # count  = countrecs(subvcf, truvcf, vtype=evtype, truthmask=False)
    # print "number of unmasked mutations in submission: " + str(count)
