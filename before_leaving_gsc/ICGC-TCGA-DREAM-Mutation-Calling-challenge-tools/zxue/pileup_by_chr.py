import os
import argparse
import Queue
from threading import Thread
import subprocess


def gen_cmds(chrs, types, ref_genome, dir_bam, dir_pileup):
    for chr_str in chrs:
        for type_ in types:
            in_bam = os.path.join(dir_bam, '{type}.bam'.format(type=type_))
            out_pileup = os.path.join(
                dir_pileup, '{type}.{chr_str}.pileup'.format(type=type_, chr_str=chr_str))
            pileup_log = out_pileup + '.log'
            cmd = construct_mpileup_cmd(
                ref_genome, in_bam, out_pileup, chr_str, pileup_log)
            yield cmd


def construct_mpileup_cmd(ref_fa, in_bam, out_pileup, chr_str, pileup_log=None):
    cmd = 'samtools mpileup -B -f {ref_fa} {in_bam} -r {chr_str}  -o {out_pileup}'.format(**locals())
    if pileup_log is not None:
        cmd += ' 2>{pileup_log}'.format(**locals())
    return cmd


def worker(q):
    while True:
        cmd = q.get()
        print cmd
        subprocess.call(cmd, shell=True)
        q.task_done()


def get_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', '--ref-genome', required=True,
                        choices=['hg19.fa', 'hs37d5.fa'],
                        help='it will look for the corresponding file in /home/zxue/ref_data')
    parser.add_argument('-t', '--types', nargs='+', required=True,
                        help='a type of normal indicates the existence of ./bam/normal.bam')
    return parser


def main():
    options = get_parser().parse_args()

    dir_ref_data = '/home/zxue/ref_data'
    # synthetic challenge set3, set4
    if options.ref_genome == 'hs37d5.fa':
        chrs = range(1, 23) + ['X', 'Y', 'MT'];
    # set5
    elif options.ref_genome == 'hg19.fa':
        chrs = ['chr{0}'.format(__) for __ in range(1, 23) + ['X', 'Y', 'M']]
    else:
        raise ValueError('reference genome is not known: {0}'.format(options.ref_genome))

        
    ref_genome = os.path.expanduser(os.path.join(dir_ref_data, options.ref_genome))

    dir_top = os.getenv('PWD')
    dir_bam = os.path.join(dir_top, 'bam')
    dir_varscan = os.path.join(dir_top, 'varscan')
    dir_pileup=os.path.join(dir_varscan, 'pileup')
    dir_somatic=os.path.join(dir_varscan, 'somatic')

    for __ in [dir_varscan, dir_pileup, dir_somatic]:
        if not os.path.exists(__):
            os.mkdir(__)

    num_worker_threads = 8
    q = Queue.Queue()
    for i in range(num_worker_threads):
         t = Thread(target=worker, args=(q,))
         t.daemon = True
         t.start()

    for item in gen_cmds(chrs, options.types, ref_genome, dir_bam, dir_pileup):
        q.put(item)

    q.join()   

if __name__ == "__main__":
    main()
