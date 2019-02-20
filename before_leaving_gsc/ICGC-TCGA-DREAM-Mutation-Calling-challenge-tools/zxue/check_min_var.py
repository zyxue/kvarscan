#!/usr/bin/env python

import sys

import pandas as pd

def p2f(x):
    # percentage string to float
    return float(x.strip('%')) / 100.




infile = sys.argv[1]            # e.g. ./varscan/somatic/chr19.indel
print infile
df = pd.read_csv(infile, sep='\t')
df = df[df.somatic_status == 'Somatic']

df['tumor_var_freq'] = df.tumor_var_freq.apply(p2f)
df['normal_var_freq'] = df.normal_var_freq.apply(p2f)

print df[['tumor_var_freq', 'normal_reads2', 'tumor_reads2', 'somatic_p_value']].describe()
