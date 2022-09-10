
import sys
import re
import os
import pandas as pd
import h5py

_, output_file_name, clues_dir, *clues_file_names = sys.argv # pylint: disable=unbalanced-tuple-unpacking

# open output file:
output_file = open(output_file_name, 'w')

# loop over base names of 
for clues_file_name in clues_file_names:
    chrom, pos, pop = re.search(r'([^_/]+)_(\d+)_([^_]+)$', clues_file_name).groups()
    h5 = h5py.File(os.path.join(clues_dir, clues_file_name) + '.h5', 'r')
    log_likelihood_ratio = h5['logLikelihoodRatios'][h5.attrs['iHat'], h5.attrs['jHat']]
    selection_coef = h5.attrs['sHat']

    print(chrom, pos, pop, log_likelihood_ratio, selection_coef, sep='\t', file=output_file)