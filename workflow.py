from gwf import Workflow, AnonymousTarget
import os, re, sys
import numpy as np
from subprocess import PIPE, Popen

gwf = Workflow()

def modpath(p, parent=None, base=None, suffix=None):
    par, name = os.path.split(p)
    name_no_suffix, suf = os.path.splitext(name)
    if type(suffix) is str:
        suf = suffix
    if parent is not None:
        par = parent
    if base is not None:
        name_no_suffix = base

    new_path = os.path.join(par, name_no_suffix + suf)
    if type(suffix) is tuple:
        assert len(suffix) == 2
        new_path, nsubs = re.subn(r'{}$'.format(suffix[0]), suffix[1], new_path)
        assert nsubs == 1, nsubs
    return new_path


def extract_window(with_indels, chrom, win_start, win_end, pop):
    output_vcf_file = f'steps/recode_vcf/{chrom}_{win_start}_{win_end}/{chrom}_{win_start}_{win_end}_{pop}.recode.vcf'
    vcf_base_name = modpath(output_vcf_file, suffix=('.recode.vcf', ''))
    
    inputs = [with_indels]
    outputs = [output_vcf_file]
    options = {
        'memory': '8g',
        'walltime': '05:00:00'
    }
    
    spec = f'''
    
    mkdir -p steps/recode_vcf/{chrom}_{win_start}_{win_end}

    vcftools --gzvcf {with_indels} --chr {chrom} --from-bp {win_start} --to-bp {win_end-1} \
        --keep ~/simons/faststorage/data/1000Genomes/metainfo/{pop}_male.txt \
        --keep ~/simons/faststorage/data/1000Genomes/metainfo/{pop}_female.txt \
        --remove-indels --remove-filtered-all --max-alleles 2 --recode --non-ref-ac-any 1 \
        --out {vcf_base_name}
    
    '''
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def vcf2sites(chrom, win_start, win_end, pop):
    vcffile=f'steps/recode_vcf/{chrom}_{win_start}_{win_end}/{chrom}_{win_start}_{win_end}_{pop}.recode.vcf'
    sitesfile=f'steps/sitesfiles/{chrom}_{win_start}_{win_end}/{chrom}_{win_start}_{win_end}_{pop}.sites'
    
    inputs = [vcffile]
    outputs = [sitesfile]
    options = {
        'memory': '8g',
        'walltime': '01:00:00'
    }
    
    spec = f'''
    
    mkdir -p steps/sitesfiles/{chrom}_{win_start}_{win_end}

    python scripts/VCFtoSITES.py {vcffile} {sitesfile}
    
    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def argsample(chrom, win_start, win_end, pop):
    sitesfile=f'steps/sitesfiles/{chrom}_{win_start}_{win_end}/{chrom}_{win_start}_{win_end}_{pop}.sites'
    times_file=f'data/tennessen_times_fine.txt'
    popsize_file=f'data/tennessen_popsize_fine.txt'
    arg_sample_file=f'steps/argweaver/argsample/argsample_{chrom}_{win_start}_{win_end}_{pop}/argsample_{chrom}_{win_start}_{win_end}_{pop}'
    argweaver_bedfile=f'steps/argweaver/argsample_{chrom}_{win_start}_{win_end}_{pop}/argsample_{chrom}_{win_start}_{win_end}_{pop}.bed.gz'

    inputs = [sitesfile]
    outputs = [argweaver_bedfile]
    options = {
        'memory': '40g',
        'walltime': '14-00:00:00'
    }

    spec = f'''

    mkdir -p steps/argweaver/argsample/argsample_{chrom}_{win_start}_{win_end}_{pop}

    mkdir -p steps/argweaver/argsample_{chrom}_{win_start}_{win_end}_{pop}
    
    arg-sample -s {sitesfile} --times-file {times_file} --popsize-file {popsize_file} -r 1e-8 -m 1.2e-8 --ntimes 20 --maxtime 200e3 \
        -c 25 -n 30000 --resume -o {arg_sample_file}
    
    ../../software/argweaver/bin/smc2bed-all {arg_sample_file}
    
    cp -p steps/argweaver/argsample/argsample_{chrom}_{win_start}_{win_end}_{pop}/argsample_{chrom}_{win_start}_{win_end}_{pop}.bed.gz {argweaver_bedfile}
    
    cp -p steps/argweaver/argsample/argsample_{chrom}_{win_start}_{win_end}_{pop}/argsample_{chrom}_{win_start}_{win_end}_{pop}.bed.gz.tbi \
        steps/argweaver/argsample_{chrom}_{win_start}_{win_end}_{pop}/argsample_{chrom}_{win_start}_{win_end}_{pop}.bed.gz.tbi
    
    cp -p steps/argweaver/argsample/argsample_{chrom}_{win_start}_{win_end}_{pop}/argsample_{chrom}_{win_start}_{win_end}_{pop}.log \
        steps/argweaver/argsample_{chrom}_{win_start}_{win_end}_{pop}/argsample_{chrom}_{win_start}_{win_end}_{pop}.log
    
    '''
    
    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)


def clues(snp_position, chrom, win_start, win_end, pop, derived_allele, derived_freq):
    bedfile=f'steps/argweaver/argsample_{chrom}_{win_start}_{win_end}_{pop}/argsample_{chrom}_{win_start}_{win_end}_{pop}.bed.gz'
    logfile=f'steps/argweaver/argsample_{chrom}_{win_start}_{win_end}_{pop}/argsample_{chrom}_{win_start}_{win_end}_{pop}.log'
    sitesfile=f'steps/sitesfiles/{chrom}_{win_start}_{win_end}/{chrom}_{win_start}_{win_end}_{pop}.sites'
    treesfile=f'$JOBDIR/{chrom}_{win_start}_{win_end}_{pop}.trees'
    cond_trans_matrix_file='data/trans_tennesen_fine.hdf5'
    clues_output_file=f'steps/clues/clues_{chrom}_{win_start}_{win_end}_{pop}/clues_{chrom}_{snp_position}_{pop}.h5'
    clues_output_base_name=modpath(clues_output_file, suffix=('.h5', ''))
    
    inputs = [bedfile]
    outputs = [clues_output_file]
    options = {
        'memory': '8g',
        'walltime': '1-00:00:00'
    }

    spec = f'''
    
    JOBDIR=/scratch/$GWF_JOBID

    mkdir -p steps/clues/clues_{chrom}_{win_start}_{win_end}_{pop}
    
    arg-summarize -a {bedfile} -r {chrom}:{snp_position}-{snp_position} -l {logfile} -E > {treesfile}
    
    python ../../software/clues/clues.py {treesfile} {cond_trans_matrix_file} {sitesfile} {derived_freq} --posn {snp_position} \
        --derivedAllele {derived_allele} --noAncientHap --approx 10000 --thin 10 --burnin 100 --output {clues_output_base_name}

    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def clues2csv(clues_files, chrom, win_start, win_end, pop):
    clues_csv_file_name = f'steps/clues/csv/clues_{chrom}_{win_start}_{win_end}_{pop}.csv'
    clues_file_base_names = ' '.join([modpath(f, parent='', suffix='') for f in clues_files])

    inputs = clues_files
    outputs = [clues_csv_file_name]
    options = {
        'memory': '1g',
        'walltime': '1:00:00'
    }

    spec = f'''

    mkdir -p steps/clues/csv

    python scripts/extract_clues_info.py {clues_csv_file_name} steps/clues/clues_{chrom}_{win_start}_{win_end}_{pop} {clues_file_base_names}

    '''

    return AnonymousTarget(inputs=inputs, outputs=outputs, options=options, spec=spec)

def execute(cmd, stdin=None):
    process = Popen(cmd.split(), stdin=PIPE, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate(stdin)
    return stdout, stderr

def read_snp_info(snp_file):
    snp_list = list()
    with open('snps.txt', 'r') as snp_file:
        for line in snp_file:
            chrom, snp_pos, derived_allele, derived_freq = line.split()
            snp_pos = int(snp_pos)
            derived_freq = float(derived_freq)
            snp_list.append((chrom, snp_pos, derived_allele, derived_freq))
    return snp_list
  
def get_single_snp(freq_data_file, chrom, pop, snp_pos):
    snp_file_name = 'snps.txt'
    execute(f"python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom} {pop} {snp_file_name} --snppos {snp_pos}")
    snp_list = read_snp_info(snp_file_name)
    return snp_list


def get_snps(freq_data_file, chrom, pop, window_start, window_end, min_freq, nr_snps):
    snp_file_name = 'snps.txt'
    execute(f"python ./scripts/get_derived_freq_data.py {freq_data_file} {chrom} {pop} {snp_file_name} --start {window_start} --end {window_end} --minfreq {min_freq} --nrsnps {nr_snps}")
    snp_list = read_snp_info(snp_file_name)
    return snp_list

freq_data_file = 'data/derived_pop_freqs.h5'
overlap_of_windows = 200000
min_freq = 0.25
nr_snps = 100

vcf_file = 'data/ALL.chr3.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz'
chrom = '3'
windows = [ [45000000,46400000],
            [46000000,47400000],
            [47000000,48400000],
            [48000000,49400000],
            [49000000,50400000],
            [50000000,51400000],
            [51000000,52400000],
            [52000000,53400000],
            [53000000,54400000],
            [54000000,55400000]]

"""
vcf_file = 'data/chr2.vcf.gz'
chrom = '2'
windows = [ [136000000,137000000],
            [136600000,137600000]]
"""

for pop in ['CEU', 'FIN', 'GBR', 'LWK', 'YRI']:
    for window_start, window_end in windows:

        # selected snps in window:
        snp_list = get_snps(freq_data_file, chrom, pop, window_start+overlap_of_windows, window_end-overlap_of_windows, min_freq, nr_snps)

        # lactase snp:
        if window_start < 136608646-overlap_of_windows and window_end > 136608646+overlap_of_windows:
            snp_list.extend(get_single_snp(freq_data_file, chrom, pop, 136608646))

        gwf.target_from_template(
            name=f'extract_window_{chrom}_{window_start}_{window_end}_{pop}',
            template=extract_window(
                with_indels=vcf_file,
                chrom = chrom,
                win_start = window_start,
                win_end = window_end,
                pop = pop
            )
        )

        gwf.target_from_template(
            name=f'reformat_vcf_{chrom}_{window_start}_{window_end}_{pop}',
            template=vcf2sites(
                chrom = chrom,
                win_start = window_start,
                win_end = window_end,
                pop = pop
            )
        )

        gwf.target_from_template(
            name=f'argsample_{chrom}_{window_start}_{window_end}_{pop}',
            template=argsample(
                chrom = chrom,
                win_start = window_start,
                win_end = window_end,
                pop = pop
            )
        )

        clues_task_list = list()
        for chrom, snp_pos, derived_allele, derived_freq in snp_list:
            clues_task = gwf.target_from_template(
                name=f'clues_{chrom}_{window_start}_{window_end}_{snp_pos}_{pop}',
                template=clues(
                    snp_position = snp_pos,
                    chrom = chrom,
                    win_start = window_start,
                    win_end = window_end,
                    pop = pop,
                    derived_allele = derived_allele,
                    derived_freq = derived_freq
                )
            )
            clues_task_list.append(clues_task)
            clues_files = [output for task in clues_task_list for output in task.outputs]

        gwf.target_from_template(
            name=f'csv_{chrom}_{window_start}_{window_end}_{pop}',
            template=clues2csv(
                clues_files = clues_files,
                chrom = chrom,
                win_start = window_start,
                win_end = window_end,
                pop = pop
            )
        )
