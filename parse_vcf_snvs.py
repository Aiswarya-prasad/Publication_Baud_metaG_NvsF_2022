import os
import sys
import time
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import combinations, chain
from collections import Counter
import math
import statistics as st
import multiprocessing
import pickle
import random

'''
parse and summarize vcf files into the files needed for figure 2
+ {sdp}_sample_var_host_{threshold}.txt > concatenated into all_samples_var_host_{threshold}.txt - 2A
+ {sdp}_cum_curve_{threshold}.txt - 2B
+ {sdp}_shared_fraction_{threshold}.txt - 2C

The paths for the files are hard coded in the script so that will need to be
changed if the script is run on a different machine
'''

snvs_analysis_directory = 'snvs' # this is where the vcf files are supposed to be

all_samples = ['N01','N02','N03','N04','N06','N07','N08','N09','N10','N11','N12','N13','N14','N15','N16',
                  'F01','F02','F03','F04','F06','F07','F08','F09','F10','F11','F12','F13','F14','F15','F16']

class GenomePosition:
    def __init__(self, position):
        self.position = position # marks the position
        self.depth = {} # dict of samples and their depth at this position
        self.samples = [] # list of samples that for this position (includes missing data)
        self.alleles = {} # dict of samples and a list of their alleles at this position
        self.count = {} # dict of samples and a dict of alleles and their counts
        self.freq = {} # dict of samples and a dict of alleles and their frequency
        self.missingfrac = None # this can only be calculated if needed after all samples have been added

    def add_sample_info(self, sample, depth, depth_diff, allele, freq, missing=None):
        '''
        samples for which allele freq is non-zero
        are added to the position object
        freq is actually the count freq_calc is the ratio of the counts of
        allele and depth
        if allele freq is missing, it is not added to the position object
        '''
        if freq == '0':
            return True
        if depth == '.' and freq == '.':
            self.samples.append(sample)
            self.alleles[sample] = ['.']
            self.missing_frac = missing
            self.freq[sample] = {}
            self.freq[sample]['.'] = -1
            self.count[sample] = {}
            self.count[sample]['.'] = -1
            return True
        else:
            depth = int(depth)
            freq_calc = float(freq)/int(depth-depth_diff)
            if sample not in self.samples:
                self.samples.append(sample)
            self.depth[sample] = depth
            self.missing_frac = missing
            if sample in self.freq.keys():
                self.freq[sample][allele] = float(freq_calc)
            else:
                self.freq[sample] = {}
                self.freq[sample][allele] = float(freq_calc)
            if sample in self.count.keys():
                self.count[sample][allele] = int(freq)
            else:
                self.count[sample] = {}
                self.count[sample][allele] = int(freq)
            if sample in self.alleles.keys():
                if allele not in self.alleles[sample]:
                    self.alleles[sample].append(allele)
            else:
                self.alleles[sample] = [allele]
            return True

    def count_allele(self, sample, allele):
        if allele not in self.alleles[sample]:
            return 0
        return self.count[sample][allele]

    def freq_allele(self, sample, allele):
        if allele not in self.alleles[sample]:
            return 0
        return self.freq[sample][allele]
 
    def depth_allele(self, sample, allele):
        if allele not in self.alleles[sample]:
            return 0
        return self.depth[sample][allele]

    def __str__(self):
        return f'position: {self.position}\n' + \
                f'alleles: \n ' + \
                f'{self.count}\n'
    
    def __repr__(self):
        return f'position: {self.position}\n' + \
                f'alleles: \n ' + \
                f'{self.count}\n'
 
    def is_polymorphic(self, sample, cutoff=0.1, assymetric=False, consider_counts=True):
        '''
        if there is no cutoff, consider a position polymorphic if 
        there are at least 2 alleles present and if there is a cutoff,
        consider a position polymorphic if there are at least 2 alleles
        and both have a frequency greater than the cutoff
        This is not the same as the is_polymorphic subroutin in filter_snvs.pl
        because it considers polymorphic positions that are informative i.e.
        not all 0 or 1 for the allele across samples.
        Also, those that are above the cutoff but only supported by 1 read
        will be removed (set consider_counts False to avoid this)
        '''
        if self.alleles[sample] == ['.']:
            return False
        if round(sum(self.freq[sample].values()), 7) != 1:
            print(f'Allele freqs do not add up to 1 for {self.position} in {sample}')
            print(f'Allele freqs: {self.freq[sample]}')
            print(f'Allele counts: {self.freq[sample]}')
        # now that we know stuff adds up to 1, we want to know that there is
            # polymorphism at this position
            # regardless of how many alleles there are, at least one allele
            # should have a freq >= cutoff and there should be more than 1 allele
            # with a non-zero freq - this is not explicitly checked because if freq
            # is 0, it is not added to the position object
        if assymetric:
            if len(self.alleles[sample]) > 1:
                # if there is more than 1 allele and there is no cutoff, return True
                if cutoff == 0:
                    return True
                else:
                    if [x for x in self.freq[sample].values()][1] >= cutoff and [x for x in self.freq[sample].values()][0] >= 0.05:
                        if consider_counts:
                            # This checks if at least 2 or more alleles in this position have a count > 1
                            if sum([self.count[sample][an_allele] > 1 for an_allele in self.alleles[sample]]) >= 2:
                                # print(f'count of minor allele is {self.count[sample][an_allele]}/{self.depth[sample]} for {self.position} in {sample}')
                                return True
                            else:
                                return False
                        return True
                    else:
                        return False
            else:
                # print('not counting as polymorphic')
                return False
        else:
            if len(self.alleles[sample]) > 1:
                # if there is more than 1 allele and there is no cutoff, return True
                if cutoff == 0:
                    return True
                else:
                    if sum([x >= cutoff for x in self.freq[sample].values()]) >= 2:
                        if consider_counts:
                            # This checks if at least 2 or more alleles in this position have a count > 1
                            if sum([self.count[sample][an_allele] > 1 for an_allele in self.alleles[sample]]) >= 2:
                                # print(f'count of minor allele is {self.count[sample][an_allele]}/{self.depth[sample]} for {self.position} in {sample}')
                                return True
                            else:
                                return False
                        return True
                    else:
                        return False
            else:
                # print('not counting as polymorphic')
                return False

    def alleles_detected(self, sample, cutoff=0.1):
        '''
        return a list of alleles detected for a sample at this position
        at a frequency >= cutoff
        '''
        alleles_detected = []
        for a in self.alleles[sample]:
            if self.freq[sample][a] >= cutoff and self.count[sample][a] > 1:
                alleles_detected.append(a)
        return alleles_detected

    def get_dominant_allele(self, sample):
        '''
        return the allele with the highest frequency
        '''
        # if any alleles have exactly the same freq
        # print a warning message about them
        # and return the first one
        if '.' in self.alleles[sample]:
            return ''
        dom_allele = max(self.freq[sample], key=self.freq[sample].get)
        if len(set(self.freq[sample].values())) != len(self.freq[sample].values()) and self.freq[sample][dom_allele] in [self.freq[sample][x] for x in self.freq[sample].keys() if x != dom_allele]:
            print(f'Warning: multiple alleles with the same frequency at {self.position} in {sample}')
            print(f'Allele freqs: {self.freq[sample]}')
            print(f'Allele counts: {self.count[sample]}')
            print(f'choosing as dominant: {dom_allele}')
        return dom_allele
    
    def get_dominant_allele_relaxed(self, sample):
        '''
        return the allele with the highest frequency
        if there are two alleles with the same freq
        and that happens to be those of the highest frequnecy
        consider both to be dominant!
        '''
        if '.' in self.alleles[sample]:
            return []
        dom_allele = max(self.freq[sample], key=self.freq[sample].get)
        dom_alleles = []
        if len(set(self.freq[sample].values())) != len(self.freq[sample].values()) and self.freq[sample][dom_allele] in [self.freq[sample][x] for x in self.freq[sample].keys() if x != dom_allele]:
            dom_alleles = [x for x in self.freq[sample].keys() if self.freq[sample][x] == self.freq[sample][dom_allele]]
        else:
            dom_alleles = [dom_allele]
        return dom_alleles
    
    def polymorphic_for(self, cutoff=0.1):
        '''
        return a list of samples for which the position is polymorphic
        '''
        return [sample for sample in self.samples if self.is_polymorphic(sample, cutoff=cutoff)]


def populate_position_objects(vcf_file, all_samples, missing_fraction_cutoff=0.1):
        # vcf_file = f'{snvs_analysis_directory}/{sdp}.vcf.gz'
        position_objs = {}
        with gzip.open(vcf_file, 'rb') as fh:
            for line in fh:
                if line.startswith(b'##'):
                    continue
                if line.startswith(b'#CHROM'):
                    samples = line.strip().split(b'\t')[9:]
                    samples = [x.decode() for x in samples]
                    sample_indices = [9+i for i, x in enumerate(samples) if x in all_samples]
                else:
                    line = line.strip().split(b'\t')
                    line = [x.decode() for x in line]
                    # parse position and allele info from the vcf file line
                    position = line[1]
                    ref_allele = line[3]
                    alt_allele = line[4]
                    # filter for missing fraction
                    tot = 0
                    missing = 0
                    for i in range(9, len(line)):
                        if line[i].split(':')[5] == '.':
                            missing += 1
                        tot += 1
                    missing_fraction = round(missing / tot, 2)
                    if missing_fraction >= missing_fraction_cutoff:
                        continue
                    for i in range(9, len(line)):
                    # for i in sample_indices:
                        sample = samples[i-9]
                        if position not in position_objs:
                            position_objs[position] = GenomePosition(position)
                        # 6th (index 5) field is AO for the allele count,
                        # 2nd (index 1) field is DP for depth at the position
                        # missing positions are represented by a dot and not added to the position object
                        # in KE's workflow positions with >= 0.1 missing fraction are removed
                        # this keeps the data comparable so "oversequenced" samples are not overestimated
                        # in the number of polymorphic sites
                        # mark positions with too much (>=0.1) missing data this is the fractions of samples
                        # that have a dot in the AO field for that position
                        # Sometimes allelic freqs do not add up to depth. This could result in
                        # a position being considered polymorphic when it is not. This is because
                        # even when the other allele count is 0, the freq of the allele considered
                        # would not be 1 but something less than 1!
                        # for this, define a new attribute called depth_diff
                        # use depth-depth_diff as the depth to calculate the allele freqs
                        depth_diff = 0
                        if line[i].split(':')[1] != '.':
                            ad_nums = [int(x) for x in line[i].split(':')[2].split(',')]
                            depth_diff = int(line[i].split(':')[1]) - sum(ad_nums)
                        position_objs[position].add_sample_info(sample, line[i].split(':')[1], depth_diff, ref_allele, line[i].split(':')[3], missing_fraction)
                        position_objs[position].add_sample_info(sample, line[i].split(':')[1], depth_diff, alt_allele, line[i].split(':')[5], missing_fraction)
        return position_objs

core_lengths_df = pd.read_csv('04_RebuttalAnalyses/05_CoreGenomeLengthsInfo/core_lengths_final_red.txt', sep='\t', header=0)

sdp_list = ["bapis","bifido_1.1","bifido_1.2","bifido_1.3","bifido_1.4",
                         "bifido_1.5","bifido_2","firm4_1","firm4_2","firm5_1",
                         "firm5_2","firm5_3","firm5_4","fper_1","gilli_1",
                         "gilli_2","snod_1"]

for sdp in sdp_list:
    print('Working on', sdp)
    positions = populate_position_objects(f'{snvs_analysis_directory}/{sdp}.vcf.gz', all_samples)
    if not os.path.exists('{snvs_analysis_directory}/positions_dicts/'):
        with open(f'{snvs_analysis_directory}/{sdp}_positions_dict.pkl', 'wb') as fh_out:
            pickle.dump(positions, fh_out)

for sdp in sdp_list:
    print('Working on', sdp)
    if os.path.exists('{snvs_analysis_directory}/{sdp}_positions_dict.pkl'):
        positions = pickle.load(open(f'{snvs_analysis_directory}/{sdp}_positions_dict.pkl', 'rb'))

    # read core lengths
    core_lengths_df[core_lengths_df['Strain_id'] == sdp]['Tot_length(bp)'].values[0]

    # make the {sdp}_sample_var_host.txt file for all of these cases
    with open(f'{snvs_analysis_directory}/filter_0/{sdp}_sample_var_host.txt', 'w') as fh_out:
        fh_out.write('phylo\tstrain\thost\tcolony\tsample\tnb_var_pos\tfraction\n')
        phylo = sdp.split('_')[0]
        strain = sdp
        for sample in all_samples:
            host = sample[0]
            colony = sample[1:]
            nb_var_pos = 0
            for loc in positions:
                if sample in positions[loc].samples:
                    if positions[loc].is_polymorphic(sample, cutoff=0):
                        nb_var_pos += 1
            fraction = round(nb_var_pos / core_lengths_df[core_lengths_df['Strain_id'] == sdp]['Tot_length(bp)'].values[0] * 100, 2)
            fh_out.write(f'{phylo}\t{strain}\t{host}\t{colony}\t{sample}\t{nb_var_pos}\t{fraction}\n')

    with open(f'{snvs_analysis_directory}/filter_01/{sdp}_sample_var_host.txt', 'w') as fh_out:
        fh_out.write('phylo\tstrain\thost\tcolony\tsample\tnb_var_pos\tfraction\n')
        phylo = sdp.split('_')[0]
        strain = sdp
        for sample in all_samples:
            host = sample[0]
            colony = sample[1:]
            nb_var_pos = 0
            for loc in positions:
                if sample in positions[loc].samples:
                    if positions[loc].is_polymorphic(sample, cutoff=0.01):
                        nb_var_pos += 1
            fraction = round(nb_var_pos / core_lengths_df[core_lengths_df['Strain_id'] == sdp]['Tot_length(bp)'].values[0] * 100, 2)
            fh_out.write(f'{phylo}\t{strain}\t{host}\t{colony}\t{sample}\t{nb_var_pos}\t{fraction}\n')
    
    with open(f'{snvs_analysis_directory}/filter_1/{sdp}_sample_var_host.txt', 'w') as fh_out:
        fh_out.write('phylo\tstrain\thost\tcolony\tsample\tnb_var_pos\tfraction\n')
        phylo = sdp.split('_')[0]
        strain = sdp
        for sample in all_samples:
            host = sample[0]
            colony = sample[1:]
            nb_var_pos = 0
            for loc in positions:
                if sample in positions[loc].samples:
                    if positions[loc].is_polymorphic(sample, cutoff=0.1):
                        nb_var_pos += 1
            fraction = round(nb_var_pos / core_lengths_df[core_lengths_df['Strain_id'] == sdp]['Tot_length(bp)'].values[0] * 100, 2)
            fh_out.write(f'{phylo}\t{strain}\t{host}\t{colony}\t{sample}\t{nb_var_pos}\t{fraction}\n')

    with open(f'{snvs_analysis_directory}/filter_05/{sdp}_sample_var_host.txt', 'w') as fh_out:
        fh_out.write('phylo\tstrain\thost\tcolony\tsample\tnb_var_pos\tfraction\n')
        phylo = sdp.split('_')[0]
        strain = sdp
        for sample in all_samples:
            host = sample[0]
            colony = sample[1:]
            nb_var_pos = 0
            for loc in positions:
                if sample in positions[loc].samples:
                    if positions[loc].is_polymorphic(sample, cutoff=0.05):
                        nb_var_pos += 1
            fraction = round(nb_var_pos / core_lengths_df[core_lengths_df['Strain_id'] == sdp]['Tot_length(bp)'].values[0] * 100, 2)
            fh_out.write(f'{phylo}\t{strain}\t{host}\t{colony}\t{sample}\t{nb_var_pos}\t{fraction}\n')


for sdp in sdp_list:
    df_plot_vars_0 = pd.read_csv(f'{snvs_analysis_directory}/filter_0/{sdp}_sample_var_host.txt', sep='\t', header=0)
    df_plot_vars_01 = pd.read_csv(f'{snvs_analysis_directory}/filter_01/{sdp}_sample_var_host.txt', sep='\t', header=0)
    df_plot_vars_1 = pd.read_csv(f'{snvs_analysis_directory}/filter_1/{sdp}_sample_var_host.txt', sep='\t', header=0)
    df_plot_vars_05 = pd.read_csv(f'{snvs_analysis_directory}/filter_05/{sdp}_sample_var_host.txt', sep='\t', header=0)
    df_plot_vars_KE = pd.read_csv(f'snvs/{sdp}_sample_var_host.txt', sep='\t', header=None)
    df_plot_vars_KE.columns = ['phylo', 'strain', 'host', 'colony', 'sample', 'nb_var_pos', 'fraction']
    # only keep the samples that are in all_samples
    df_plot_vars_KE =  df_plot_vars_KE[df_plot_vars_KE['sample'].isin(all_samples)]
    # combine the dataframes by filter condition
    df_plot_vars_0['cutoff'] = '0'
    df_plot_vars_01['cutoff'] = '0.01'
    df_plot_vars_1['cutoff'] = '0.1'
    df_plot_vars_05['cutoff'] = '0.05'
    df_plot_vars_KE['cutoff'] = 'KE'
    df_plot_vars = pd.concat([df_plot_vars_0, df_plot_vars_01, df_plot_vars_1, df_plot_vars_05, df_plot_vars_KE])

    df_plot_vars['nb_var_pos'] = df_plot_vars['nb_var_pos'].astype(int) 
    df_plot_vars['fraction'] = df_plot_vars['fraction'].astype(float)
    df_plot_vars['cutoff'] = df_plot_vars['cutoff'].astype(str)
    # only keep those with non-zero nb_var_pos
    df_plot_vars = df_plot_vars[df_plot_vars['nb_var_pos'] > 0]
    # plt the data - grouby F and N host next to each other and cutoff on x axis and fraction on y axis
    plt.style.use('ggplot')
    # add a marker to show the mean for each group
    # get the mean for each group
    means = df_plot_vars.groupby(['cutoff', 'host'])['fraction'].mean().reset_index()
    means['cutoff_index'] = means['cutoff'].apply(lambda x: ['0', '0.01', '0.1', '0.05', 'KE'].index(x))
    fig, ax = plt.subplots(figsize=(10, 10))
    sns.stripplot(x='cutoff', y='fraction', hue='host', 
                data=df_plot_vars, ax=ax, jitter=0.1, dodge=True,
                s = 8, alpha=0.85, zorder=1)
    ax.set_xlabel('Cutoff')
    ax.set_ylabel('Fraction variable sites')
    # specify colors '#D55E00' for N and '#56B4E9' for F
    ax.legend(loc='upper left', title='Host', labels=['N', 'F'])
    ax.set_title(f'{sdp}')
    plt.ylim(0, math.ceil(max(df_plot_vars['fraction']) + 1))
    # Add a marker for each mean, with dodge for different hosts
    for i, row in means.iterrows():
        cutoff_index = row['cutoff_index']
        cutoff = row['cutoff']
        host = row['host']
        fraction = row['fraction']
        dodge_amount = -0.2 if host == 'N' else 0.2  # Determine dodge amount based on host
        ax.scatter(x=cutoff_index + dodge_amount, y=fraction,
                   marker='+', color='black', s=100, zorder=1)
    plt.savefig(f'{snvs_analysis_directory}/figures/{sdp}_var_sites_vs_cutoff.png', dpi=300)

# def calc_shared_KE_approach(positions, sample1, sample2, cutoff_detect):
#     distance = 1
#     shared = 0
#     total = 0
#     for loc in positions:
#         if any([positions[loc].is_polymorphic(sam, cutoff=cutoff_detect) for sam in all_samples if sam in positions[loc].samples]):
#             if sample1 in positions[loc].samples or sample2 in positions[loc].samples:
#             # if an allele is present in at least 1 sample, consider it for the total
#                 total += 1
#                 alleles_1 = set(positions[loc].alleles_detected(sample1, cutoff=cutoff_detect))
#                 alleles_2 = set(positions[loc].alleles_detected(sample2, cutoff=cutoff_detect))
#                 alleles_1.discard('.')
#                 alleles_2.discard('.')
#                 # if either sample contains an allele that is not present in the other
#                 # consider it not shared
#                 if len(alleles_1) > 0 and len(alleles_2) > 0:
#                     if len(alleles_1.difference(alleles_2)) == 0 or len(alleles_2.difference(alleles_1)) == 0:
#                         shared += 1
#     if total > 0:
#         distance = 1 - shared / total
#         if shared == 0:
#             print(f'setting 0 for {sample1} and {sample2} cut off at {cutoff_detect} because shared is 0!!!!!!!')
#     else:
#         distance = 0
#         # print(f'setting 0 for {sample1} and {sample2} cut off at {cutoff_detect} because total is 0')
#     return (shared, total)

def calc_shared_dominant(positions, sample1, sample2, cutoff_detect):
    distance = 1
    shared = 0
    total = 0
    for loc in positions:
        if sample1 in positions[loc].samples and sample2 in positions[loc].samples:
            total += 1
            dom_alleles_1 = set(positions[loc].get_dominant_allele_relaxed(sample1))
            dom_alleles_2 = set(positions[loc].get_dominant_allele_relaxed(sample2))
            if len(dom_alleles_1.intersection(dom_alleles_2)) >= 1:
                # and all alleles in the intersection are above the cutoff in each sample
                if all([positions[loc].freq_allele(sample1, x) >= cutoff_detect for x in dom_alleles_1.intersection(dom_alleles_2)]) and all([positions[loc].freq_allele(sample2, x) >= cutoff_detect for x in dom_alleles_1.intersection(dom_alleles_2)]):
                    shared += 1
                else:
                    print(f'for {loc} in {sample1} and {sample2}, dominant allele is {dom_alleles_1} and {dom_alleles_2} but alleles detected are {set(positions[loc].alleles_detected(sample1, cutoff=0.1))} and {set(positions[loc].alleles_detected(sample2, cutoff=0.1))}')
            # dom_allele_1 = positions[loc].get_dominant_allele(sample1)
            # dom_allele_2 = positions[loc].get_dominant_allele(sample2)
            # if dom_allele_1 == dom_allele_2:
            #     if dom_allele_1 in set(positions[loc].alleles_detected(sample1, cutoff=cutoff_detect)) and dom_allele_2 in set(positions[loc].alleles_detected(sample2, cutoff=cutoff_detect)):
            #         shared += 1
            #     else:
            #         print(f'for {loc} in {sample1} and {sample2}, dominant allele is {dom_allele_1} and {dom_allele_2} but alleles detected are {set(positions[loc].alleles_detected(sample1, cutoff=0.1))} and {set(positions[loc].alleles_detected(sample2, cutoff=0.1))}')
    if total > 0:
        distance = 1 - shared / total
        if shared == 0:
            print(f'setting 0 for {sample1} and {sample2} cut off at {cutoff_detect} because shared is 0!!!!!!!')
    else:
        distance = 0
        # print(f'setting 0 for {sample1} and {sample2} cut off at {cutoff} because total')
    return (shared, total)

# # /home/aiswarya/mnt/nas_recherche/general_data/D2c/gbaud/20190625_gbaud_beemetagenomics/{snvs_analysis_directory}/bapis_shared_fraction_filt_01.txt
# # were generated using the following function

# WE USE THIS FOR THE PUBLICATION

# def calc_shared_KE_approach(positions, sample1, sample2, cutoff_detect):
#     shared = 0
#     total = 0
#     for loc in positions:
#         if sample1 in positions[loc].samples and sample2 in positions[loc].samples:
#             # for outputs labelled _corrected this condition is and
#             # before, it was or (for the various cutoffs)
#             if positions[loc].is_polymorphic(sample1, cutoff=cutoff_detect) or positions[loc].is_polymorphic(sample2, cutoff=cutoff_detect):
#                 total += 1
#                 alleles_1 = set(positions[loc].alleles_detected(sample1, cutoff=cutoff_detect))
#                 alleles_2 = set(positions[loc].alleles_detected(sample2, cutoff=cutoff_detect))
#                 alleles_1.discard('.')
#                 alleles_2.discard('.')
#                 # for outputs labelled _corrected this condition is >= 1
#                 # before, it was > 1
#                 to_print_1 = [f'{x}:{positions[loc].freq_allele(sample1, x)}' for x in alleles_1]
#                 to_print_2 = [f'{x}:{positions[loc].freq_allele(sample2, x)}' for x in alleles_2]
#                 print(f'{loc}\t{sample1}\t{sample2}\t{to_print_1}\t{to_print_2}')
#                 if len(alleles_1.intersection(alleles_2)) > 1:
#                     shared += 1

# # for _corrected we use:

# def calc_shared_KE_approach(positions, sample1, sample2, cutoff_detect):
#     shared = 0
#     total = 0
#     for loc in positions:
#         if sample1 in positions[loc].samples and sample2 in positions[loc].samples:
#             # for outputs labelled _corrected this condition is and
#             # before, it was or (for the various cutoffs)
#             if positions[loc].is_polymorphic(sample1, cutoff=cutoff_detect) or positions[loc].is_polymorphic(sample2, cutoff=cutoff_detect):
#                 total += 1
#                 alleles_1 = set(positions[loc].alleles_detected(sample1, cutoff=cutoff_detect))
#                 alleles_2 = set(positions[loc].alleles_detected(sample2, cutoff=cutoff_detect))
#                 alleles_1.discard('.')
#                 alleles_2.discard('.')
#                 # for outputs labelled _corrected this condition is >= 1
#                 # before, it was > 1
#                 to_print_1 = [f'{x}:{positions[loc].freq_allele(sample1, x)}' for x in alleles_1]
#                 to_print_2 = [f'{x}:{positions[loc].freq_allele(sample2, x)}' for x in alleles_2]
#                 print(f'{loc}\t{sample1}\t{sample2}\t{to_print_1}\t{to_print_2}')
#                 if len(alleles_1.intersection(alleles_2)) >= 1:
#                     shared += 1
#    return (shared, total)

# # for _final we use:

# def calc_shared_KE_approach(positions, sample1, sample2, cutoff_detect):
#     shared = 0
#     total = 0
#     for loc in positions:
#         if any([positions[loc].is_polymorphic(sam, cutoff=cutoff_detect) for sam in all_samples if sam in positions[loc].samples]):
#             total += 1
#             # if position is present in both sample, it can be compared
#             # if is not poly in at least 1 sample, consider it for counting
#             alleles_1 = set(positions[loc].alleles_detected(sample1, cutoff=cutoff_detect))
#             alleles_2 = set(positions[loc].alleles_detected(sample2, cutoff=cutoff_detect))
#             alleles_1.discard('.')
#             alleles_2.discard('.')
#             # avoid cases where there is a missing allele in either sample
#             # nothing left after discarding missing alleles
#             if len(alleles_1) == 0 or len(alleles_2) == 0:
#                 continue
#             if positions[loc].is_polymorphic(sample1, cutoff=cutoff_detect) and positions[loc].is_polymorphic(sample2, cutoff=cutoff_detect):
#                 shared += 1
#     return (shared, total)


# for final_all we use:

def calc_shared_KE_approach(positions, sample1, sample2, cutoff_detect):
    shared = 0
    total = 0
    for loc in positions:
        if sample1 in positions[loc].samples and sample2 in positions[loc].samples:
            # if position is present in both sample, it can be compared
            # if is not poly in at least 1 sample, consider it for counting
            alleles_1 = set(positions[loc].alleles_detected(sample1, cutoff=cutoff_detect))
            alleles_2 = set(positions[loc].alleles_detected(sample2, cutoff=cutoff_detect))
            alleles_1.discard('.')
            alleles_2.discard('.')
            # avoid cases where there is a missing allele in either sample
            # nothing left after discarding missing alleles
            if len(alleles_1) == 0 or len(alleles_2) == 0:
                continue
            if positions[loc].is_polymorphic(sample1, cutoff=cutoff_detect) or positions[loc].is_polymorphic(sample2, cutoff=cutoff_detect):
                total += 1
                if positions[loc].is_polymorphic(sample1, cutoff=cutoff_detect) and positions[loc].is_polymorphic(sample2, cutoff=cutoff_detect):
                    shared += 1
                '''
                It has never happened that there is only 1 shared allele
                because if both the sites are polymorphic, they are always
                polymorphic for at least 2 alleles if this is not true
                that line will be printed.
                '''
                if alleles_1.intersection(alleles_2) == 1:
                    print(f'loc: {loc}, sample1: {sample1}, sample2: {sample2}, alleles_1: {[f"{al}:{positions[loc].freq_allele(sample1, al)}" for al in alleles_1]}, alleles_2: {[f"{al}:{positions[loc].freq_allele(sample2, al)}" for al in alleles_2]}, shared: {alleles_1.intersection(alleles_2)}')
    return (shared, total)

# Now that we have the distance matrices, we can plot the cumulative curves
# for each sdp and separate by host
'''
reworking distance matrix calc. Instead of computing distance matrix, make a list like
KEs file ...shared_fraction.txt
each line has 
Sample1	Sample2	Nb_shared_alleles	Nb_scored_alleles	Shared_fraction

then use snvs/distance_matrix.R to write the matrices

run the following for each sdp
cd snvs
Rscript distance_matrix.R ${sdp}_shared_fraction_filt_0.txt ${sdp}_dist_matrix_filt_0.txt
Rscript distance_matrix.R ${sdp}_shared_fraction_filt_01.txt ${sdp}_dist_matrix_filt_01.txt
Rscript distance_matrix.R ${sdp}_shared_fraction_filt_05.txt ${sdp}_dist_matrix_filt_05.txt
Rscript distance_matrix.R ${sdp}_shared_fraction_filt_1.txt ${sdp}_dist_matrix_filt_1.txt

Rscript distance_matrix.R ${sdp}_shared_fraction_dom_0.txt ${sdp}_dist_matrix_dom_0.txt
# cutoff makes not much sense for dominant alleles
'''

# for sdp in ["bapis","bifido_1.1","bifido_1.2","bifido_1.3","bifido_1.4","bifido_1.5","bifido_2","firm4_1","firm4_2"]:
# for sdp in ["firm5_1", "firm5_2","firm5_3","firm5_4","fper_1","gilli_1", "gilli_2","snod_1"]:
for sdp in sdp_list:
    print('Working on', sdp)
    if os.path.exists(f'{snvs_analysis_directory}/{sdp}_positions_dict.pkl'):
        positions = pickle.load(open(f'{snvs_analysis_directory}/{sdp}_positions_dict.pkl', 'rb'))
    sample_pairs = list(combinations(all_samples, 2))
    # for each sample pair, calculate the distance
    # and write to a file
    # for each cutoff
    for cutoff, name in zip([0.1], ['1']):
    # for cutoff, name in zip([0.01, 0, 0.05, 0.1], ['01', '0', '05', '1']):
        print('Working on cutoff', cutoff)
        if True:
        # with open(f'{snvs_analysis_directory}/{sdp}_shared_fraction_filt_{name}_corrected.txt', 'w') as fh_out:
            # fh_out.write('Sample1\tSample2\tNb_shared_alleles\tNb_scored_alleles\tShared_fraction\n')
            for sample1, sample2 in sample_pairs:
                shared, total = calc_shared_KE_approach(positions, sample1, sample2, cutoff_detect = cutoff)
                if shared > 0 and total > 0:
                    shared_fraction = shared / total
                    # fh_out.write(f'{sample1}\t{sample2}\t{shared}\t{total}\t{shared_fraction}\n')
        if cutoff == 0:
            print('Working on dom')
            with open(f'{snvs_analysis_directory}/{sdp}_shared_fraction_dom_{name}.txt', 'w') as fh_out:
                fh_out.write('Sample1\tSample2\tNb_shared_alleles\tNb_scored_alleles\tShared_fraction\n')
                for sample1, sample2 in sample_pairs:
                    shared, total = calc_shared_dominant(positions, sample1, sample2, cutoff_detect = cutoff)
                    if shared > 0 and total > 0:
                        shared_fraction = shared / total
                        fh_out.write(f'{sample1}\t{sample2}\t{shared}\t{total}\t{shared_fraction}\n')
 
# how often does it happen that the sites are not polymorphic but they do not contain the same
# allele between two samples?

# for sdp in sdp_list:
#     if os.path.exists(f'{snvs_analysis_directory}/{sdp}_positions_dict.pkl'):
#         positions = pickle.load(open(f'{snvs_analysis_directory}/{sdp}_positions_dict.pkl', 'rb'))
#     sample_pairs = list(combinations(all_samples, 2))
#     cutoff = 0.01
#     name = '01'
    
# for every position, if the allele is present in at least 1 sample in a pair
# consider it for the total
# if it is present in both samples, consider it for the shared
# if it is present in only 1 sample, do not consider it for the shared
# if it is present in neither sample, do not consider it for the total

for sdp in sdp_list:
    print('Working on', sdp)
    if os.path.exists(f'{snvs_analysis_directory}/{sdp}_positions_dict.pkl'):
        positions = pickle.load(open(f'{snvs_analysis_directory}/{sdp}_positions_dict.pkl', 'rb'))
    # sample_pairs = list(combinations(all_samples, 2))
    # sample_pairs = [x for x in sample_pairs if coverage_in_sample(x[0], sdp) and coverage_in_sample(x[1], sdp)]
    # for each sample pair, calculate the distance
    # and write to a file
    # for each cutoff
    for cutoff, name in zip([0.01], ['01']):
    # for cutoff, name in zip([0.01, 0, 0.05, 0.1], ['01', '0', '05', '1']):
        print('Working on cutoff', cutoff)
        if True:
        # with open(f'{snvs_analysis_directory}/{sdp}_shared_fraction_filt_{name}_corrected.txt', 'w') as fh_out:
            # fh_out.write('Sample1\tSample2\tNb_shared_alleles\tNb_scored_alleles\tShared_fraction\n')
            for sample1, sample2 in sample_pairs:
                inter_geq1 = 0
                inter_g1 = 0
                poly = 0
                poly_both_int0 = 0
                poly_both = 0
                present = 0
                for loc in positions:
                    if sample1 in positions[loc].samples and sample2 in positions[loc].samples:
                        present += 1
                        if positions[loc].is_polymorphic(sample1, cutoff=cutoff) or positions[loc].is_polymorphic(sample2, cutoff=cutoff):
                            poly += 1
                            alleles_1 = set(positions[loc].alleles_detected(sample1, cutoff=cutoff))
                            alleles_2 = set(positions[loc].alleles_detected(sample2, cutoff=cutoff))
                            if len(alleles_1.intersection(alleles_2)) > 1:
                                inter_g1 += 1
                            if len(alleles_1.intersection(alleles_2)) >= 1:
                                inter_geq1 += 1
                            if positions[loc].is_polymorphic(sample1, cutoff=cutoff) and positions[loc].is_polymorphic(sample2, cutoff=cutoff):
                                poly_both += 1
                                if len(alleles_1.intersection(alleles_2)) == 0:
                                    poly_both_int0 += 1
                print(f'{sample1}\t{sample2}\t{present:,}: present')
                print(f'{sample1}\t{sample2}\t{poly:,}: poly')
                print(f'{sample1}\t{sample2}\t{inter_geq1:,}: inter_geq1')
                print(f'{sample1}\t{sample2}\t{poly_both:,}: poly_both')
                print(f'{sample1}\t{sample2}\t{poly_both_int0:,}: poly_both_int0')
                print(f'{sample1}\t{sample2}\t{inter_g1:,}: inter_g1')

                                # print(f'loc: {loc}, sample1: {sample1}, sample2: {sample2}, alleles_1: {[f"{al}:{positions[loc].freq_allele(sample1, al)}" for al in alleles_1]}, alleles_2: {[f"{al}:{positions[loc].freq_allele(sample2, al)}" for al in alleles_2]}, shared: {alleles_1.intersection(alleles_2)}')

# make fig 2b data for cumulative curves
'''
Per host, for each sample size, track the positions that are detected as polymorphic
and get the fraction that they make of the total core genome length
Then do this for increasing sample sizes using a random subset of samples
from 1 to 15 and for as many curves as needed
In the output have columns for
host curve_id nb_samples fraction_variable sdp
for each cutoff
'''

def coverage_in_sample(sample, sdp, limit = 20):
    with open(f'snvs/{sdp}_corecov_coord.txt', 'r') as cov_fh:
        for line in cov_fh:
            line = line.strip().split('\t')
            if line[1] == sample:
                if float(line[2]) >= limit:
                    return True

# read core lengths
core_lengths_df = pd.read_csv('04_RebuttalAnalyses/05_CoreGenomeLengthsInfo/core_lengths_final.txt', sep='\t', header=0)
n_samples = [x for x in all_samples if x[0] == 'N']
f_samples = [x for x in all_samples if x[0] == 'F']
nb_curves = 10
for sdp in sdp_list:
    print('Working on', sdp)
    if os.path.exists(f'{snvs_analysis_directory}/{sdp}_positions_dict.pkl'):
        positions = pickle.load(open(f'{snvs_analysis_directory}/{sdp}_positions_dict.pkl', 'rb'))
    # get the core length for the sdp
    print('Getting core length')
    core_length = core_lengths_df[core_lengths_df['Strain_id'] == sdp]['Tot_length(bp)'].values[0]
    df_cum_curve = pd.DataFrame(columns=['host', 'curve_id', 'nb_samples', 'fraction_variable', 'sdp', 'filt_type'])
    for host in ['N', 'F']:
        print('Working on host', host)
        if host == 'N':
            samples = [x for x in n_samples if coverage_in_sample(x, sdp)]
            # only consider samples that have at least 20x coverage of that SDP
        else:
            samples = [x for x in f_samples if coverage_in_sample(x, sdp)]
            # only consider samples that have at least 20x coverage of that SDP
        for cutoff in ['0', '0.01', '0.05', '0.1']:
            print('Working on cutoff', cutoff)
            for curve_id in range(nb_curves):
                print('Working on curve', curve_id)
                for nb_samples in range(1, len(samples)):
                    print('Working on nb_samples', nb_samples)
                    # get a random sample of nb_samples from samples
                    samples_subset = random.sample(samples, nb_samples)
                    # get the positions that are polymorphic in at least one of the samples
                    locs_var = set()
                    for sample in samples_subset:
                        for loc in positions:
                            if sample in positions[loc].samples:
                                if positions[loc].is_polymorphic(sample, cutoff=float(cutoff)):
                                    if loc not in locs_var:
                                        locs_var.add(loc)
                    # get the fraction of the core genome that is variable
                    fraction_variable = round(len(locs_var) / core_length * 100, 2)
                    df_cum_curve = df_cum_curve._append({'host': host, 'curve_id': curve_id, 'nb_samples': nb_samples, 'fraction_variable': fraction_variable, 'sdp': sdp, 'filt_type': cutoff}, ignore_index=True)
    # for each cutoff write the df as a file and exlude the cutoff column
    for cutoff in ['0', '0.01', '0.05', '0.1']:
        df_cum_curve_filt = df_cum_curve[df_cum_curve['filt_type'] == cutoff]
        df_cum_curve_filt = df_cum_curve_filt.drop(columns=['filt_type'])
        df_cum_curve_filt.to_csv(f'{snvs_analysis_directory}/figures/{sdp}_cum_curve_filt_{cutoff}.txt', sep='\t', index=False)

# Are there allles that are completely fixed in both samples but they are different in each.
for sdp in sdp_list:
        if os.path.exists(f'{snvs_analysis_directory}/{sdp}_positions_dict.pkl'):
            positions = pickle.load(open(f'{snvs_analysis_directory}/{sdp}_positions_dict.pkl', 'rb'))
        print('Working on', sdp)
        sample_pairs = list(combinations(all_samples, 2))
        sample_pairs = [x for x in sample_pairs if coverage_in_sample(x[0], sdp) and coverage_in_sample(x[1], sdp)]
        cutoff = 0.01
        with open(f'{snvs_analysis_directory}/{sdp}_are_fixed_allles_unique.txt', 'w') as fh_out:
            for sample1, sample2 in sample_pairs:
                fixed_and_shared = 0
                fixed_and_unique = 0
                total_seen = 0
                for loc in positions:
                    if sample1 in positions[loc].samples and sample2 in positions[loc].samples:
                        total_seen += 1
                        if positions[loc].is_polymorphic(sample1, cutoff=cutoff) or positions[loc].is_polymorphic(sample2, cutoff=cutoff):
                            continue
                        else:
                            alleles_1 = set(positions[loc].alleles_detected(sample1, cutoff=cutoff))
                            alleles_2 = set(positions[loc].alleles_detected(sample2, cutoff=cutoff))
                            if len(alleles_1.intersection(alleles_2)) == 0:
                                fixed_and_unique += 1
                            else:
                                fixed_and_shared += 1
                print(f'{sample1}\t{sample2}\t{total_seen:,}: total_seen')
                print(f'{sample1}\t{sample2}\t{fixed_and_shared:,}: fixed_and_shared')
                print(f'{sample1}\t{sample2}\t{fixed_and_unique:,}: fixed_and_unique')
                fh_out.write(f'{sample1}\t{sample2}\t{total_seen:,}: total_seen')
                fh_out.write(f'{sample1}\t{sample2}\t{fixed_and_shared:,}: fixed_and_shared')
                fh_out.write(f'{sample1}\t{sample2}\t{fixed_and_unique:,}: fixed_and_unique')