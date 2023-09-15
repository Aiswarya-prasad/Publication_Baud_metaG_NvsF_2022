import os
import pickle
import sys
import glob
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.ticker import FuncFormatter


'''
Please specify for each species/phylotype how large the reference set of core genes was 
that was used for strain calling. How long was the concatenated length of all used 
genes in each case? What fraction of the ref genome did that correspond to?

Make a summary table for this information.
'''

os.getcwd()
beebiome_db_path = 'path/to/beebiome/db'
figpath = 'figures'
genome_path = f'{beebiome_db_path}/fna_files/'
workingdir = 'core_genome_lengths_summary'
# os.system(f'bash calc_core_lengths.sh path/to/bed_filt_files > {workingdir}/core_lengths_reduced_filt_bed.txt')
df_core = pd.read_csv(f'{workingdir}/core_lengths_reduced_filt_bed.txt', sep='\t')
# os.system(f'bash calc_core_lengths.sh path/to/bed_red_files > {workingdir}/table_core_lengths_raw.txt')
df_core_values = pd.read_csv(f'{workingdir}/table_core_length_raw.txt', sep='\t')
# replace df_core Tot_length(bp) with df_core_values Tot_length(bp) by matching new_id
# becauase the original df_core Tot_length(bp) was calculated from bed file with low coverage genes removed!
# append the column name of Tot_length(bp) from df_core_values with _filt
df_core = df_core.rename(columns={'Tot_length(bp)': 'Tot_length_filt(bp)'})
df_core = df_core.merge(df_core_values[['new_id', 'Tot_length(bp)']], on='new_id', how='left')
# drop bifido_1_cerana
df_core = df_core[df_core['Strain_id'] != 'bifido_1_cerana']
genomes = df_core['new_id'].values
for genome in genomes:
    genome_file = os.path.join(genome_path, f'{genome}.fna')
    genome_length = 0
    for record in SeqIO.parse(genome_file, 'fasta'):
        genome_length += len(record.seq)
    df_core.loc[df_core['new_id'] == genome, 'genome_length'] = genome_length
df_core['core_genome_percent'] = df_core['Tot_length(bp)']/df_core['genome_length']*100
df_core['core_genome_filt_percent'] = df_core['Tot_length_filt(bp)']/df_core['genome_length']*100
df_core.to_csv(f'{workingdir}/core_lengths_final.txt', sep='\t', index=False)
# plot a bar graph showing genome length and core genome length as stacked bars
# format y axis to show in millions
def millions_formatter(x, pos):
    return f'{x / 1000000}M'

# only keep sdps of interest
sdp_list = ["bapis","bifido_1.1","bifido_1.2","bifido_1.3","bifido_1.4",
                         "bifido_1.5","bifido_2","firm4_1","firm4_2","firm5_1",
                         "firm5_2","firm5_3","firm5_4","fper_1","gilli_1",
                         "gilli_2", "gilli_3","snod_1"]
plt.figure(figsize=(10, 10))
plt.style.use('ggplot')
sns.set_context('paper', font_scale=1.5)
# make stacked bar plot of genome length and core genome length and core genome length filtered
df_core_plot = df_core
# df_core_plot = df_core[df_core['Strain_id'].isin(sdp_list)]


ax = sns.barplot(y='Strain_id', x='genome_length', data=df_core_plot, color='#1f78b4')
sns.barplot(y='Strain_id', x='Tot_length(bp)', data=df_core_plot, color='#b2df8a')
sns.barplot(y='Strain_id', x='Tot_length_filt(bp)', data=df_core_plot, color='#33a02c')
sns.set(rc={'axes.facecolor':'#f2f2f2', 'figure.facecolor':'white'})
plt.ylabel('Species')
plt.xlabel('Genome length (bp)')
# add legend for color
colors = {'Core genome length': '#b2df8a', 'Core genome length filtered': '#33a02c', 'Total genome length': '#1f78b4'}
labels = list(colors.keys())
handles = [plt.Rectangle((0,0),1,1, color=colors[label]) for label in labels]
plt.legend(handles, labels, loc='right')
ax.xaxis.set_major_formatter(FuncFormatter(millions_formatter))
plt.savefig(f'{figpath}/SuppFig10-pre-CoreGenomeLengths.pdf', format='pdf', bbox_inches='tight')
