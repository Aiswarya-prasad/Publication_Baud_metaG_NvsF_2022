import os
import sys
import glob
import pickle
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from scipy.stats import kruskal, wilcoxon, pearsonr

os.getcwd()
# set the right paths
my_figpath = 'figures/'

samples = ['N01','N02','N03','N04','N06','N07','N08','N09','N10','N11','N12','N13','N14','N15','N16', 'F01','F02','F03','F04','F06','F07','F08','F09','F10','F11','F12','F13','F14','F15','F16']
phylo_list = ["api","bapis","bifido","bom","com","firm4","firm5","fper","gilli","lkun","snod"]
treatments = ["Nurses","Foragers"]
sdp_list = ["bifido_1.1", "bifido_1.2","bifido_1.3","bifido_1.4",
                         "bifido_1.5","bifido_2","firm4_1","firm4_2","firm5_1",
                         "firm5_2","firm5_3","firm5_4","gilli_1","gilli_2",
                         "snod_1","bapis", "fper_1", "com_1"]
# read all files functional_gene_content/*_sigOG_loctags.txt
sigOGs = set()
for filepath in glob.glob('functional_gene_content/*_sigOG_loctags.txt'):
    if 'unassigned' in filepath:
        continue                                                                                                                                                                                                                                                   
    df_sigOGs_new = pd.read_csv(filepath, sep='\t', header=0)
    sigOGs.update(df_sigOGs_new.iloc[:,0])
# Read in the data
og_blasts = pd.read_csv('functional_gene_content/all_filt_ffn.blastn.1.filt.SDP.OG', sep='\t', header=None)
og_blasts.columns = ['orf', 'contig', 'gene', 'strain', 'phylo', 'sdp', 'OG']
dict_ogs = {}
dict_ogs = og_blasts.set_index('orf').to_dict()['OG']
set_signif_ogs = set([x for x in pd.read_csv('functional_gene_content/allsignifOGs.txt', header=None)[0]])

for sample in samples:
    df_coverages_dfs = pd.DataFrame()
    if os.path.exists(f'functional_gene_content/{sample}_orfcov.txt') == False:
        continue
    df_coverages_dfs = pd.read_csv(f'functional_gene_content/{sample}_orfcov.txt', sep='\t', header=0)
    # filter to keep only those records for whom the orfs that are there in og_blasts
    df_coverages_dfs = df_coverages_dfs[df_coverages_dfs['#rname'].isin(dict_ogs.keys())]
    df_coverages_dfs['length'] = df_coverages_dfs['endpos'] - df_coverages_dfs['startpos']
    df_coverages_dfs['sample'] = sample
    df_coverages_dfs['adj_depth'] = df_coverages_dfs['meandepth'] / df_coverages_dfs['length']
    # # add og information to column
    # df_coverages_dfs['OG'] = df_coverages_dfs['#rname'].apply(lambda x: dict_ogs[x])
    # # summarise in terms of length and coverage
    # df_coverages_dfs = df_coverages_dfs.groupby(['OG', 'sample']).agg({'length': 'sum', 'adj_depth': 'sum', 'meandepth': 'sum'}).reset_index()
    # add significance information
    # df_coverages_dfs['significance'] = df_coverages_dfs['OG'].apply(lambda x: 'significant' if x in set_signif_ogs else 'non-significant')
    df_coverages_dfs['significance'] = df_coverages_dfs['#rname'].apply(lambda x: 'significant' if dict_ogs[x] in set_signif_ogs else 'non-significant')
    # df_coverages_dfs_combined = pd.concat([df_coverages_dfs_combined, df_coverages_dfs])
    # if figure files exist, skip the loop!
    if not os.path.exists(my_figpath + f'adj_depth_vs_length5000_scatterplot_{sample}.png'):
        # scatter plot of coverage vs length for significant on top of non-significant ogs (blue) using plt
        plt.figure(figsize=(10, 6))  # Adjust figure size if needed
        plt.style.use('ggplot')
        plt.scatter(x=df_coverages_dfs[df_coverages_dfs['significance'] == 'non-significant']['length'],
                    y=df_coverages_dfs[df_coverages_dfs['significance'] == 'non-significant']['adj_depth'],
                    c='blue', label='non-significant', s=3, alpha=0.5)
        plt.scatter(x=df_coverages_dfs[df_coverages_dfs['significance'] == 'significant']['length'],
                    y=df_coverages_dfs[df_coverages_dfs['significance'] == 'significant']['adj_depth'],
                    c='red', label='significant', s=3, alpha=0.5)
        plt.xlim(0, 5000)
        plt.xlabel('Gene Length')
        plt.ylabel('Ratio of mean of coverage per base for the ORF and length')
        plt.title(f'Scatterplot of Coverage vs Length (sample: {sample})')
        plt.legend()
        plt.savefig(my_figpath + f'adj_depth_vs_length5000_scatterplot_{sample}.png', dpi=300)
    if not os.path.exists(my_figpath + f'adj_depth_vs_length_scatterplot_{sample}.png'):
        # scatter plot of coverage vs length for significant on top of non-significant ogs (blue) using plt
        plt.figure(figsize=(10, 6))  # Adjust figure size if needed
        plt.style.use('ggplot')
        plt.scatter(x=df_coverages_dfs[df_coverages_dfs['significance'] == 'non-significant']['length'],
                    y=df_coverages_dfs[df_coverages_dfs['significance'] == 'non-significant']['adj_depth'],
                    c='blue', label='non-significant', s=3, alpha=0.5)
        plt.scatter(x=df_coverages_dfs[df_coverages_dfs['significance'] == 'significant']['length'],
                    y=df_coverages_dfs[df_coverages_dfs['significance'] == 'significant']['adj_depth'],
                    c='red', label='significant', s=3, alpha=0.5)
        plt.xlabel('Gene Length')
        plt.ylabel('Ratio of mean of coverage per base for the ORF and length')
        plt.title(f'Scatterplot of Coverage vs Length (sample: {sample})')
        plt.legend()
        plt.savefig(my_figpath + f'adj_depth_vs_length_scatterplot_{sample}.png', dpi=300)
    if not os.path.exists(my_figpath + f'coverage_vs_length5000_scatterplot_{sample}.png'):
        # scatter plot of coverage vs length for significant on top of non-significant ogs (blue) using plt
        plt.figure(figsize=(10, 6))  # Adjust figure size if needed
        plt.style.use('ggplot')
        plt.scatter(x=df_coverages_dfs[df_coverages_dfs['significance'] == 'non-significant']['length'],
                    y=df_coverages_dfs[df_coverages_dfs['significance'] == 'non-significant']['meandepth'],
                    c='blue', label='non-significant', s=3, alpha=0.5)
        plt.scatter(x=df_coverages_dfs[df_coverages_dfs['significance'] == 'significant']['length'],
                    y=df_coverages_dfs[df_coverages_dfs['significance'] == 'significant']['meandepth'],
                    c='red', label='significant', s=3, alpha=0.5)
        # limit x axis to 5000
        plt.xlim(0, 5000)
        plt.xlabel('Gene Length')
        plt.ylabel('Mean of coverage per base for the ORF')
        plt.title(f'Scatterplot of Coverage vs Length (sample: {sample})')
        plt.legend()
        plt.savefig(my_figpath + f'coverage_vs_length5000_scatterplot_{sample}.png', dpi=300)
    if not os.path.exists(my_figpath + f'coverage_vs_length1000_scatterplot_{sample}.png'):
        # scatter plot of coverage vs length for significant on top of non-significant ogs (blue) using plt
        plt.figure(figsize=(10, 6))  # Adjust figure size if needed
        plt.style.use('ggplot')
        plt.scatter(x=df_coverages_dfs[df_coverages_dfs['significance'] == 'non-significant']['length'],
                    y=df_coverages_dfs[df_coverages_dfs['significance'] == 'non-significant']['meandepth'],
                    c='blue', label='non-significant', s=3, alpha=0.5)
        plt.scatter(x=df_coverages_dfs[df_coverages_dfs['significance'] == 'significant']['length'],
                    y=df_coverages_dfs[df_coverages_dfs['significance'] == 'significant']['meandepth'],
                    c='red', label='significant', s=3, alpha=0.5)
        # limit x axis to 1000
        plt.xlim(300, 1000)
        r, p = pearsonr(x=df_coverages_dfs[df_coverages_dfs['length'] <= 1000]['length'], 
                        y=df_coverages_dfs[df_coverages_dfs['length'] <= 1000]['meandepth'])
        pval = '< 0.05' if p < 0.05 else '> 0.05'
        plt.text(550, 200, f'Pearson\'s r = {r:.4f},\np = {pval}')
        plt.xlabel('Gene Length')
        plt.ylabel('Mean of coverage per base for the ORF')
        plt.title(f'Scatterplot of Coverage vs Length (sample: {sample})')
        plt.legend()
        plt.savefig(my_figpath + f'coverage_vs_length1000_scatterplot_{sample}.png', dpi=600)
    if not os.path.exists(my_figpath + f'coverage_vs_length_scatterplot_{sample}.png'):
        # scatter plot of coverage vs length for significant on top of non-significant ogs (blue) using plt
        plt.figure(figsize=(10, 6))  # Adjust figure size if needed
        plt.style.use('ggplot')
        plt.scatter(x=df_coverages_dfs[df_coverages_dfs['significance'] == 'non-significant']['length'],
                    y=df_coverages_dfs[df_coverages_dfs['significance'] == 'non-significant']['meandepth'],
                    c='blue', label='non-significant', s=3, alpha=0.5)
        plt.scatter(x=df_coverages_dfs[df_coverages_dfs['significance'] == 'significant']['length'],
                    y=df_coverages_dfs[df_coverages_dfs['significance'] == 'significant']['meandepth'],
                    c='red', label='significant', s=3, alpha=0.5)
        plt.xlabel('Gene Length')
        plt.ylabel('Mean of coverage per base for the ORF')
        plt.title(f'Scatterplot of Coverage vs Length (sample: {sample})')
        plt.legend()
        plt.savefig(my_figpath + f'coverage_vs_length_scatterplot_{sample}.png', dpi=300)
    # if not os.path.exists(my_figpath + f'coverage_vs_length_boxplot_{sample}.png'):
    #     # scatter plot of coverage vs length for significant on top of non-significant ogs (blue) using plt
    #     plt.figure(figsize=(10, 6))  # Adjust figure size if needed
    #     sig_lengths = df_coverages_dfs[df_coverages_dfs['significance'] == 'significant']['length']
    #     non_sig_lengths = df_coverages_dfs[df_coverages_dfs['significance'] == 'non-significant']['length']
    #     _, p_value = kruskal(sig_lengths, non_sig_lengths)
    #     p_string = '< 0.05' if p_value < 0.05 else '> 0.05'
    #     plt.figure(figsize=(10, 6))
    #     # plt.boxplot([sig_lengths, non_sig_lengths], labels=['significant', 'non-significant'])
    #     sns.boxplot(x='significance', y='length', data=df_coverages_dfs, showfliers = False)
    #     plt.xlabel('Significance')
    #     plt.ylabel('Gene Length')
    #     plt.title(f'Boxplot of Gene Length vs Significance (sample: {sample})')
    #     # Display p-vlue on the plot
    #     plt.annotate(f'p-value {p_string}', xy=(0.5, 0.9), xycoords='axes fraction', ha='center')
    #     plt.savefig(my_figpath + f'length_boxplot_{sample}.png', dpi=300)

df_coverages_dfs_combined = pd.DataFrame()
for sample in samples:
    df_coverages_dfs = pd.DataFrame()
    if os.path.exists(f'functional_gene_content/{sample}_orfcov.txt') == False:
        continue
    df_coverages_dfs = pd.read_csv(f'functional_gene_content/{sample}_orfcov.txt', sep='\t', header=0)
    # filter to keep only those records for whom the orfs that are there in og_blasts
    df_coverages_dfs = df_coverages_dfs[df_coverages_dfs['#rname'].isin(dict_ogs.keys())]
    df_coverages_dfs['length'] = df_coverages_dfs['endpos'] - df_coverages_dfs['startpos']
    df_coverages_dfs['sample'] = sample
    df_coverages_dfs['adj_depth'] = df_coverages_dfs['meandepth'] / df_coverages_dfs['length']
    # # add og information to column
    # df_coverages_dfs['OG'] = df_coverages_dfs['#rname'].apply(lambda x: dict_ogs[x])
    # # summarise in terms of length and coverage
    # df_coverages_dfs = df_coverages_dfs.groupby(['OG', 'sample']).agg({'length': 'sum', 'adj_depth': 'sum', 'meandepth': 'sum'}).reset_index()
    # add significance information
    # df_coverages_dfs['significance'] = df_coverages_dfs['OG'].apply(lambda x: 'significant' if x in set_signif_ogs else 'non-significant')
    df_coverages_dfs['significance'] = df_coverages_dfs['#rname'].apply(lambda x: 'significant' if dict_ogs[x] in set_signif_ogs else 'non-significant')
    df_coverages_dfs['host'] = df_coverages_dfs['sample'].apply(lambda x: 'Nurse' if x.startswith('N') else 'Forager')
    df_coverages_dfs_combined = pd.concat([df_coverages_dfs_combined, df_coverages_dfs])

# # pickle it
# df_coverages_dfs_combined.to_pickle(my_figpath + '/df_coverages_dfs_combined.pkl')
# # read it
# df_coverages_dfs_combined = pd.read_pickle(my_figpath + '/df_coverages_dfs_combined.pkl')

plt.figure(figsize=(10, 6))  # Adjust figure size if needed
plt.style.use('ggplot')
# change is in to startwith
plt.scatter(x=df_coverages_dfs_combined[(df_coverages_dfs_combined['significance'] == 'non-significant') & (df_coverages_dfs_combined['length'] <= 5000) & (df_coverages_dfs_combined['sample'].isin(['N01','N02','N03','N04','N06','N07','N08','N09','N10','N11','N12','N13','N14','N15','N16']))]['length'],
            y=df_coverages_dfs_combined[(df_coverages_dfs_combined['significance'] == 'non-significant') & (df_coverages_dfs_combined['length'] <= 5000) & (df_coverages_dfs_combined['sample'].isin(['N01','N02','N03','N04','N06','N07','N08','N09','N10','N11','N12','N13','N14','N15','N16']))]['adj_depth'],
            c='#fdbf6f', label='non-significant', s=3, alpha=0.5)
plt.scatter(x=df_coverages_dfs_combined[(df_coverages_dfs_combined['significance'] == 'significant') & (df_coverages_dfs_combined['length'] <= 5000) & (df_coverages_dfs_combined['sample'].isin(['N01','N02','N03','N04','N06','N07','N08','N09','N10','N11','N12','N13','N14','N15','N16']))]['length'],
            y=df_coverages_dfs_combined[(df_coverages_dfs_combined['significance'] == 'significant') & (df_coverages_dfs_combined['length'] <= 5000) & (df_coverages_dfs_combined['sample'].isin(['N01','N02','N03','N04','N06','N07','N08','N09','N10','N11','N12','N13','N14','N15','N16']))]['adj_depth'],
            c='#ff7f00', label='significant', s=3, alpha=0.5)
plt.scatter(x=df_coverages_dfs_combined[(df_coverages_dfs_combined['significance'] == 'non-significant') & (df_coverages_dfs_combined['length'] <= 5000) & (df_coverages_dfs_combined['sample'].isin(['F01','F02','F03','F04','F06','F07','F08','F09','F10','F11','F12','F13','F14','F15','F16']))]['length'],
            y=df_coverages_dfs_combined[(df_coverages_dfs_combined['significance'] == 'non-significant') & (df_coverages_dfs_combined['length'] <= 5000) & (df_coverages_dfs_combined['sample'].isin(['F01','F02','F03','F04','F06','F07','F08','F09','F10','F11','F12','F13','F14','F15','F16']))]['adj_depth'],
            c='#a6cee3', label='non-significant', s=3, alpha=0.5)
plt.scatter(x=df_coverages_dfs_combined[(df_coverages_dfs_combined['significance'] == 'significant') & (df_coverages_dfs_combined['length'] <= 5000) & (df_coverages_dfs_combined['sample'].isin(['F01','F02','F03','F04','F06','F07','F08','F09','F10','F11','F12','F13','F14','F15','F16']))]['length'],
            y=df_coverages_dfs_combined[(df_coverages_dfs_combined['significance'] == 'significant') & (df_coverages_dfs_combined['length'] <= 5000) & (df_coverages_dfs_combined['sample'].isin(['F01','F02','F03','F04','F06','F07','F08','F09','F10','F11','F12','F13','F14','F15','F16']))]['adj_depth'],
            c='#1f78b4', label='non-significant', s=3, alpha=0.5)
plt.xlabel('Gene Length')
plt.ylabel('Ratio of mean of coverage per base for the ORF and length')
plt.title(f'Scatterplot of Coverage vs Length (sample: {sample})')
plt.legend()
plt.savefig(my_figpath + f'adj_depth_vs_length5000_scatterplot_nurse_forager.png', dpi=300)
# plt.savefig(my_figpath + f'adj_depth_vs_length_scatterplot_nurse_forager.png', dpi=300)

df_coverages_dfs_combined_info = df_coverages_dfs_combined.copy()
# add host column showimng whether the sample is nurse or forager
df_coverages_dfs_combined_info['host'] = df_coverages_dfs_combined_info['sample'].apply(lambda x: 'Nurse' if x.startswith('N') else 'Forager')
# scatter plot of coverage vs length for significant on top of non-significant ogs (blue) using plt higlighting also whether the sample is nurse or forager binary
# limit to length 5000 and facet sig and non-sig
plt.figure(figsize=(10, 6))  # Adjust figure size if needed
plt.style.use('ggplot')
g = sns.FacetGrid(df_coverages_dfs_combined_info[df_coverages_dfs_combined_info['length'] <= 5000], col='host', row='significance', hue='host', height=6)
g.map(plt.scatter, 'length', 'adj_depth', alpha=0.3, s=3)
g.add_legend()
plt.savefig(my_figpath + f'adj_depth_vs_length5000_scatterplot_facet.png', dpi=300)

plt.figure(figsize=(10, 6))  # Adjust figure size if needed
plt.style.use('ggplot')
g = sns.FacetGrid(df_coverages_dfs_combined_info[df_coverages_dfs_combined_info['length'] <= 5000], col='host', row='significance', hue='host', height=6)
g.map(plt.scatter, 'length', 'meandepth', alpha=0.3, s=3)
g.add_legend()
plt.savefig(my_figpath + f'coverage_vs_length5000_scatterplot_facet.png', dpi=300)

# plot a histogram of gene lengths of significant and non-significant ogs on top of each other
plt.figure(figsize=(10, 6))  # Adjust figure size if needed
plt.style.use('ggplot')
plt.hist(df_coverages_dfs_combined[df_coverages_dfs_combined['significance'] == 'non-significant']['length'], bins=100, alpha=0.5, label='non-significant', color='blue')
plt.hist(df_coverages_dfs_combined[df_coverages_dfs_combined['significance'] == 'significant']['length'], bins=100, alpha=0.5, label='significant', color='red')
plt.xlabel('Gene Length')
plt.ylabel('Frequency')
plt.xscale("log")
plt.title(f'Histogram of Gene Lengths')
plt.legend()
plt.savefig(my_figpath + f'length_histogram.png', dpi=300)

plt.figure(figsize=(10, 6))  # Adjust figure size if needed
plt.style.use('ggplot')
sns.displot(df_coverages_dfs_combined[df_coverages_dfs_combined['length'] <= 5000], x = 'length', y = 'coverage', hue = 'significance', kind = 'kde', fill = True)
plt.xlabel('Gene Length')
plt.ylabel('Frequency')
plt.xscale("log")
plt.title(f'Histogram of Gene Lengths')
plt.legend()
plt.savefig(my_figpath + f'length_coverage_kde.png', dpi=300)

plt.figure(figsize=(10, 6))  # Adjust figure size if needed
plt.style.use('ggplot')
plt.hist(df_coverages_dfs_combined[df_coverages_dfs_combined['significance'] == 'significant']['length'], bins=100, alpha=0.5, label='significant', color='red')
plt.xlabel('Gene Length')
plt.ylabel('Frequency')
plt.title(f'Histogram of Gene Lengths')
plt.legend()
plt.savefig(my_figpath + f'length_histogram_sig.png', dpi=300)

plt.figure(figsize=(10, 6))  # Adjust figure size if needed
plt.style.use('ggplot')
plt.hist(df_coverages_dfs_combined[df_coverages_dfs_combined['significance'] == 'non-significant']['length'], bins=100, alpha=0.5, label='non-significant', color='blue')
plt.xlabel('Gene Length')
plt.ylabel('Frequency')
plt.title(f'Histogram of Gene Lengths')
plt.legend()
plt.savefig(my_figpath + f'length_histogram_nonsig.png', dpi=300)

# boxplot of genelengths for significant and non-significant ogs colored by host
# do a wilcox and add the point behind the boxplot
plt.figure(figsize=(10, 6))  # Adjust figure size if needed
plt.style.use('ggplot')
plt.scatterplot(x='significance', y='length', hue='host', data=df_coverages_dfs_combined.sample(10000), alpha=0.5, s=3, x_jitter=0.2)
sns.boxplot(x='significance', y='length', hue='significance', data=df_coverages_dfs_combined, showfliers = False)
k_res = kruskal(df_coverages_dfs_combined[df_coverages_dfs_combined['significance'] == 'significant']['length'], df_coverages_dfs_combined[df_coverages_dfs_combined['significance'] == 'non-significant']['length'])
plt.annotate(f'p-value = {k_res.pvalue}', xy=(0.5, 0.9), xycoords='axes fraction', ha='center')
plt.xlabel('Significance')
plt.ylabel('Gene Length')
plt.title(f'Boxplot of Gene Length vs Significance')
plt.savefig(my_figpath + f'length_boxplot.png', dpi=300)

# plot all the significant ones next to non-significant in a violin plot
plt.figure(figsize=(10, 6))  # Adjust figure size if needed
plt.style.use('ggplot')
sns.violinplot(x='significance', y='length', hue='significance', data=df_coverages_dfs_combined[df_coverages_dfs_combined['length'] <= 100000], inner=None, color='lightgray')
# sns.stripplot(x='significance', y='length', data=df_coverages_dfs_combined[df_coverages_dfs_combined['length'] <= 100000], jitter=True, size=1.5)
plt.xlabel('Significance')
plt.ylabel('Gene Length')
plt.title(f'Violinplot of Gene Length vs Significance')
plt.savefig(my_figpath + f'length_violinplot.png', dpi=300)

plt.figure(figsize=(10, 6))  # Adjust figure size if needed
plt.style.use('ggplot')
sns.violinplot(x='significance', y='length', hue='significance', data=df_coverages_dfs_combined[df_coverages_dfs_combined['length'] <= 5000], inner=None, color='lightgray')
# sns.stripplot(x='significance', y='length', data=df_coverages_dfs_combined[df_coverages_dfs_combined['length'] <= 100000], jitter=True, size=1.5)
plt.xlabel('Significance')
plt.ylabel('Gene Length')
plt.title(f'Violinplot of Gene Length vs Significance')
plt.savefig(my_figpath + f'length_violinplot5000.png', dpi=300)

with open(my_figpath + 'pearson_values_300_to_1000.txt', 'w') as f:
    f.write(f'sample\tr\tp\n')
    for sample in samples:
        df_coverages_dfs = pd.DataFrame()
        if os.path.exists(f'functional_gene_content/{sample}_orfcov.txt') == False:
            continue
        df_coverages_dfs = pd.read_csv(f'functional_gene_content/{sample}_orfcov.txt', sep='\t', header=0)
        # filter to keep only those records for whom the orfs that are there in og_blasts
        df_coverages_dfs = df_coverages_dfs[df_coverages_dfs['#rname'].isin(dict_ogs.keys())]
        df_coverages_dfs['length'] = df_coverages_dfs['endpos'] - df_coverages_dfs['startpos']
        df_coverages_dfs['sample'] = sample
        df_coverages_dfs['adj_depth'] = df_coverages_dfs['meandepth'] / df_coverages_dfs['length']
        # # add og information to column
        # df_coverages_dfs['OG'] = df_coverages_dfs['#rname'].apply(lambda x: dict_ogs[x])
        # # summarise in terms of length and coverage
        # df_coverages_dfs = df_coverages_dfs.groupby(['OG', 'sample']).agg({'length': 'sum', 'adj_depth': 'sum', 'meandepth': 'sum'}).reset_index()
        # add significance information
        # df_coverages_dfs['significance'] = df_coverages_dfs['OG'].apply(lambda x: 'significant' if x in set_signif_ogs else 'non-significant')
        df_coverages_dfs['significance'] = df_coverages_dfs['#rname'].apply(lambda x: 'significant' if dict_ogs[x] in set_signif_ogs else 'non-significant')
        # df_coverages_dfs_combined = pd.concat([df_coverages_dfs_combined, df_coverages_dfs])
        # if figure files exist, skip the loop!
        r, p = pearsonr(x=df_coverages_dfs[df_coverages_dfs['length'] <= 1000]['length'], 
                        y=df_coverages_dfs[df_coverages_dfs['length'] <= 1000]['meandepth'])
        f.write(f'{sample}\t{r}\t{p}\n')

for sample in ['N03', 'F07']:
    df_coverages_dfs = pd.DataFrame()
    if os.path.exists(f'functional_gene_content/{sample}_orfcov.txt') == False:
        continue
    df_coverages_dfs = pd.read_csv(f'functional_gene_content/{sample}_orfcov.txt', sep='\t', header=0)
    # filter to keep only those records for whom the orfs that are there in og_blasts
    df_coverages_dfs = df_coverages_dfs[df_coverages_dfs['#rname'].isin(dict_ogs.keys())]
    df_coverages_dfs['length'] = df_coverages_dfs['endpos'] - df_coverages_dfs['startpos']
    df_coverages_dfs['sample'] = sample
    df_coverages_dfs['adj_depth'] = df_coverages_dfs['meandepth'] / df_coverages_dfs['length']
    # # add og information to column
    # df_coverages_dfs['OG'] = df_coverages_dfs['#rname'].apply(lambda x: dict_ogs[x])
    # # summarise in terms of length and coverage
    # df_coverages_dfs = df_coverages_dfs.groupby(['OG', 'sample']).agg({'length': 'sum', 'adj_depth': 'sum', 'meandepth': 'sum'}).reset_index()
    # add significance information
    # df_coverages_dfs['significance'] = df_coverages_dfs['OG'].apply(lambda x: 'significant' if x in set_signif_ogs else 'non-significant')
    df_coverages_dfs['significance'] = df_coverages_dfs['#rname'].apply(lambda x: 'significant' if dict_ogs[x] in set_signif_ogs else 'non-significant')    
    # scatter plot of coverage vs length for significant on top of non-significant ogs (blue) using plt
    plt.figure(figsize=(10, 6))  # Adjust figure size if needed
    plt.style.use('ggplot')
    plt.scatter(x=df_coverages_dfs[df_coverages_dfs['significance'] == 'non-significant']['length'],
                y=df_coverages_dfs[df_coverages_dfs['significance'] == 'non-significant']['meandepth'],
                c='blue', label='non-significant', s=3, alpha=0.5)
    plt.scatter(x=df_coverages_dfs[df_coverages_dfs['significance'] == 'significant']['length'],
                y=df_coverages_dfs[df_coverages_dfs['significance'] == 'significant']['meandepth'],
                c='red', label='significant', s=3, alpha=0.5)
    # limit x axis to 1000
    plt.xlim(300, 1000)
    plt.ylim(0, 500)
    r, p = pearsonr(x=df_coverages_dfs[df_coverages_dfs['length'] <= 1000]['length'], 
                    y=df_coverages_dfs[df_coverages_dfs['length'] <= 1000]['meandepth'])
    pval = '< 0.05' if p < 0.05 else '> 0.05'
    plt.text(550, 200, f'Pearson\'s r = {r:.4f},\np = {pval}')
    plt.xlabel('Gene Length')
    plt.ylabel('Mean of coverage per base for the ORF')
    plt.title(f'Scatterplot of Coverage vs Length (sample: {sample})')
    plt.legend()
    plt.savefig('figures' + f'coverage_vs_length1000_scatterplot_{sample}.png', dpi=600)