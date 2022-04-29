# -*- coding: utf-8 -*-

import itertools

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

tools = ['seqgen', 'indelible', 'alf', 'pyvolve', 'asymmetree']
tree_sizes = [20, 50, 100]
lengths = [100, 500, 2000]

df = pd.read_csv('results/times_intel.csv')

fig, axs = plt.subplots(2, 3, sharex=True, sharey=True)
fig.set_size_inches(10,8)
fs = 12

for i, j in itertools.product(range(2), range(3)):
    
    df_filtered = df[(df['mode'] == i+1) & (df['N'] == tree_sizes[j])]
    
    axs[i][j].set_yscale('log')
    
    for t, tool in enumerate(tools):
        medians = []
        
        for l in lengths:
            df_filtered2 = df_filtered[(df_filtered['tool'] == tool) &
                                       (df_filtered['l'] == l)]
            medians.append(df_filtered2['time'].median())
        
        print(tool, medians)
        axs[i][j].plot(np.arange(-0.3 + 0.15 * t, 5.0, 1.0)[:3],
                       medians, zorder=1)

    sns.boxplot(y='time', x='l', 
                data=df_filtered, 
                hue='tool', hue_order=tools,
                ax=axs[i][j],
                fliersize=2, linewidth=0.7, zorder=3)
    
    for k in range(2):
        axs[i][j].axvline(k+0.5, linewidth=0.75, color='gray')
    
    axs[i][j].set_xlabel('seq. length' if i == 1 else '', fontsize=fs)
    axs[i][j].set_ylabel('time [s]' if j == 0 else '', fontsize=fs)
    axs[i][j].tick_params(axis='both', which='major', labelsize=fs)
    
    if i == 0:
        axs[i][j].set_title(f'no. of leaves: {tree_sizes[j]}')
    
    if not (i == 0 and j == 0):
        axs[i][j].legend().set_visible(False)


plt.tight_layout()
plt.savefig('results/runtime_boxplot.pdf')