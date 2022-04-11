#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


__author__ = 'David Schaller'


plot_dir = 'results'
if not os.path.exists(plot_dir):
    os.makedirs(plot_dir)
    
quotient = True

filename = 'boxplots.pdf'

df = pd.read_csv('results/distances.csv')
df.sort_values('sim_ID', inplace=True)
df.set_index('sim_ID', inplace=True)
df.replace({'r': 'replacing', 'a': 'additive', 'False': 'no bias'}, inplace=True)


fig, axs = plt.subplots(1, 1, sharex=True, sharey=False,)
                        #gridspec_kw={'height_ratios': [2, 1]})
axs = [axs]
fig.set_size_inches(10,7)
fs = 18

# customized boxplot properties
bp_props = {'boxprops': dict(linewidth=1, edgecolor='black'),
            'flierprops': dict(markerfacecolor='black',
                                markeredgecolor='black',
                                markersize=1),
            'medianprops': dict(linewidth=1.5,
                                color='black',
                                solid_capstyle='butt'),
            'whiskerprops': dict(linewidth=1, color='black',
                                  solid_capstyle='butt'),
            'capprops': dict(linewidth=1, color='black')
            }

sns.boxplot(x='bias_mode', y='distance', hue='type', data=df, ax=axs[0],
            palette=sns.color_palette("Set2"),
            **bp_props)



axs[0].set_xlabel('transfer distance bias', fontsize=fs)
axs[0].set_ylabel('elapsed time since divergence', fontsize=fs)
axs[0].legend(title='', fontsize=fs-4, loc='upper right',
              fancybox=True, shadow=True, )
axs[0].tick_params(axis='both', which='major', labelsize=fs-2)
    

plt.tight_layout()
plt.savefig( os.path.join(plot_dir, filename) )