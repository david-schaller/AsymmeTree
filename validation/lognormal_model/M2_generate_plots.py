# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms


df = pd.read_csv('testfile_scenarios.csv')

variances = sorted(df['variance'].unique())

df['lograte'] = np.log2(df['rate'])

fig, axs = plt.subplots(2, 3, sharex=False, sharey=False)
axs = axs.flatten()
fig.set_size_inches(16,7)
fs = 16

for i, var in enumerate(variances):
    
    df_var = df[df['variance'] == var]
    
    mean = np.mean(df_var['rate'])
    
    axs[i].hist(df_var['rate'], bins=100, density=True, alpha=0.5)
    axs[i].hist(df_var['lograte'], bins=100, density=True, alpha=0.5)
    
    axs[i].tick_params(axis='both', labelsize=fs)
    
    axs[i].axvline(mean, color='r')
    
    # axs[i].set_title(f'β = {var}', fontsize=fs)
    axs[i].text(0.05, 0.95, f'β = {var}',
                horizontalalignment='left',
                verticalalignment='top',
                transform=axs[i].transAxes,
                fontsize=fs)
    
    if i >= 2:
        axs[i].text(mean + 0.05 + i * 0.03, 0.95, f'mean = {mean:.3f}',
                    horizontalalignment='left',
                    verticalalignment='top',
                    transform=transforms.blended_transform_factory(
                        axs[i].transData, axs[i].transAxes),
                    fontsize=fs)
    else:
        axs[i].text(mean - 0.07, 0.95, f'mean = {mean:.3f}',
                    horizontalalignment='right',
                    verticalalignment='top',
                    transform=transforms.blended_transform_factory(
                        axs[i].transData, axs[i].transAxes),
                    fontsize=fs)


plt.tight_layout()
plt.savefig('lognormal_bins.pdf')
