# -*- coding: utf-8 -*-

import os
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt


directory = 'results'

death_rates = [0.0, 0.2, 0.5, 0.8]

R_ages, py_ages = [], []

for tool, ages in zip(('TreeSim', 'py'), (R_ages, py_ages)):
    
    for d in death_rates:
        ages.append([])
        filename = os.path.join(directory, f'{tool}_age{d}.txt')
        with open(filename, 'r') as f:
            for line in f.readlines():
                line = line.strip()
                if line:
                    ages[-1].append(float(line))
    
        print(len(ages[-1]))


fig, axs = plt.subplots(1, 4, sharex=False, sharey=True)
fig.set_size_inches(14,5)
bins = 20
colors = ['b', 'darkorange']
fs = 14

for i, d in enumerate(death_rates):
    
    mean1 = str(round(np.mean(R_ages[i]),2))
    mean2 = str(round(np.mean(py_ages[i]),2))

    axs[i].hist([R_ages[i], py_ages[i]], label=[f'TreeSim\n(mean {mean1})',
                                                f'AsymmeTree\n(mean {mean2})'],
                 bins=bins, alpha=0.7, color=colors, density=False)
    
    pval = stats.mannwhitneyu(R_ages[i], py_ages[i], alternative="two-sided")
    print(pval)
    
    axs[i].text(0.9, 0.65, 'p-val. ' + '{:.3g}'.format(pval[1]),
                    horizontalalignment='right', verticalalignment='top',
                    transform=axs[i].transAxes, fontsize=fs)
    axs[i].set_title(f'spec. rate 1.0, ext. rate {d}', fontsize=fs)
    axs[i].set_xlabel('tree age', fontsize=fs)
    axs[i].tick_params(axis='both', which='major', labelsize=fs)
    axs[i].legend(fontsize=fs)

plt.tight_layout()
plt.savefig(os.path.join(directory, 'species_age_histo.pdf'))
