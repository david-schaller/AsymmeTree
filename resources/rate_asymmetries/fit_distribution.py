#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import csv
import numpy as np
from scipy import stats  

import matplotlib
matplotlib.rcParams['text.usetex'] = True
matplotlib.rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern Sans serif']})
matplotlib.rc('text.latex', preamble=r'\usepackage{cmbright}')
import matplotlib.pyplot as plt

values = []

with open('S_cerevisiae_ByrneWolfe2007.csv', newline='') as csvfile:
    reader = csv.reader(csvfile, delimiter=',', quotechar='|')
    for row in reader:
        values.append(float(row[1]))

values = np.asarray(values)

print(values)

fig, ax = plt.subplots(1, 1, sharex=False, sharey=False)
fig.set_size_inches(10,5)
fs = 18

# plot normed histogram
plt.hist(values, density=True, bins = 80, color="lightgrey", rwidth=1.01)

xt = plt.xticks()[0]  
xmin, xmax = min(xt), max(xt)  
lnspc = np.linspace(xmin, xmax, len(values))


a, loc, scale = stats.gamma.fit(values)  
pdf_gamma = stats.gamma.pdf(lnspc, a, loc, scale)
ax.plot(lnspc, pdf_gamma, label="fitted Gamma distribution\n(shape={:.3}, scale={:.4}, loc={:.3})".format(a, scale, loc))

print("shape={:.3}, scale={:.4}, loc={:.3}".format(a, scale, loc))

ax.set_xlabel(r"$R'$-value", fontsize=fs)
ax.set_ylabel(r"Frequency / density", fontsize=fs)
ax.set_xticks([i for i in range(0,17)])
ax.set_xlim(0,16)
ax.tick_params(axis='both', labelsize=fs-2)

plt.legend(fontsize=fs, labelspacing=1.0)
plt.tight_layout()
#plt.savefig('fitted_distribution.pdf')
plt.show()  