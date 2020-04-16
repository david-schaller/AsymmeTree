# -*- coding: utf-8 -*-

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv('results/all_distances.csv')

distances = sorted(df['d'].unique().tolist())

models = ['JC69', 'K80', 'WAG', 'JTT']

fig, axs = plt.subplots(4, 1, sharex=True, sharey=True)
axs = axs.flatten()
fig.set_size_inches(8,12)
fs = 12



for i in range(len(models)):
    
    model = models[i]
    
    df1 = df[(df['model'] == model) & (df['d_hat'] < 6.0)]

    sns.boxplot(y='d_hat', x='d', 
                data=df1, 
                hue='length',
                ax=axs[i],
                fliersize=2, linewidth=0.7)
    
    axs[i].set_ylabel('estimated distance', fontsize=fs)
    if i > 2:
        axs[i].set_xlabel('true distance', fontsize=fs)
    else:
        axs[i].set_xlabel('')
    
    xw = 1 / (len(distances))
    for j in range(len(distances)):
        axs[i].axhline(y=distances[j], xmin=(j)*xw, xmax=(j+1)*xw,
                       color='r', linewidth=0.7)
        
    axs[i].text(0.95, 0.05, model,
                horizontalalignment='right',
                verticalalignment='bottom',
                transform=axs[i].transAxes,
                fontsize=fs)


plt.tight_layout()
plt.savefig("results/distance_boxplot.pdf", dpi=450)