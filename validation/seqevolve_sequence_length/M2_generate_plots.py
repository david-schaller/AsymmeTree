# -*- coding: utf-8 -*-

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


df = pd.read_csv('results/all_distances.csv')

distances = sorted(df['d'].unique().tolist())

models = ['JC69', 'K80', 'WAG', 'JTT']

fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)
axs = axs.T 
axs = axs.flatten()
fig.set_size_inches(10,10)
fs = 14



for i in range(len(models)):
    
    ax_twin = axs[i].twinx()
    ax_twin.axhline(y=1, xmin=0, xmax=1,
                    color='grey', linewidth=0.7)
    
    model = models[i]
    
    df1 = df[(df['model'] == model)]
    
    df_na = df1.groupby(by=['model', 'd', 'length'], as_index=False).count()
    df_na['proportion_nan'] = 1 - df_na['d_hat'] / df_na['i']
    
    df1 = df1[(df1['d_hat'] < 6.0)]

    sns.boxplot(y='d_hat', x='d', 
                data=df1, 
                hue='length',
                ax=axs[i],
                fliersize=2, linewidth=0.7)
    
    axs[i].set_ylim(0, 7)
    axs[i].set_ylabel('estimated distance' if i in (0,1) else '', 
                      fontsize=fs)
    
    if i in (1,3):
        axs[i].set_xlabel('true distance', fontsize=fs)
    else:
        axs[i].set_xlabel('')
        
    if i == 0:
        axs[i].set_title('nucleotide models', fontsize=fs)
    elif i == 2:
        axs[i].set_title('amino acid models', fontsize=fs)
        
    axs[i].set_yticks([0, 1, 2, 3, 4, 5, 6])
    axs[i].tick_params(axis='both', which='major', labelsize=fs)
    
    xw = 1 / (len(distances))
    for j in range(len(distances)):
        axs[i].axhline(y=distances[j], xmin=(j)*xw, xmax=(j+1)*xw,
                       color='r', linewidth=0.7)
        
    axs[i].text(0.5, 0.95, model,
                horizontalalignment='center',
                verticalalignment='top',
                transform=axs[i].transAxes,
                fontsize=fs)
    if i == 0:
        handles, hlabels = axs[i].get_legend_handles_labels()
        hlabels = [f'{int(s):,}' for s in hlabels]
        axs[i].legend(handles, hlabels,
                      title='seq. length', loc='upper left', fontsize=fs,
                      title_fontsize=fs,)
    else:
        axs[i].legend().set_visible(False)
    
    sns.barplot(y='proportion_nan', x='d',
                data=df_na,
                hue='length',
                ax=ax_twin,
                palette='pastel')
    
    ax_twin.legend().set_visible(False)
    ax_twin.set_ylim(0, 7)
    ax_twin.set_ylabel('Proportion NaN ' \
                       if i in (2,3) else '', 
                       fontsize=fs, loc='bottom')
    ax_twin.set_yticks([0,1])
    ax_twin.tick_params(axis='both', which='major', labelsize=fs, colors='grey')
    ax_twin.spines['right'].set_color('grey')
    ax_twin.yaxis.label.set_color('grey')


plt.tight_layout()
plt.savefig("results/distance_boxplot.pdf", dpi=450)
