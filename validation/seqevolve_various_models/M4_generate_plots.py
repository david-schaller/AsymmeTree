# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt


df = pd.read_csv('results/all_distances.csv')


fig, axs = plt.subplots(2, 2, sharex=True, sharey=True)
axs = axs.flatten()
fig.set_size_inches(12,10)
fs = 12

models = df.columns[1:]
for i in range(len(models)):
    
    model = models[i]
    
    df_nonan = df.dropna(subset=[model])
    X = df_nonan['true'].values.reshape(-1, 1)
    Y = df_nonan[model].values.reshape(-1, 1)

    axs[i].scatter(df['true'], df[model], marker=".")
    axs[i].set_xlim(0.0, max(df['true']))
    axs[i].set_ylim(0.0, 5.5)
    axs[i].tick_params(axis='both', which='major', labelsize=fs)
    axs[i].tick_params(axis='both', which='major', labelsize=fs)
    if i % 2 == 0:
        axs[i].set_ylabel('estimated distance', fontsize=fs)
    if i > 1:
        axs[i].set_xlabel('true distance', fontsize=fs)
    
    lr = LinearRegression()
    lr.fit(X, Y)
    axs[i].text(0.05, 0.95, "{}\n\nslope {:.5}\nintercept {:.5}".format(model, lr.coef_[0, 0], lr.intercept_[0]),
                horizontalalignment='left',
                verticalalignment='top',
                transform=axs[i].transAxes,
                fontsize=fs)
    
    line_x = np.array([0.0, max(df['true'])]).reshape(-1, 1)
    line_y = lr.predict(line_x)
    axs[i].plot(line_x, line_y, 'darkred')

plt.tight_layout()
plt.savefig("results/distance_scatter.png", dpi=450)