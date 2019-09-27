# -*- coding: utf-8 -*-
"""
Created on Wed Sep 25 13:57:00 2019

@author: karam
"""
import pandas as pd

filename = 'Metabolism_Genes_7v1.txt'
genedf = pd.read_table(filename)

actidict = {}

for index,row in genedf.iterrows():
    activities = row['Disease or Function']
    for activity in activities.split(','):
        if activity in list(actidict.keys()):
            for effect in row['Effect on Disease or Function'].split(','):
                if effect == 'increases':
                    actidict[activity][1] += 1
                if effect == 'decreases':
                    actidict[activity][0] += 1
        else:
            actidict[activity] = [0,0]
            for effect in row['Effect on Disease or Function'].split(','):
                if effect == 'increases':
                    actidict[activity][1] += 1
                if effect == 'decreases':
                    actidict[activity][0] += 1

acts = list(actidict.keys())
upvals = [value[0] for value in list(actidict.values())]
downvals = [value[1] for value in list(actidict.values())]

import numpy as np
import matplotlib.pyplot as plt

indexes = range(len(acts))
p1 = plt.bar(indexes, downvals, 0.35)
p2 = plt.bar(indexes, upvals, 0.35,
             bottom=downvals)

plt.ylabel('Number of Hits')
plt.title('Affected Pathways')
plt.xticks(indexes, acts, rotation = 90)
plt.yticks(np.arange(0, max(max(upvals),max(downvals)), 10))
plt.legend((p1[0], p2[0]), ('Down', 'Up'))
plt.tick_params(labelsize = 4)

plt.savefig('plot_of_{0}.jpg'.format(filename),bbox_inches='tight',dpi=600)