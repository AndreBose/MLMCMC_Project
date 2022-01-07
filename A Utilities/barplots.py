#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 14:18:46 2021

@author: andreaboselli
"""

# %% Import libraries

import numpy             as np
import matplotlib.pyplot as plt



# %% Set tables to be converted in barplots

#                          \mu,      \theta
ESS_table     = np.array([[483.3679 ,701.9119 ],  # MH
                          [367.3521 ,384.8069 ],  # DEMZ
                          [1107.439 ,938.6840 ],  # MLDA_without_VR_nsubs_5
                          [1471.502 ,1511.576 ]]) # MLDA_without_VR_nsubs_20

ESS_sec_table = np.array([[0.104096 ,0.151161 ],  # MH
                          [0.075778 ,0.079379 ],  # DEMZ
                          [0.270917 ,0.229634 ],  # MLDA_without_VR_nsubs_5
                          [0.283208 ,0.290921 ]]) # MLDA_without_VR_nsubs_20

variables_names =         ['mu'     ,'theta'  ]

methods_names   = ['Metropolis',
                   'DEMetropolisZ',
                   'MLDA without VR (nsubs = 5)',
                   'MLDA without VR (nsubs = 20)',]

methods_colors  = ['deepskyblue',
                   'green',
                   'lightcoral',
                   'brown']



# %% Convert ESS_table into a barplot

nvars    = len(variables_names)
nmethods = len(methods_names)

width = 0.10                                    # the width of the bars
x = np.arange(nvars) * (width * nmethods * 1.7) # the labels locations


fig1, ax1 = plt.subplots(dpi=300)

for i in range(nmethods):
    
    ax1.bar(x      = x + i*width - (nmethods-1)*width/2, 
            height = ESS_table[i], 
            width  = width, 
            color  = methods_colors[i],
            label  = methods_names [i])


plt.xticks(x, variables_names)
plt.ylim(0,np.max(ESS_table)*1.5)
ax1.set_ylabel('ESS')
ax1.set_title ('ESS')


ax1.legend(framealpha=1, loc='best')
fig1.tight_layout()
ax1.grid(axis='y', linewidth=0.3)

plt.show()



# %% Convert ESS_sec_table into a barplot


fig2, ax2 = plt.subplots(dpi=300)

for i in range(nmethods):
    
    ax2.bar(x      = x + i*width - (nmethods-1)*width/2, 
            height = ESS_sec_table[i], 
            width  = width, 
            color  = methods_colors[i],
            label  = methods_names [i])


plt.xticks(x, variables_names)
plt.ylim(0,np.max(ESS_sec_table)*1.5)
ax2.set_ylabel('ESS/sec')
ax2.set_title ('ESS/sec')


ax2.legend(framealpha=1, loc='best')
fig2.tight_layout()
ax2.grid(axis='y', linewidth=0.3)

plt.show()