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
                          [465.1011 ,579.1599 ]]) # MLDA_without_VR

ESS_sec_table = np.array([[0.118379 ,0.171902 ],  # MH
                          [0.059717 ,0.074362 ]]) # MLDA_without_VR

variables_names =         ['mu'     ,'theta'  ]

methods_names   = ['Metropolis',
                   'MLDA without VR']

methods_colors  = ['deepskyblue',
                   'red']



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
ax1.set_ylabel('ESS')
ax1.set_title ('ESS')


ax1.legend(framealpha=1)
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
ax2.set_ylabel('ESS/sec')
ax2.set_title ('ESS/sec')


ax2.legend(framealpha=1)
fig2.tight_layout()
ax2.grid(axis='y', linewidth=0.3)

plt.show()