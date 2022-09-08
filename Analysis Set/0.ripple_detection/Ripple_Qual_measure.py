# -*- coding: utf-8 -*-
"""
Created on Thu Sep  8 10:32:58 2022

@author: Jaemin
"""

import pandas as pd
import matplotlib.pyplot as plt


R_dir = 'D:/HPC-SWR project/Processed Data/ripples_mat/R2_300'

R_df = pd.read_excel(f'{R_dir}/RipplesTable_Behav_CA1_speed_filtered.xlsx')
R_df.set_index(keys=['ID'],inplace=True,drop=True)

R_df_t = R_df[(R_df['ensemble']>2) & (R_df['manual']>0)]

plt.style.use('default')
plt.rcParams['figure.figsize'] = (8,6)
plt.rcParams['font.size']=12

fig,ax = plt.subplots()

ax.violinplot([R_df_t.RippleDuration[R_df_t['experimenter']=='LSM'],
            R_df_t.RippleDuration[R_df_t['experimenter']=='SEB'],
            R_df_t.RippleDuration[R_df_t['experimenter']=='JS']],
           showmeans=True)
ax.set_xlabel('Dataset')
# ax.set_xticklabel(['LSM','Delcasso','Jhoseph'])

plt.show()


fig,ax = plt.subplots()

ax.boxplot([R_df_t.MeanFreq[R_df_t['experimenter']=='LSM'],
            R_df_t.MeanFreq[R_df_t['experimenter']=='SEB'],
            R_df_t.MeanFreq[R_df_t['experimenter']=='JS']],
           notch=True,whis=2.5)
ax.set_xlabel('Dataset')
# ax.set_xticklabel(['LSM','Delcasso','Jhoseph'])

plt.show()
