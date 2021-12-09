# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 10:43:10 2021

1. Reactivated group의 RDI 분포 확인
@author: user
"""

#%% import libraries
import numpy as np 
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

#%% import data
ROOT_data = 'D:/HPC-SWR project/Processed Data'

thisRID = 561
thisSID = 1
thisRegion = 'CA1'

df_rip = pd.read_excel(f'{ROOT_data}/RipplesTable_r{thisRID}_s{str(thisSID).zfill(2)}_{thisRegion}.xlsx')
df_unit = pd.read_excel(f'{ROOT_data}/UnitsTable_r{thisRID}_s{str(thisSID).zfill(2)}_{thisRegion}.xlsx')
df_react = pd.read_excel(f'{ROOT_data}/ReactTable_r{thisRID}_s{str(thisSID).zfill(2)}_{thisRegion}.xlsx')

#%% data processing
valid_rips = df_rip.RipID[df_rip.Filter>0]

df_react_valid = df_react[df_react['RipID'].isin(valid_rips) & df_react['Type']>0]
df_unit_valid = df_unit[(df_unit.Type>0)]


df_unit_valid['TT'] = df_unit_valid['TT'].astype(str).str.zfill(2)
df_unit_valid['Unit'] = df_unit_valid.Unit.astype(str).str.zfill(2)
df_unit_valid['TT-Unit'] = df_unit_valid['TT'] +  '-' + df_unit_valid['Unit']
df_unit_valid_dv = df_unit_valid[(df_unit_valid.PeakArea==3)]

df_react_valid['TT'] = df_react_valid['TT'].astype(str).str.zfill(2)
df_react_valid['UnitID'] = df_react_valid.UnitID.astype(str).str.zfill(2)
df_react_valid['TT-Unit'] = df_react_valid['TT'] +  '-' + df_react_valid['UnitID']
df_react_valid_dv = df_react_valid[df_react_valid['TT-Unit'].isin(df_unit_valid_dv['TT-Unit'])]

df_react_pivot_dv = pd.pivot_table(df_react_valid_dv, index='TT-Unit', columns = 'RipID',values='Region',aggfunc='count')
df_react_pivot_dv = ~df_react_pivot_dv.isna()*1

df_react_pivot_dv_rdi =df_react_pivot_dv.mul(df_unit_valid_dv.RDI_ZB.values,axis=0)
    

#%%
unit=list()
for i in range(df_react_pivot_dv.shape[0]):
    temp = df_react_pivot_dv.iloc[i,:]
    u = np.where(temp)
    unit.append(u)
#%%
pivot_overlap=pd.DataFrame(index = df_react_pivot_dv.index, columns=df_react_pivot_dv.index)
for i in range(len(unit)):
    for j in range(len(unit)):
        pivot_overlap.iat[i,j] = len(np.intersect1d(unit[i],unit[j]))
        pivot_overlap.iat[j,i] = len(np.intersect1d(unit[i],unit[j]))
    # pivot_overlap.iat[i,i] = 0

#%%
# data plotting
# sns.scatterplot(df_unit_valid_dv.RDI_ZB,df_unit_valid_dv.RDI_PM, hue = df_unit_valid.NumRipples_1,legend=0)
# sns.color_palette("tab10")
# plt.show()

pivot_overlap = pivot_overlap.astype('int64')
sns.heatmap(pivot_overlap, cmap ='Blues', linewidths = 0.30, annot = True)
plt.figure(figsize=(60,80))
plt.show()

sns.heatmap(df_react_pivot_dv_rdi, cmap ='viridis', annot = False)
