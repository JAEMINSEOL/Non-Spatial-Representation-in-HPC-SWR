# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 17:00:11 2022
3. Activated cells의 co-activity rate 계산
@author: JM_Seol
"""

#%% import libraries
import numpy as np 
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import itertools as ite
import Mod_SWR as swr
#%% import data
ROOT_data = 'D:/HPC-SWR project/Processed Data'
thisRID=561
thisSID='all'
thisRegion='CA1'
   
df_rip_valid = pd.read_excel(f'{ROOT_data}/RipplesTable_r{thisRID}_{thisSID}_{thisRegion}_v.xlsx')
df_unit_valid = pd.read_excel(f'{ROOT_data}/UnitsTable_r{thisRID}_{thisSID}_{thisRegion}_v.xlsx')
df_act_valid = pd.read_excel(f'{ROOT_data}/ActTable_r{thisRID}_{thisSID}_{thisRegion}_v.xlsx')
dat = pd.read_excel(f'{ROOT_data}/ReactTable_r{thisRID}_all_{thisRegion}_2.xlsx')

#%%
df_unit_valid_stem = df_unit_valid[((df_unit_valid.PeakArea==3) | (df_unit_valid.PeakArea==2))]

df_act_valid_stem = df_act_valid[df_act_valid['TT-Unit'].isin(df_unit_valid_stem['TT-Unit'])]

df_act_pivot_stem = pd.pivot_table(df_act_valid_stem, index='TT-Unit', columns = 'RipID',values='Region',aggfunc='count')
df_act_pivot_stem = ~df_act_pivot_stem.isna()*1
df_act_pivot_stem[df_act_pivot_stem==0]=np.nan

df_act_pivot_stem_rdi =df_act_pivot_stem.mul(df_unit_valid_stem.RDI_ZB.values,axis=0)

#%%

df_act_comb = pd.DataFrame(columns=[0,1])
for thisRipID in df_act_valid['RipID'].unique():
    temp = df_act_valid[df_act_valid['RipID']==thisRipID] 
    for subset in ite.combinations(temp['TT-Unit'],2):       
        df_act_comb=df_act_comb.append(pd.DataFrame(list(subset)).T)

typ = ['RDI_ZB','RDI_PM']
for i in range(len(df_act_comb)):
    u1 = df_act_comb.iloc[i,0]
    u2 = df_act_comb.iloc[i,1]
    x = [df_unit_valid.loc[df_unit_valid['TT-Unit']==u1,typ[0]].iloc[0], 
            df_unit_valid.loc[df_unit_valid['TT-Unit']==u2,typ[0]].iloc[0]]
    y= [df_unit_valid.loc[df_unit_valid['TT-Unit']==u1,typ[1]].iloc[0],
        df_unit_valid.loc[df_unit_valid['TT-Unit']==u2,typ[1]].iloc[0]]
    plt.plot(x, y, marker = 'o',color='k',alpha=0.03)
plt.xlabel("RDI_ZB")
plt.ylabel("RDI_PM")
plt.show()


#
unit=list()
for i in range(df_act_pivot_stem.shape[0]):
    temp = df_act_pivot_stem.iloc[i,:]
    u = np.where(~np.isnan(temp))
    unit.append(u)
#
pivot_overlap=pd.DataFrame(index = df_act_pivot_stem.index, columns=df_act_pivot_stem.index)
for i in range(len(unit)):
    for j in range(len(unit)):
        pivot_overlap.iat[i,j] = len(np.intersect1d(unit[i],unit[j]))
        pivot_overlap.iat[j,i] = len(np.intersect1d(unit[i],unit[j]))
    # pivot_overlap.iat[i,i] = 0
    
#%% data plotting
# sns.scatterplot(df_unit_valid_dv.RDI_ZB,df_unit_valid_dv.RDI_PM, hue = df_unit_valid.NumRipples_1,legend=0)
# sns.color_palette("tab10")
# plt.show()

# unit co-reactivate heatmap
pivot_overlap = pivot_overlap.astype('int64')
mask = np.triu(np.ones_like(pivot_overlap, dtype=np.bool))
sns.heatmap(pivot_overlap, cmap ='Blues', linewidths = 0.30, annot = False, mask=mask)
        
