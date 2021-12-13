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
df_react_pivot_dv[df_react_pivot_dv==0]=np.nan

df_react_pivot_dv_rdi =df_react_pivot_dv.mul(df_unit_valid_dv.RDI_ZB.values,axis=0)

df_rip['meanRDI_ZB'] = 0
df_rip['meanRDI_PM'] = 0
df_rip['meanRDI_LR'] = 0
for i in range(df_rip.shape[0]):
    thisRipID = df_rip.RipID[i]
    thisUnitID = df_react_valid_dv['TT-Unit'][df_react_valid_dv['RipID']==thisRipID]
    thisUnits = df_unit_valid_dv[df_unit_valid_dv['TT-Unit'].isin(thisUnitID)]
    df_rip['meanRDI_ZB'].iloc[i] = np.mean(thisUnits.RDI_ZB)
    df_rip['meanRDI_PM'].iloc[i] = np.mean(thisUnits.RDI_PM)
    df_rip['meanRDI_LR'].iloc[i] = np.mean(thisUnits.RDI_LR)
    
df_rip_valid = df_rip[df_rip['RipID'].isin(df_react_pivot_dv_rdi.columns)]

#
unit=list()
for i in range(df_react_pivot_dv.shape[0]):
    temp = df_react_pivot_dv.iloc[i,:]
    u = np.where(temp)
    unit.append(u)
#
pivot_overlap=pd.DataFrame(index = df_react_pivot_dv.index, columns=df_react_pivot_dv.index)
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
sns.heatmap(pivot_overlap, cmap ='Blues', linewidths = 0.30, annot = True, mask=mask)



df_unit_valid_dv=df_unit_valid_dv.set_axis(df_unit_valid_dv['TT-Unit'],axis=0)
df_rip_valid=df_rip_valid.set_axis(df_rip_valid['RipID'],axis=0)
# RDI distribution plot



x = df_unit_valid_dv.loc[df_react_valid_dv['TT-Unit'],['RDI_ZB','RDI_PM','RDI_LR']]
x=x.reset_index(drop=True)
y = df_react_valid_dv['RipID']
y=y.reset_index(drop=True)
h = df_rip_valid.loc[df_react_valid_dv['RipID'],['Context','meanRDI_ZB','meanRDI_PM','meanRDI_LR']]
h=h.reset_index(drop=True)
dat = pd.concat([x,y,h],axis=1)
#%%
dat['rank']=dat[['meanRDI_ZB','RipID']].apply(tuple,axis=1).rank(method='dense')
dat=dat.sort_values('rank')

plt.figure(figsize=(20,10))

ax=sns.scatterplot(y='RDI_ZB',x='rank',hue='Context',data=dat,palette='tab10')
handles, labels  =  ax.get_legend_handles_labels()
ax.legend(handles, ['Zebra','Pebbles','Bamboo','Mountain'])
ax.set_xlabel("Ripples")
ax.set_xticks(dat['rank'].unique())
ax.set_xticklabels(dat['RipID'].unique(), fontsize=8,rotation=90)
ax.set_ylim(-1,1)
ax.grid(which="major",alpha=0.5)
ax.grid()

