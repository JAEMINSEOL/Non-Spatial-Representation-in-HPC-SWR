# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 10:43:10 2021

1. Reactivated group의 RDI 분포 확인

21.12.20 - Wilcoxon signed ranksum test 결과 추가
@author: user
"""

#%% import libraries
import numpy as np 
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import itertools as ite
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
    thisUnitID = df_react_valid['TT-Unit'][df_react_valid['RipID']==thisRipID]
    thisUnits = df_unit_valid[df_unit_valid['TT-Unit'].isin(thisUnitID)]
    df_rip['meanRDI_ZB'].iloc[i] = np.mean(thisUnits.RDI_ZB)
    df_rip['meanRDI_PM'].iloc[i] = np.mean(thisUnits.RDI_PM)
    df_rip['meanRDI_LR'].iloc[i] = np.mean(thisUnits.RDI_LR)
    
df_rip_valid = df_rip[df_rip['Filter']>0]

#
unit=list()
for i in range(df_react_pivot_dv.shape[0]):
    temp = df_react_pivot_dv.iloc[i,:]
    u = np.where(~np.isnan(temp))
    unit.append(u)
#
pivot_overlap=pd.DataFrame(index = df_react_pivot_dv.index, columns=df_react_pivot_dv.index)
for i in range(len(unit)):
    for j in range(len(unit)):
        pivot_overlap.iat[i,j] = len(np.intersect1d(unit[i],unit[j]))
        pivot_overlap.iat[j,i] = len(np.intersect1d(unit[i],unit[j]))
    # pivot_overlap.iat[i,i] = 0


#%% make react sheet & save
df_unit_valid_dv=df_unit_valid_dv.set_axis(df_unit_valid_dv['TT-Unit'],axis=0)
df_unit_valid=df_unit_valid.set_axis(df_unit_valid['TT-Unit'],axis=0)
df_rip_valid=df_rip_valid.set_axis(df_rip_valid['RipID'],axis=0)
df_react_valid=df_react_valid[df_react_valid['RipID'].isin(df_rip_valid['RipID'])]
# RDI distribution plot



x = df_unit_valid.loc[df_react_valid['TT-Unit'],['TT-Unit','PeakArea','RDI_ZB','RDI_PM','RDI_LR']]
x=x.reset_index(drop=True)
y = df_react_valid['RipID']
y=y.reset_index(drop=True)
h = df_rip_valid.loc[df_react_valid['RipID'],['Context','meanRDI_ZB','meanRDI_PM','meanRDI_LR']]
h=h.reset_index(drop=True)
dat = pd.concat([x,y,h],axis=1)
dat['RDI_LR']=dat['RDI_LR'].fillna(0)
dat['meanRDI_ZB']=dat['meanRDI_ZB'].fillna(0)
dat['meanRDI_PM']=dat['meanRDI_PM'].fillna(0)
dat['meanRDI_LR']=dat['meanRDI_LR'].fillna(0)

dat['rank']=dat[['RipID']].apply(tuple,axis=1).rank(method='dense')
dat=dat.sort_values('rank')

# dat.to_excel(f'{ROOT_data}/ReactTable_r{thisRID}_s{str(thisSID).zfill(2)}_{thisRegion}_2.xlsx')

#%% Wilcoxon signed rank sum test
Targ='LR'
dat['p_Wcx']=1
for thisRID in dat.RipID.unique():
    temp = dat[dat['RipID']==thisRID]
    # temp = temp.drop_duplicates(subset='TT-Unit')
    p=sp.stats.wilcoxon(temp[f'RDI_{Targ}'])
    dat['p_Wcx'][dat['RipID']==thisRID]=p[1]

#%% plotting

plt.figure(figsize=(40,10))
dat_part = dat[(dat['Context']!=0)]
dat_part['rank']=dat_part[[f'meanRDI_{Targ}','RipID']].apply(tuple,axis=1).rank(method='dense')
temp = dat_part
temp['p_Wcx'] = temp['p_Wcx']<0.05

ax=sns.scatterplot(y=f'RDI_{Targ}',x='rank',data=dat_part[dat_part['PeakArea']==3],color='k')
ax=sns.scatterplot(y=f'RDI_{Targ}',x='rank',data=dat_part[dat_part['PeakArea']==2],color='k',marker='x')
ax=sns.scatterplot(y='p_Wcx',x='rank',data=temp[temp.p_Wcx==True],color='r',marker='*',s=160)
handles, labels  =  ax.get_legend_handles_labels()
ax.legend(handles, ['Stbox','Stem','Dv','Arm'])
ax.set_xlabel("Ripples")
ax.set_xticks(dat_part['rank'].unique())
ax.set_xticklabels(dat_part['RipID'].unique(), fontsize=7,rotation=90)
ax.set_ylim(-1,1.1)
ax.set_xlim(0.5,dat_part['rank'].max()+0.5)
ax.grid(which="major",alpha=0.5)
ax.grid()

ax.axhline(0, color='k', linewidth=1,alpha=1)
cxtshade = dat_part[dat_part['Context']==1]['rank'].unique().astype(int)
for i in cxtshade:
    ax.axvline(i, color='r', linewidth=5,alpha=0.3)
    
cxtshade = dat_part[dat_part['Context']==2]['rank'].unique().astype(int)
for i in cxtshade:
    ax.axvline(i, color='b', linewidth=5,alpha=0.3)
    
cxtshade = dat_part[dat_part['Context']==3]['rank'].unique().astype(int)
for i in cxtshade:
    ax.axvline(i, color='y', linewidth=5,alpha=0.3)
    
cxtshade = dat_part[dat_part['Context']==4]['rank'].unique().astype(int)
for i in cxtshade:
    ax.axvline(i, color='g', linewidth=5,alpha=0.3)
    
#%% 


df_react_comb = pd.DataFrame(columns=[0,1])
for thisRipID in df_react_valid_dv['RipID'].unique():
    temp = df_react_valid_dv[df_react_valid_dv['RipID']==thisRipID] 
    for subset in ite.combinations(temp['TT-Unit'],2):       
        df_react_comb=df_react_comb.append(pd.DataFrame(list(subset)).T)

typ = ['RDI_ZB','RDI_PM']
for i in range(len(df_react_comb)):
    u1 = df_react_comb.iloc[i,0]
    u2 = df_react_comb.iloc[i,1]
    x = [df_unit_valid.loc[df_unit_valid['TT-Unit']==u1,typ[0]].iloc[0], 
            df_unit_valid.loc[df_unit_valid['TT-Unit']==u2,typ[0]].iloc[0]]
    y= [df_unit_valid.loc[df_unit_valid['TT-Unit']==u1,typ[1]].iloc[0],
        df_unit_valid.loc[df_unit_valid['TT-Unit']==u2,typ[1]].iloc[0]]
    plt.plot(x, y, marker = 'o',color='k',alpha=0.03)
plt.xlabel("RDI_ZB")
plt.ylabel("RDI_PM")
plt.show()        
# ax=sns.scatterplot(y='RDI_LR',x='rank')
#%% data plotting
# sns.scatterplot(df_unit_valid_dv.RDI_ZB,df_unit_valid_dv.RDI_PM, hue = df_unit_valid.NumRipples_1,legend=0)
# sns.color_palette("tab10")
# plt.show()

# unit co-reactivate heatmap
pivot_overlap = pivot_overlap.astype('int64')
mask = np.triu(np.ones_like(pivot_overlap, dtype=np.bool))
sns.heatmap(pivot_overlap, cmap ='Blues', linewidths = 0.30, annot = True, mask=mask)


#%%






#%%

