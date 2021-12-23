# -*- coding: utf-8 -*-
"""
Created on Thu Dec  9 10:43:10 2021

2. Activated group의 RDI 분포 확인
input: ReactTable_valid, RipplesTable_valid, UnitsTable_valid
output: Active group detecting plots
@author: JM_Seol
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
thisRID=561
thisSID='all'
thisRegion='CA1'
   
df_rip_valid = pd.read_excel(f'{ROOT_data}/RipplesTable_r{thisRID}_{thisSID}_{thisRegion}_v.xlsx')
df_unit_valid = pd.read_excel(f'{ROOT_data}/UnitsTable_r{thisRID}_{thisSID}_{thisRegion}_v.xlsx')
df_act_valid = pd.read_excel(f'{ROOT_data}/ActTable_r{thisRID}_{thisSID}_{thisRegion}_v.xlsx')

#%% data processing
df_unit_valid_dv = df_unit_valid[(df_unit_valid.PeakArea==3)]

df_act_valid_dv = df_act_valid[df_act_valid['TT-Unit'].isin(df_unit_valid_dv['TT-Unit'])]

df_act_pivot_dv = pd.pivot_table(df_act_valid_dv, index='TT-Unit', columns = 'RipID',values='Region',aggfunc='count')
df_act_pivot_dv = ~df_act_pivot_dv.isna()*1
df_act_pivot_dv[df_act_pivot_dv==0]=np.nan

df_act_pivot_dv_rdi =df_act_pivot_dv.mul(df_unit_valid_dv.RDI_ZB.values,axis=0)


#
unit=list()
for i in range(df_act_pivot_dv.shape[0]):
    temp = df_act_pivot_dv.iloc[i,:]
    u = np.where(~np.isnan(temp))
    unit.append(u)
#
pivot_overlap=pd.DataFrame(index = df_act_pivot_dv.index, columns=df_act_pivot_dv.index)
for i in range(len(unit)):
    for j in range(len(unit)):
        pivot_overlap.iat[i,j] = len(np.intersect1d(unit[i],unit[j]))
        pivot_overlap.iat[j,i] = len(np.intersect1d(unit[i],unit[j]))
    # pivot_overlap.iat[i,i] = 0


#%% make react sheet & save

# RDI distribution plot
x = df_unit_valid.loc[df_act_valid['TT-Unit'],['TT-Unit','PeakArea','RDI_ZB','RDI_PM','RDI_LR']]
x=x.reset_index(drop=True)
y = df_act_valid.loc[:,['RipID','SpkTime']]
y=y.reset_index(drop=True)
h = df_rip_valid.loc[df_act_valid['RipID'],['Context','meanRDI_ZB','meanRDI_PM','meanRDI_LR']]
h=h.reset_index(drop=True)
dat = pd.concat([x,y,h],axis=1)
dat['RDI_LR']=dat['RDI_LR'].fillna(0)
dat['meanRDI_ZB']=dat['meanRDI_ZB'].fillna(0)
dat['meanRDI_PM']=dat['meanRDI_PM'].fillna(0)
dat['meanRDI_LR']=dat['meanRDI_LR'].fillna(0)

dat['rank']=dat[['RipID']].apply(tuple,axis=1).rank(method='dense')
dat=dat.sort_values('rank')

dat.to_excel(f'{ROOT_data}/ReactTable_r{thisRID}_all_{thisRegion}_2.xlsx')


#%% Wilcoxon signed rank sum test & plotting
Targ='PM'
y=[[0,0],[0,0],[0,0],[0,0]]
for c in range(1,5):
    dat_part = dat[(dat['Context']==3) & (dat['PeakArea']==3)]
    dat_part['rank']=dat_part[[f'meanRDI_{Targ}','RipID']].apply(tuple,axis=1).rank(method='dense')
    
    dat_part['p_Wcx']=1
    dat_part['p_Wcx2']=1
    for thisRID in dat_part.RipID.unique():
        try:
            temp = dat_part[dat_part['RipID']==thisRID]
            # temp = temp.drop_duplicates(subset='TT-Unit')
            p=sp.stats.wilcoxon(temp[f'RDI_{Targ}'],alternative='greater')
            dat_part['p_Wcx'][dat_part['RipID']==thisRID]=p[1]
            p=sp.stats.wilcoxon(temp[f'RDI_{Targ}'],alternative='less')
            dat_part['p_Wcx2'][dat_part['RipID']==thisRID]=p[1]
        except:
            print(1)
    
    temp = dat_part
    temp['p_Wcx'] = temp['p_Wcx']<0.05
    temp['p_Wcx2'] = temp['p_Wcx2']<0.05
    t=temp[temp['p_Wcx']<2]
    y[c-1][0] = len(t['RipID'].unique())
    t=temp[temp['p_Wcx2']==1]
    y[c-1][1] = len(t['RipID'].unique())
y=pd.DataFrame(y)
#%%
plt.figure(figsize=(40,10))
ax=sns.scatterplot(y=f'RDI_{Targ}',x='rank',data=dat_part[dat_part['PeakArea']==3],color='k')
ax=sns.scatterplot(y=f'RDI_{Targ}',x='rank',data=dat_part[dat_part['PeakArea']==2],color='k',marker='x')
ax=sns.scatterplot(y='p_Wcx2',x='rank',data=temp[temp.p_Wcx2==True],color='r',marker='*',s=160)
ax=sns.scatterplot(y='p_Wcx',x='rank',data=temp[temp.p_Wcx==True],color='b',marker='*',s=160)
handles, labels  =  ax.get_legend_handles_labels()
ax.legend(handles, ['Stbox','Stem','Dv','Arm'])
ax.set_xlabel("Ripples")
ax.set_xticks(dat_part['rank'].unique())
ax.set_xticklabels(dat_part['RipID'].unique(), fontsize=7,rotation=90)
ax.set_xlim(-1,1)
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


df_act_comb = pd.DataFrame(columns=[0,1])
for thisRipID in df_act_valid_dv['RipID'].unique():
    temp = df_act_valid_dv[df_act_valid_dv['RipID']==thisRipID] 
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

