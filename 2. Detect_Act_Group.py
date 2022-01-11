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
import Mod_SWR as swr
#%% import data
ROOT_data = 'D:/HPC-SWR project/Processed Data'
thisRID=561
thisSID='all'
thisRegion='CA1'
   
df_rip_valid = pd.read_excel(f'{ROOT_data}/RipplesTable_r{thisRID}_{thisSID}_{thisRegion}_v.xlsx')
df_unit_valid = pd.read_excel(f'{ROOT_data}/UnitsTable_r{thisRID}_{thisSID}_{thisRegion}_v.xlsx')
df_act_valid = pd.read_excel(f'{ROOT_data}/ActTable_r{thisRID}_{thisSID}_{thisRegion}_v.xlsx')


#%% make react sheet & save

# RDI distribution plot
df_unit_valid = df_unit_valid.rename(index=df_unit_valid['TT-Unit'])
x = df_unit_valid.loc[df_act_valid['TT-Unit'],
                      ['TT-Unit','PeakArea','RDI_ZB','RDI_PM','RDI_LR']]
x=x.reset_index(drop=True)
y = df_act_valid.loc[:,['RipID','SpkTime']]
y=y.reset_index(drop=False)

df_rip_valid = df_rip_valid.rename(index=df_rip_valid['RipID'])
h = df_rip_valid.loc[df_act_valid['RipID'],
                     ['Context','meanRDI_ZB','meanRDI_PM','meanRDI_LR']]
h=h.reset_index(drop=False)
dat = pd.concat([x,y,h],axis=1)
dat['RDI_LR']=dat['RDI_LR'].fillna(0)
dat['meanRDI_ZB']=dat['meanRDI_ZB'].fillna(0)
dat['meanRDI_PM']=dat['meanRDI_PM'].fillna(0)
dat['meanRDI_LR']=dat['meanRDI_LR'].fillna(0)

dat['rank']=dat[['RipID']].apply(tuple,axis=1).rank(method='dense')
dat=dat.sort_values('rank')

dat.to_excel(f'{ROOT_data}/ReactTable_r{thisRID}_all_{thisRegion}_2.xlsx')

#%% scene RDI distribution for each scene
TargRDI=['ZB','PM']
clist = ['#5AB7D4','#F79334','#00506A','#9A4700']
CxtList = ['Zebra', 'Pebbles', 'Bamboo', 'Mountains']
datasets={}
f,axes=plt.subplots(2,2,figsize=(10,8),sharex=True,sharey=True)
for Cxt in range(1,5):
    
    dat_part = dat[(dat['Context']==Cxt) & ((dat['PeakArea']==2) | (dat['PeakArea']==3))]
    dat_part = dat_part.drop_duplicates(['TT-Unit'])
    
    p=sp.stats.wilcoxon(dat_part[f'RDI_{TargRDI[np.mod(Cxt+1,2)]}'])
    sns.histplot(dat_part[f'RDI_{TargRDI[np.mod(Cxt+1,2)]}'],binrange=(-1.2,1.2),bins=12,
                 color=clist[Cxt-1],ax=axes[divmod(Cxt-1,2)], stat='probability')
    axes[divmod(Cxt-1,2)].plot((0,0),(0,0.5),color='k')
    axes[divmod(Cxt-1,2)].set_title(f'{CxtList[Cxt-1]}, n={len(dat_part)}, p={round(p[1],3)}')
    axes[divmod(Cxt-1,2)].set_ylim(0,0.4)
    
    datasets[Cxt] = dat_part[f'RDI_{TargRDI[np.mod(Cxt+1,2)]}']
plt.suptitle('RDI distribution for activated units in each scene ripple')

opts = {'title_hist': 'Unit distribution histogram','title_dist': 'CDF plot (K-S test '}
d1= datasets[1]
d2 = datasets[3]
opts_indiv = {**opts, **{'bins':20,'range':(-1.2,1.2),'fontsize':15,'xlab': 'RDI_ZB','cmap' : [clist[0],clist[2]],
                         'legend':[f'{CxtList[0]} (n={len(d1)})',f'{CxtList[2]} (n={len(d2)})']}}
axes = swr.DrawDist_2samp(True, d1,d2,**opts_indiv)

d1= datasets[2]
d2 = datasets[4]
opts_indiv = {**opts, **{'bins':20,'range':(-1.2,1.2),'fontsize':15,'xlab': 'RDI_PM','cmap' : [clist[1],clist[3]],
                         'legend':[f'{CxtList[1]} (n={len(d1)})',f'{CxtList[3]} (n={len(d2)})']}}
axes = swr.DrawDist_2samp(True, d1,d2,**opts_indiv)
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

