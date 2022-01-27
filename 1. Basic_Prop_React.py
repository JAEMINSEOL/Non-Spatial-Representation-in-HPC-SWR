# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 17:28:19 2021

1. Reactivated cell들이 갖는 특성 확인
input: ReactTable_valid, RipplesTable_valid, UnitsTable_valid
output: unit property plots

@author: JM_Seol
"""

#%% import libraries
import numpy as np 
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import Mod_SWR as swr
#%% import data
ROOT_data = 'D:/HPC-SWR project/Processed Data'
thisRID=561
thisSID='all'
thisRegion='CA1'

df_rip = pd.read_excel(f'{ROOT_data}/RipplesTable_r{thisRID}_{thisSID}_{thisRegion}.xlsx')
df_unit = pd.read_excel(f'{ROOT_data}/UnitsTable_r{thisRID}_{thisSID}_{thisRegion}.xlsx')
df_act = pd.read_excel(f'{ROOT_data}/ActTable_r{thisRID}_{thisSID}_{thisRegion}.xlsx')

df_rip_valid = pd.read_excel(f'{ROOT_data}/RipplesTable_r{thisRID}_{thisSID}_{thisRegion}_v.xlsx')
df_unit_valid = pd.read_excel(f'{ROOT_data}/UnitsTable_r{thisRID}_{thisSID}_{thisRegion}_v.xlsx')
df_act_valid = pd.read_excel(f'{ROOT_data}/ActTable_r{thisRID}_{thisSID}_{thisRegion}_v.xlsx')

df_clust_summ = pd.read_excel('D:/HPC-SWR project/Information Sheet/ClusterSummary_SWR.xlsx')

df_unit = df_unit[(df_unit.Type!=0) & ((df_unit.PeakArea==2) | (df_unit.PeakArea==3))]
df_unit_valid = df_unit_valid[(df_unit_valid.Type!=0) & ((df_unit_valid.PeakArea==2) | (df_unit_valid.PeakArea==3))]
#%% RDI, field peak position distribution
opts = {'legend':[f'All (n={len(df_unit)})',f'Activated (n={len(df_unit_valid)})'],
        'title_hist': 'Unit distribution histogram','title_dist': 'CDF plot (K-S test ',
        'cmap' : ['k','r']}

# SI score (1D out)
prop = 'SpaInfoScore1D out'
d1 = swr.Exp_from_ClustSum(df_clust_summ,thisRID,df_unit, prop)
d2 = swr.Exp_from_ClustSum(df_clust_summ,thisRID,df_unit_valid, prop)
opts_indiv = {**opts, **{'bins':10,'range':(0,2.5),'fontsize':15}}
axes = swr.DrawDist_2samp(False, d1,d2,**opts_indiv)

# Mean firing rate
prop = 'onmazeAvgFR1D out'
d1 = swr.Exp_from_ClustSum(df_clust_summ,thisRID,df_unit, prop)
d2 = swr.Exp_from_ClustSum(df_clust_summ,thisRID,df_unit_valid, prop)
opts_indiv = {**opts, **{'bins':24,'range':(0,12),'fontsize':15,'xlab': 'Mean FR (1D out)'}}
axes = swr.DrawDist_2samp(False, d1,d2,**opts_indiv)

# Peak firing rate
prop = 'onmazeMaxFR1D out'
d1 = swr.Exp_from_ClustSum(df_clust_summ,thisRID,df_unit, prop)
d2 = swr.Exp_from_ClustSum(df_clust_summ,thisRID,df_unit_valid, prop)
opts_indiv = {**opts, **{'bins':35,'range':(0,35),'fontsize':15,'xlab': 'Peak FR (1D out)'}}
axes = swr.DrawDist_2samp(False, d1,d2,**opts_indiv)

# Number of Spike
prop = 'numOfSpk1D out'
d1 = swr.Exp_from_ClustSum(df_clust_summ,thisRID,df_unit, prop)
d2 = swr.Exp_from_ClustSum(df_clust_summ,thisRID,df_unit_valid, prop)
opts_indiv = {**opts, **{'bins':22,'range':(0,2200),'fontsize':15,'xlab': 'Num. of spikes'}}
axes = swr.DrawDist_2samp(False, d1,d2,**opts_indiv)

# Sparsity
prop = 'Sparsity'
d1 = swr.Exp_from_ClustSum(df_clust_summ,thisRID,df_unit, prop)
d2 = swr.Exp_from_ClustSum(df_clust_summ,thisRID,df_unit_valid, prop)
opts_indiv = {**opts, **{'bins':20,'range':(0,1),'fontsize':15,'xlab': 'Sparsity'}}
axes = swr.DrawDist_2samp(False, d1,d2,**opts_indiv)

# Peak to Valley
prop = ' Peak (peak to valley)'
d1 = swr.Exp_from_ClustSum(df_clust_summ,thisRID,df_unit, prop)
d2 = swr.Exp_from_ClustSum(df_clust_summ,thisRID,df_unit_valid, prop)
opts_indiv = {**opts, **{'bins':20,'range':(0,400),'fontsize':15,'xlab': 'Peak to valley'}}
axes = swr.DrawDist_2samp(False, d1,d2,**opts_indiv)

# Peak Position
opts_indiv = {**opts, **{'bins':21,'range':(0,42),'fontsize':15,'xlab': 'Peak Position (bin)'}}
axes = swr.DrawDist_2samp(False, df_unit.iloc[:,9],df_unit_valid.iloc[:,9],**opts_indiv)

# RDIs
for i in range(10,13):
    opts_indiv = {**opts, **{'bins':24,'range':(-1.2, 1.2),'fontsize':15}}
    axes = swr.DrawDist_2samp(True, df_unit.iloc[:,i],df_unit_valid.iloc[:,i],**opts_indiv)
#%% Number of Ripples for each scene
temp=list(np.zeros(5))
for i in range(1,6):
    temp[i-1] = sum(df_rip_valid.Session==i)
temp=pd.DataFrame(data=temp)
plt.figure(figsize=(10,8))
ax = sns.barplot(data=temp.T)
ax.set_xticklabels(['Session1','Session2','Session3','Session4','Session5'])
plt.ylabel('Ripples')
plt.title('Number of Ripples in Rat561 Sessions')
#%% Ripple participation rate (all scenes)
temp = df_unit_valid.loc[:,'NumRipples_1':'NumRipples_4'].sum(axis=1)
temp = pd.concat([temp,df_unit_valid.Session],axis=1)
temp = temp.rename(columns={0:'ReactRips'})
for thisSID in temp.Session.unique():
    temp.loc[temp['Session']==thisSID,'TotRips']=sum(df_rip_valid.Session==thisSID)

plt.figure(figsize=(10,8))
ax = sns.histplot(temp.ReactRips / temp.TotRips,stat='probability',bins=10,binrange=(0, 0.5))
ax.set_xlabel('Ripple Participation Rate')
plt.rc('font',size=15)
plt.title('Participation in Ripple')
#%% Ripple participation rate (each scene)
clist = ['#5AB7D4','#F79334','#00506A','#9A4700']
CxtList = ['Zebra', 'Pebbles', 'Bamboo', 'Mountains']

f,axes=plt.subplots(2,2,figsize=(10,8),sharex=True,sharey=True)
for Cxt in range(1,5):
    temp = df_unit_valid.loc[:,['TT-Unit', f'NumRipples_{Cxt}']]
    temp = pd.concat([temp,df_unit_valid.Session],axis=1)
    
    for thisSID in temp.Session.unique():
        temp.loc[temp['Session']==thisSID,'TotRips']=sum((df_rip_valid.Session==thisSID) & (df_rip_valid.Context==Cxt))


    sns.histplot(temp[f'NumRipples_{Cxt}'] / temp.TotRips,stat='probability',
                 bins=12,binrange=(0, 0.6), color=clist[Cxt-1],ax=axes[divmod(Cxt-1,2)])
    axes[divmod(Cxt-1,2)].set_title(CxtList[Cxt-1])
    
    
plt.xlabel('')
plt.suptitle('Participation in Ripple for Each Scene')
f.text(0.5, 0.04, 'Ripple Participation Rate', ha='center')



#%% Ripple participation rate (each scene) - CDF plot + K-S test
clist = ['#5AB7D4','#F79334','#00506A','#9A4700']
CxtList = ['Zebra', 'Pebbles', 'Bamboo', 'Mountains']

for Cxt1 in range(1,5):
    for Cxt2 in range(1,5):
        if Cxt1>=Cxt2:
            continue
        
        temp = df_unit_valid.loc[:,['TT-Unit', f'NumRipples_{Cxt1}']]
        temp1 = pd.concat([temp,df_unit_valid.Session],axis=1)
        for thisSID in temp1.Session.unique():
            temp1.loc[temp1['Session']==thisSID,'TotRips']=sum((df_rip_valid.Session==thisSID) & (df_rip_valid.Context==Cxt1))
        
        temp = df_unit_valid.loc[:,['TT-Unit', f'NumRipples_{Cxt2}']]
        temp2 = pd.concat([temp,df_unit_valid.Session],axis=1)
        for thisSID in temp2.Session.unique():
            temp2.loc[temp2['Session']==thisSID,'TotRips']=sum((df_rip_valid.Session==thisSID) & (df_rip_valid.Context==Cxt2))

        d1 = temp1[f'NumRipples_{Cxt1}'] / temp1.TotRips
        d2 = temp2[f'NumRipples_{Cxt2}'] / temp2.TotRips
        
        opts = {'legend':[f'{CxtList[Cxt1-1]}',f'{CxtList[Cxt2-1]}'],
                'title_hist': 'Unit distribution histogram','title_dist': 'CDF plot (K-S test ',
                'cmap' : [clist[Cxt1-1],clist[Cxt2-1]]}
        opts_indiv = {**opts, **{'bins':12,'range':(0,0.6),'fontsize':15,'xlab': 'Ripple Participation Rate'}}
        axes = swr.DrawDist_2samp(False, d1,d2,**opts_indiv)


#%% ripple occurrence numbers for each trials 



for s in df_rip_valid.Session.unique():
    
    
    df_TrialTime = pd.read_excel(f'D:/HPC-SWR project/Information Sheet/TrialTime_561-0{s}.xlsx')
    df_TrialTime.columns=['s1', 's2', 's3', 's4', 'fw', 'stbox']
    ITI_Time = df_TrialTime['s1'].iloc[1:].reset_index(drop=True) -  df_TrialTime['stbox'].iloc[0:-1].reset_index(drop=True)
    df_rip_temp = df_rip_valid[df_rip_valid.Session==s]
    m= len(df_TrialTime)
    temp = df_rip_temp.Trial.unique()
    rips_trial = pd.Series(0, index=range(m))
    
    for t in range(m):      
        rips_trial.iloc[t] = sum(df_rip_temp.Trial==t+1)
    rips_trial.index=range(1,m+1)
    fig, ax = plt.subplots(1,3,figsize=(14,10))
    ax[0].barh(range(1,m+1),rips_trial, color="k")
    ax[0].invert_yaxis()
    ax[0].set_yticks(range(1,m+1),labels = range(1,m+1),size=8)
    ax[0].set_ylabel('trial')
    ax[0].set_xlabel('ripples')
    ax[0].set_title(f'# of ripples_total, 561-0{s}')
    
       
    ax[1].barh(range(1,m),ITI_Time, color="k")
    ax[1].invert_yaxis()
    ax[1].set_yticks(range(1,m),labels = range(1,m),size=8)
    ax[1].set_ylabel('trial')
    ax[1].set_xlabel('ITI duration')
    ax[1].set_title(f'ITI duration(s), 561-0{s}')
    
    ax[2].barh(range(1,m),rips_trial.iloc[0:-1].reset_index(drop=True)/ITI_Time, color="k")
    ax[2].invert_yaxis()
    ax[2].set_yticks(range(1,m),labels = range(1,m),size=8)
    ax[2].set_ylabel('trial')
    ax[2].set_xlabel('ripple rate')
    ax[2].set_title(f'ripple rate(ripple/s), 561-0{s}')
    
    fig, ax = plt.subplots(1,4,figsize=(18,10), sharex=True)
    plt.suptitle(f'# of ripples_total, 561-0{s}')
    for c in range(0,4):
        if c==2:
            c2=1
        else:
            if c==1:
                c2=2
            else:
                c2=c
        df_rip_temp = df_rip_valid[(df_rip_valid.Session==s) & (df_rip_valid.Context==c2+1)]
        temp = df_rip_temp.Trial.unique()
        rips_trial = []
        for t in temp:       
            rips_trial.append(sum(df_rip_temp.Trial==t))
        ax[c].barh(range(len(temp)),rips_trial, color=clist[c2])
        ax[c].invert_yaxis()
        ax[c].set_yticks(range(len(temp)),labels = temp)
        ax[c].set_ylabel('trial')
        ax[c].set_xlabel('ripples')
        ax[c].set_title(f'{CxtList[c2]}')
