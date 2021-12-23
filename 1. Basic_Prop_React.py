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

df_unit = df_unit[df_unit.Type!=0]
#%% RDI, field peak position distribution
f,axes=plt.subplots(2,2,figsize=(10,8),sharey=True)
for i in range(1,5):
   sns.histplot(df_unit.iloc[:,i+8], stat='probability',color='k',ax=axes[divmod(i-1,2)],binwidth=3 if i==1 else 0.1)
   sns.histplot(df_unit_valid.iloc[:,i+8],stat='probability', color='r',ax=axes[divmod(i-1,2)],binwidth=3 if i==1 else 0.1)
   if i>1:
       plt.xlim(-1,1)
#%%
temp=list(np.zeros(5))
for i in range(1,6):
    temp[i-1] = sum(df_rip_valid.Session==i)
temp=pd.DataFrame(data=temp)
plt.figure(figsize=(10,8))
ax = sns.barplot(data=temp.T)
ax.set_xticklabels(['Session1','Session2','Session3','Session4','Session5'])
plt.ylabel('Ripples')
#%% Ripple participation rate (all scenes)
temp = df_unit_valid.loc[:,'NumRipples_1':'NumRipples_4'].sum(axis=1)
temp = pd.concat([temp,df_unit_valid.Session],axis=1)
temp = temp.rename(columns={0:'ReactRips'})
for thisSID in temp.Session.unique():
    temp.loc[temp['Session']==thisSID,'TotRips']=sum(df_rip_valid.Session==thisSID)

plt.figure(figsize=(10,8))
ax = sns.histplot(temp.ReactRips / temp.TotRips,stat='probability')
ax.set_xlabel('Participation in Ripple')

#%% Ripple participation rate (each scene)
clist = ['#5AB7D4','#F79334','#00506A','#9A4700']
f,axes=plt.subplots(2,2,figsize=(10,8),sharex=True,sharey=True)
for Cxt in range(1,5):
    temp = df_unit_valid.loc[:,['TT-Unit', f'NumRipples_{Cxt}']]
    temp = pd.concat([temp,df_unit_valid.Session],axis=1)
    
    for thisSID in temp.Session.unique():
        temp.loc[temp['Session']==thisSID,'TotRips']=sum((df_rip_valid.Session==thisSID) & (df_rip_valid.Context==Cxt))


    sns.histplot(temp[f'NumRipples_{Cxt}'] / temp.TotRips,stat='probability', color=clist[Cxt-1],ax=axes[divmod(Cxt-1,2)])
    
plt.xlabel('Participation in Ripple')


    
