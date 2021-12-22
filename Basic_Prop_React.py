# -*- coding: utf-8 -*-
"""
Created on Tue Dec 21 17:28:19 2021

0. Reactivated cell들이 갖는 특성 확인

@author: user
"""

#%% import libraries
import numpy as np 
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import itertools as ite

#%%

temp = df_unit_valid.loc[:,'NumRipples_1':'NumRipples_4'].sum(axis=1)
temp = pd.concat([temp,df_unit_valid.Session],axis=1)
temp = temp.rename(columns={0:'ReactRips'})
for thisSID in temp.Session.unique():
    temp.loc[temp['Session']==thisSID,'TotRips']=sum(df_rip_valid.Session==thisSID)
#%%
plt.figure(figsize=(10,8))
ax = sns.histplot(temp.ReactRips / temp.TotRips,stat='probability')
ax.set_xlabel('Participation in Ripple')

#%%
clist = ['r','b','y','g']
f,axes=plt.subplots(2,2,figsize=(10,8),sharex=True,sharey=True)
for Cxt in range(1,5):
    temp = df_unit_valid.loc[:,f'NumRipples_{Cxt}']
    temp = pd.concat([temp,df_unit_valid.Session],axis=1)
    
    for thisSID in temp.Session.unique():
        temp.loc[temp['Session']==thisSID,'TotRips']=sum((df_rip_valid.Session==thisSID) & (df_rip_valid.Context==Cxt))


    sns.histplot(temp[f'NumRipples_{Cxt}'] / temp.TotRips,stat='probability', color=clist[Cxt-1],ax=axes[divmod(Cxt-1,2)])
    
plt.xlabel('Participation in Ripple')
