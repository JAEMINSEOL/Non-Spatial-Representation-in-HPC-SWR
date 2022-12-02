# -*- coding: utf-8 -*-
"""
Created on Thu Nov 17 16:44:09 2022

@author: Jaemin
"""
#%% import libraries
import numpy as np 
import scipy as sp
from scipy import stats
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import itertools as ite
import Mod_SWR as swr
from matplotlib.colors import ListedColormap
from matplotlib import cm
import re
from PIL import ImageColor as imgC

#%% import data
ROOT_data = 'D:/HPC-SWR project/Processed Data'
ROOT_info = 'D:/HPC-SWR project/Information Sheet'
thisRID=561
thisSID='all'
thisRegion='CA1'
   
df_rip_valid = pd.read_excel(f'{ROOT_data}/RipplesTable_{thisRegion}_ForAnalysis.xlsx')
df_unit_valid = pd.read_excel(f'{ROOT_data}/UnitsTable_{thisRegion}_ForAnalysis.xlsx')
df_act_valid = pd.read_excel(f'{ROOT_data}/ReactTable_{thisRegion}.xlsx')
df_unit_comb=pd.read_excel(f'{ROOT_data}/ReactPair.xlsx')
df_session_list = pd.read_excel(f'{ROOT_info}/SessionList_SWR.xlsx')

Color_List = {"Zebra" : '#CF404D', 'Bamboo' : '#D06000', 'Pebbles' : '#00846D', 'Mountains' : '#404F78',
              'Left' : "#872300", 'Right' : '#005500', 'Forest' : '#009834', 'City' : '#008644'}
#%%
df_session_list['mean_RDI_LScene']=np.nan
df_session_list['std_RDI_LScene']=np.nan
df_session_list['mean_RDI_RScene']=np.nan
df_session_list['std_RDI_RScene']=np.nan
df_session_list['mean_RDI_LR']=np.nan
df_session_list['std_RDI_LR']=np.nan

for index, thisSession in df_session_list.iterrows():
    df=df_unit_valid[(thisSession['rat']==df_unit_valid['rat']) & (thisSession['session']==df_unit_valid['session'])]
    df_session_list['mean_RDI_LScene'][index] = np.nanmean(df['RDI_LScene'])
    df_session_list['std_RDI_LScene'][index] = np.nanstd(df['RDI_LScene'])
    df_session_list['mean_RDI_RScene'][index] = np.nanmean(df['RDI_RScene'])
    df_session_list['std_RDI_RScene'][index] = np.nanstd(df['RDI_RScene'])
    df_session_list['mean_RDI_LR'][index] = np.nanmean(df['RDI_LR'])
    df_session_list['std_RDI_LR'][index] = np.nanstd(df['RDI_LR'])

    
#%%
df_rip_valid['mRDI_L']=np.nan
df_rip_valid['mRDI_R']=np.nan
df_rip_valid['mRDI_C']=np.nan

df_rip_valid['zRDI_L']=np.nan
df_rip_valid['zRDI_R']=np.nan
df_rip_valid['zRDI_C']=np.nan

df_rip_valid['tRDI_L']=np.nan
df_rip_valid['tRDI_R']=np.nan
df_rip_valid['tRDI_C']=np.nan

df_rip_valid['pRDI_L']=np.nan
df_rip_valid['pRDI_R']=np.nan
df_rip_valid['pRDI_C']=np.nan

for index, thisRip in df_rip_valid.iterrows():
    
    df=df_act_valid[(thisRip['ID']==df_act_valid['RippleID'])]
    df_unit_ripple = df_unit_valid[df_unit_valid['ID'].isin(df['UnitID'])]
    thisSession = df_session_list[(thisRip['rat']==df_session_list['rat']) & (thisRip['session']==df_session_list['session'])]
    df_unit_session = df_unit_valid[(df_unit_valid['rat'].isin(thisSession['rat'])) & (df_unit_valid['session'].isin(thisSession['session']))]

    df_rip_valid['mRDI_L'][index] = np.nanmean(df_unit_ripple['RDI_LScene'])
    df_rip_valid['mRDI_R'][index] = np.nanmean(df_unit_ripple['RDI_RScene'])
    df_rip_valid['mRDI_C'][index] = np.nanmean(df_unit_ripple['RDI_LR'])
    
    df_rip_valid['zRDI_L'][index] = (df_rip_valid['mRDI_L'][index]-thisSession['mean_RDI_LScene'])/thisSession['std_RDI_LScene']
    df_rip_valid['zRDI_R'][index] = (df_rip_valid['mRDI_R'][index]-thisSession['mean_RDI_RScene'])/thisSession['std_RDI_RScene']
    df_rip_valid['zRDI_C'][index] = (df_rip_valid['mRDI_C'][index]-thisSession['mean_RDI_LR'])/thisSession['std_RDI_LR']

    [df_rip_valid['tRDI_L'][index],df_rip_valid['pRDI_L'][index]] = stats.ttest_ind(df_unit_ripple['RDI_LScene'], df_unit_session['RDI_LScene'], equal_var=True, alternative='two-sided',nan_policy='omit')
    [df_rip_valid['tRDI_R'][index],df_rip_valid['pRDI_R'][index]] = stats.ttest_ind(df_unit_ripple['RDI_RScene'], df_unit_session['RDI_RScene'], equal_var=True, alternative='two-sided',nan_policy='omit')
    [df_rip_valid['tRDI_C'][index],df_rip_valid['pRDI_C'][index]] = stats.ttest_ind(df_unit_ripple['RDI_LR'], df_unit_session['RDI_LR'], equal_var=True, alternative='two-sided',nan_policy='omit')
    
df_rip_valid.to_excel(f'{ROOT_data}/RipplesTable_{thisRegion}_ForAnalysis.xlsx')    
    
