# -*- coding: utf-8 -*-
"""
Created on Mon Mar 14 17:10:51 2022

@author: Jaemin
"""

#%% import libraries
import numpy as np 
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import Mod_SWR as swr
from sklearn.linear_model import LinearRegression
import statsmodels.api as sm

#%% import data
ROOT_data = 'D:/HPC-SWR project/Processed Data'
thisRegion='CA1'

df_rip = pd.read_excel(f'{ROOT_data}/RipplesTable_Ensemble_CA1_10s.xlsx')
df_unit = pd.read_excel(f'{ROOT_data}/UnitsTable_RPR_CA1.xlsx')
df_act = pd.read_excel(f'{ROOT_data}/ReactTable_CA1_10s.xlsx')
df_clust_summ = pd.read_excel('D:/HPC-SWR project/Information Sheet/ClusterSummary_SWR.xlsx')

df_unit['RPR_ZB'] = df_unit['RipPartRate_Z'] - df_unit['RipPartRate_B']
df_unit['RPR_PM'] = df_unit['RipPartRate_P'] - df_unit['RipPartRate_M']
df_unit['RPR_LR'] = df_unit['RipPartRate_L'] - df_unit['RipPartRate_R']

#%%
crit_RPR = 0.05
crit_RDI = 0
boot_corrcoef = pd.DataFrame()

df_unit_valid = df_unit[(df_unit.RipPartRate_all>=crit_RPR)]


for i in range(10000):
    boot_sample = df_unit_valid.sample(200,replace = True)
    
    
    list1 = ['ZB','PM','LR']
    list2 = ['Zebra','Bamboo','Pebbles','Mountains','Left','Right']
    
    temp_boot = list()
    for i in range(3):
        
        boot_sample = boot_sample[(boot_sample[f'RDI_{list1[i]}'].abs()>=crit_RDI)]
        
        X = boot_sample[f'RDI_{list1[i]}']
        y = boot_sample[f'RPR_{list1[i]}']
        
        y[y.isna()] = 0
        X[X.isna()] = 0
        
        
        line_fitter = LinearRegression()
        line_fitter.fit(X.values.reshape(-1,1), y)
    
        X2 = sm.add_constant(X)
        est = sm.OLS(y, X2)
        est2 = est.fit()
        
        temp_boot.append(np.corrcoef(X,y)[0,1])
        
    
    boot_corrcoef = boot_corrcoef.append(pd.DataFrame(temp_boot).T,ignore_index = True)
    
sns.distplot(boot_corrcoef.iloc[:,2],bins=60)
    