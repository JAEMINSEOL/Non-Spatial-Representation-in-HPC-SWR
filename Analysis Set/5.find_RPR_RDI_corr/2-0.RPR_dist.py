# -*- coding: utf-8 -*-
"""
Created on Wed Mar 16 14:58:35 2022

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
crit_RPR = 0.0
crit_RDI = 0

df_unit_valid = df_unit[(df_unit.RipPartRate_all>crit_RPR)]

#%%
sns.histplot(df_unit_valid.RPR_LR,stat='probability',binwidth=0.05, binrange = [-1 ,1])
plt.xlabel('RPR')

#%%
df_unit_valid_temp = df_unit_valid.iloc[101:-1,:]

from scipy.stats import ks_2samp
ks_2samp(df_unit_valid.RipPartRate_all,df_unit_valid_temp.RipPartRate_all)

bw_values =  [None, 0.1, 0.01]

# generate a list of kde estimators for each bw
kde = sp.stats.gaussian_kde(df_unit_valid.RipPartRate_all,bw_method=0.1)


x1=df_unit_valid.RPR_ZB
x1=x1.dropna(how='any')
count, bins_count = np.histogram(x1, bins=100)
pdf = count / sum(count)
cdf = np.cumsum(pdf)
plt.plot(bins_count[1:], cdf, label="CDF")

x2=df_unit_valid.RPR_PM
x2=x2.dropna(how='any')
count, bins_count = np.histogram(x2, bins=100)
pdf = count / sum(count)
cdf = np.cumsum(pdf)
plt.plot(bins_count[1:], cdf, label="CDF")

x3=df_unit_valid.RPR_LR
x3=x3.dropna(how='any')
count, bins_count = np.histogram(x3, bins=100)
pdf = count / sum(count)
cdf = np.cumsum(pdf)
plt.plot(bins_count[1:], cdf, label="CDF")



ks_2samp(x1,x3)