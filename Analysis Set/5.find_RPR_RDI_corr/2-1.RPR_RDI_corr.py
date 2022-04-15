# -*- coding: utf-8 -*-
"""
Created on Wed Feb  9 14:01:42 2022
2-1. unit RDI와 ripple part rate difference의 상관관계 분
@author: JM_Seol
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

df_unit['RPR_ZB'] = (df_unit['RipPartRate_Z'] - df_unit['RipPartRate_B']) / pd.concat([df_unit['RipPartRate_Z'],df_unit['RipPartRate_B']],axis=1).max(axis=1)
df_unit['RPR_PM'] = (df_unit['RipPartRate_P'] - df_unit['RipPartRate_M'])/ pd.concat([df_unit['RipPartRate_P'],df_unit['RipPartRate_M']],axis=1).max(axis=1)
df_unit['RPR_LR'] = (df_unit['RipPartRate_L'] - df_unit['RipPartRate_R']) / pd.concat([df_unit['RipPartRate_L'],df_unit['RipPartRate_R']],axis=1).max(axis=1)


#%%
crit_RPR = [0.2,0.15,0.1,0.05,0]
crit_RDI = [0.2,0.15,0.1,0.05,0]

for c1 in crit_RPR:
    for c2 in crit_RDI:
        df_unit_valid = df_unit[(df_unit.RipPartRate_all>=c1)]



        f,axes=plt.subplots(1,3,figsize=(17,8))
        plt.subplots_adjust(left=0.05, bottom=0.1, right=0.95, top=0.9, wspace=0.4, hspace=0.2)
        plt.suptitle(f"rat 234-561 cells'  RDI vs. Diff.Ripple Part Rate_{c2}_RDI_{c1}_RPRall")
        
        list1 = ['ZB','PM','LR']
        list2 = ['Zebra','Bamboo','Pebbles','Mountains','Left','Right']
        
        for i in range(3):
            df_unit_valid_temp = df_unit_valid[(df_unit_valid[f'RDI_{list1[i]}'].abs()>=c2)]
            
            X = df_unit_valid_temp[f'RDI_{list1[i]}']
            y = df_unit_valid_temp[f'RPR_{list1[i]}'] 
            
            y[y.isna()] = 0
            X[X.isna()] = 0
            
            
            line_fitter = LinearRegression()
            line_fitter.fit(X.values.reshape(-1,1), y)
            
            X2 = sm.add_constant(X)
            est = sm.OLS(y, X2)
            est2 = est.fit()
            
        
            
            axes[i].scatter( x=X,y=y)
            axes[i].plot(X,line_fitter.predict(X.values.reshape(-1,1)),color='k')
            axes[i].annotate(f'CorrCoef={np.corrcoef(X,y)[0, 1].round(3)}',(0,y.max()*0.95))
            axes[i].set_xlabel(f'Rate Diff. Index ({list2[i*2]} vs. {list2[i*2+1]})')
            axes[i].set_ylabel(f'Diff. Ripple Part. Rate ({list2[i*2]} - {list2[i*2+1]}), Normalized by max')
            plt.rc('font', size=15)
        
            plt.style.use('seaborn-whitegrid')
            
            plt.savefig(f'{ROOT_data}/plots/diff_normalized/10s_ITI/corr_rdi{c2}_rpr{c1}.png', dpi=300)
            
    # plt.axis([0, 1.25, 0, 0.3])
    
plt.close('all')  
    