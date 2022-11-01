# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 17:00:11 2022
3. Activated cells의 co-activity rate 계산
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
from matplotlib.colors import ListedColormap
import re
#%% import data
ROOT_data = 'D:/HPC-SWR project/Processed Data'
thisRID=561
thisSID='all'
thisRegion='CA1'
   
df_rip_valid = pd.read_excel(f'{ROOT_data}/RipplesTable_{thisRegion}_ForAnalysis.xlsx')
df_unit_valid = pd.read_excel(f'{ROOT_data}/UnitsTable_{thisRegion}_ForAnalysis.xlsx')
df_act_valid = pd.read_excel(f'{ROOT_data}/ReactTable_{thisRegion}.xlsx')


#%%
pivot_overlap_dict={}
        
    
    df_act_comb = pd.DataFrame(columns=[0,1])
    for thisRipID in df_act_valid['RippleID'].unique():
        temp = df_act_valid[df_act_valid['RippleID']==thisRipID] 
        if len(temp)>1:
            for subset in ite.combinations(temp['UnitID'],2):       
                df_act_comb=df_act_comb.append(pd.DataFrame(list(subset)).T)
    df_act_comb.columns=['u1','u2']
    
    df_unit_comb = pd.DataFrame(columns=[0,1])
    sessions = df_rip_valid['rat'].astype(str).str.cat(df_rip_valid['session'].astype(str),sep='-').unique()

    for thisRSID in sessions:
        thisRSIDn = re.findall(r'\d+',thisRSID)
        df_unit_now = df_unit_valid[(df_unit_valid.rat==int(thisRSIDn[0])) & (df_unit_valid.session==int(thisRSIDn[1]))]
        for subset in ite.combinations(df_unit_now['ID'],2):       
            df_unit_comb=df_unit_comb.append(pd.DataFrame(list(subset)).T)
            
    df_unit_comb.reset_index(drop=True,inplace=True)
    df_unit_comb.columns=['u1','u2']
    df_unit_comb['n']=0
    df_unit_comb['n1']=0
    df_unit_comb['n2']=0
    df_unit_comb['n12']=0
    df_unit_comb['rl1']=0
    df_unit_comb['rr1']=0
    df_unit_comb['rc1']=0
    df_unit_comb['rl2']=0
    df_unit_comb['rr2']=0
    df_unit_comb['rc2']=0
    for index,UnitComb in df_unit_comb.iterrows():
        thisSID = re.findall(r'\d+',UnitComb[1])
        df_unit_comb.n[index]  = len(df_rip_valid[(df_rip_valid.rat==int(thisSID[0])) & (df_rip_valid.session==int(thisSID[1]))])
        df_unit_comb.n1[index]=len(df_act_valid[UnitComb.u1==df_act_valid.UnitID])
        df_unit_comb.n2[index]=len(df_act_valid[UnitComb.u2==df_act_valid.UnitID])
        df_unit_comb.n12[index]=len(df_act_comb[(UnitComb.u2==df_act_comb.u1) & (UnitComb.u1==df_act_comb.u2)])+len(df_act_comb[(UnitComb.u2==df_act_comb.u2) & (UnitComb.u1==df_act_comb.u1)])
        df_unit_comb.rl1[index] = df_unit_valid[df_unit_valid.ID==UnitComb.u1]['RDI_LScene']
        df_unit_comb.rr1[index] = df_unit_valid[df_unit_valid.ID==UnitComb.u1]['RDI_RScene']
        df_unit_comb.rc1[index] = df_unit_valid[df_unit_valid.ID==UnitComb.u1]['RDI_LR']
        df_unit_comb.rl2[index] = df_unit_valid[df_unit_valid.ID==UnitComb.u2]['RDI_LScene']
        df_unit_comb.rr2[index] = df_unit_valid[df_unit_valid.ID==UnitComb.u2]['RDI_RScene']
        df_unit_comb.rc2[index] = df_unit_valid[df_unit_valid.ID==UnitComb.u2]['RDI_LR']
    
    
    df_unit_comb['p1']=df_unit_comb.n1/df_unit_comb.n
    df_unit_comb['p2']=df_unit_comb.n2/df_unit_comb.n
    df_unit_comb['p12']=df_unit_comb.n12/df_unit_comb.n
   
    
   
    temp = df_unit_comb[(df_unit_comb.rl1>0) & (df_unit_comb.rl2>0)]
    temp2 = df_unit_comb[(df_unit_comb.rl1<0) & (df_unit_comb.rl2>0)]
    temp3 = df_unit_comb[(df_unit_comb.rl1<0) & (df_unit_comb.rl2<0)]

fig, ax=plt.subplots()
ax.boxplot([temp.p12,temp2.p12,temp3.p12])
    
