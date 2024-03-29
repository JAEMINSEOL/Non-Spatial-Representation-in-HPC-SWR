# -*- coding: utf-8 -*-
"""
Created on Tue Jan 11 17:00:11 2022
3. Activated cells의 co-activity rate 계산
@author: JM_Seol
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
#%% make unit combination pair table
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
df_unit_comb['n12I']=0
df_unit_comb['n12U']=0
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
    df_unit_comb.n12I[index]=len(df_act_comb[(UnitComb.u2==df_act_comb.u1) & (UnitComb.u1==df_act_comb.u2)])+len(df_act_comb[(UnitComb.u2==df_act_comb.u2) & (UnitComb.u1==df_act_comb.u1)])
    df_unit_comb.n12U[index]= df_unit_comb.n1[index]+df_unit_comb.n2[index]-df_unit_comb.n12I[index]
    df_unit_comb.rl1[index] = df_unit_valid[df_unit_valid.ID==UnitComb.u1]['RDI_LScene']
    df_unit_comb.rr1[index] = df_unit_valid[df_unit_valid.ID==UnitComb.u1]['RDI_RScene']
    df_unit_comb.rc1[index] = df_unit_valid[df_unit_valid.ID==UnitComb.u1]['RDI_LR']
    df_unit_comb.rl2[index] = df_unit_valid[df_unit_valid.ID==UnitComb.u2]['RDI_LScene']
    df_unit_comb.rr2[index] = df_unit_valid[df_unit_valid.ID==UnitComb.u2]['RDI_RScene']
    df_unit_comb.rc2[index] = df_unit_valid[df_unit_valid.ID==UnitComb.u2]['RDI_LR']


df_unit_comb['p1']=df_unit_comb.n1/df_unit_comb.n
df_unit_comb['p2']=df_unit_comb.n2/df_unit_comb.n
df_unit_comb['p12I']=df_unit_comb.n12I/df_unit_comb.n
df_unit_comb['p12U']=df_unit_comb.n12I/df_unit_comb.n12U 

    
df_unit_comb.to_excel(f'{ROOT_data}/ReactPair.xlsx')

#%% add rpr for each unit
df_unit_valid['NR']=np.nan
df_unit_valid['RPR']=0
for index,Unit in df_unit_valid.iterrows():
    df_unit_valid['NR'][index]=len(df_act_valid[Unit.ID==df_act_valid.UnitID])
    df_unit_valid['RPR'][index] = df_unit_valid['NR'][index]/len(df_rip_valid[(df_rip_valid.rat==Unit.rat) & (df_rip_valid.session==Unit.session)])
#%%

temp = df_unit_comb[(df_unit_comb.rl1>0) & (df_unit_comb.rl2>0)]
temp2 = df_unit_comb[(df_unit_comb.rl1<0) & (df_unit_comb.rl2>0)]
temp3 = df_unit_comb[(df_unit_comb.rl1<0) & (df_unit_comb.rl2<0)]

fig, ax=plt.subplots()

df = pd.DataFrame(temp.p12U)
df['cdf'] = df.rank(method = 'average', pct = True)
plt.plot(df.sort_values('p12U').p12U,df.sort_values('p12U').cdf)

df = pd.DataFrame(temp3.p12U)
df['cdf'] = df.rank(method = 'average', pct = True)
plt.plot(df.sort_values('p12U').p12U,df.sort_values('p12U').cdf)



ax.boxplot([temp.p12U,temp2.p12U,temp3.p12U])

#%% scatterplot with connected line
df=df_unit_valid
df2=df_unit_comb
for index, unitpair in df2.iterrows():
    plt.plot([unitpair.rl1, unitpair.rl2],[unitpair.p1, unitpair.p2],c=cm.gray(1-unitpair.p12U))
    
plt.scatter(df['RDI_LScene'],df['RPR'],c='r')
plt.show()

for index, unitpair in df2.iterrows():
    if unitpair.p12U>0.3:
        plt.plot([unitpair.rl1, unitpair.rl2],[unitpair.rr1, unitpair.rr2],c=cm.jet(unitpair.p12U))
    
plt.scatter(df['RDI_LScene'],df['RDI_RScene'],c=cm.gray(df.RPR))
plt.show()


        
#%% reactivation rasterplot
df=df_act_valid.copy()
df1=df_rip_valid.copy()
df2 = df_unit_valid.copy()

#%% count ensemble scene/choice selectivity for react. rasterplot
df['RipNum']=np.nan
i=1
for index, Units in df.iterrows():
    df['RipNum'][index]=i
    if index<len(df)-1:
        if df.RippleID[index] != df.RippleID[index+1]:
            i=i+1

df1['nPCs'] = np.nan      
df1['npRDI_L'] = np.nan
df1['nnRDI_L'] = np.nan
df1['npRDI_R']= np.nan
df1['nnRDI_R']= np.nan
df1['npRDI_C'] = np.nan
df1['nnRDI_C'] = np.nan

for index, Rips in df1.iterrows():
    df0=pd.merge(df[df['RippleID']==Rips['ID']],df2,how='left', left_on='UnitID', right_on='ID')
    df1['nPCs'][index]=df0['UnitID'].nunique()
    df1['npRDI_L'][index] = sum(df0.RDI_LScene>0)
    df1['nnRDI_L'][index] = sum(df0.RDI_LScene<0)
    df1['npRDI_R'][index] = sum(df0.RDI_RScene>0)
    df1['nnRDI_R'][index] = sum(df0.RDI_RScene<0)
    df1['npRDI_C'][index] = sum(df0.RDI_LR>0)
    df1['nnRDI_C'][index] = sum(df0.RDI_LR<0)


df3 = pd.merge(df,df2, how='left', left_on='UnitID', right_on='ID')

df3 = pd.merge(df3,df1, how='left', left_on='RippleID', right_on='ID')

#%%
thisParm='RDI_LScene'
thisP='RDI_L'
for index,Session in df_session_list.iterrows():
    
    
    df4=df3.copy()
    df4=df4[(Session['rat']==df4['rat_x']) & (Session['session']==df4['session_x'])]
    df4=df4[(df3['nPCs'])>=3]
    
    
    
    df4_r = df4.loc[:,['RipNum']].sort_values(by=['RipNum'],ascending=[True]).apply(tuple, axis=1)
    f, i = pd.factorize(df4_r)
    factorized = pd.Series(f + 1, df4_r.index)

    df4['rank']=factorized


    if not(df4.empty):
        plt.figure(figsize=(6,8))
        plt.scatter(df4[thisParm][df4[thisParm]>0],df4['rank'][df4[thisParm]>0],marker='|',s=5,c='r')
        plt.scatter(df4[thisParm][df4[thisParm]<0],df4['rank'][df4[thisParm]<0],marker='|',s=5,c='b')
        plt.plot([0,0],[0,max(f)],c='k',ls='--')
        plt.xlim([-2, 2])
        plt.xlabel(f'{thisParm}')
        plt.ylim([0,max(f)])
        plt.ylabel('Ripple')
        plt.title(f'{Session.rat} - {Session.session}')
        # plt.savefig(f'{ROOT_data}/plots/Reactivated Ensemble_raw/{thisP}_{Session.rat}-{Session.session}.png')
        # plt.close()

#%%
thisParm='RDI_LR'
prefix='p'

for index in ['Zebra','Pebbles','Bamboo','Mountains']:
    if index in ['Zebra','Pebbles','Bamboo','Mountains']:
        Exper='LSM'
        CxtNum=['Zebra','Pebbles','Bamboo','Mountains'].index(index)+1
    else:
        Exper='JS'
        CxtNum=['Forest','City'].index(index)+1
        
    if thisParm =='RDI_LScene':
        thisP='RDI_L'
        if Exper=='LSM':
            c1=Color_List['Zebra']; c2 = Color_List['Bamboo']; 
        elif Exper=='JS':
            c1=Color_List['Forest']; c2 = Color_List['City'];
    elif thisParm=='RDI_RScene':
        thisP='RDI_R'; c1=Color_List['Pebbles']; c2 = Color_List['Mountains'];
    elif thisParm=='RDI_LR':
        thisP='RDI_C'; c1=Color_List['Left']; c2 = Color_List['Right'];
        
    df4=df3[(Exper==df3['experimenter_x']) & (df3['nPCs']>=3)]
    df4['pp']=df4[f'np{thisP}']/(df4[f'np{thisP}']+df4[f'nn{thisP}'])
    df4_r = df4.loc[:,['pp',f'np{thisP}',f'nn{thisP}','RipNum']].sort_values(by=['pp',f'np{thisP}',f'nn{thisP}','RipNum'],ascending=[True,True,False,False]).apply(tuple, axis=1)
    f, i = pd.factorize(df4_r)
    factorized = pd.Series(f + 1, df4_r.index)

    df4['rank']=factorized


    if not(df4.empty):
        # df4=df4[abs(df4[f't{thisP}'])>1]
        
        df4_r = df4.loc[:,[f'{prefix}{thisP}','RipNum']].sort_values(by=[f'{prefix}{thisP}','RipNum'],ascending=[True,True]).apply(tuple, axis=1)

        f, i = pd.factorize(df4_r)
        factorized = pd.Series(f + 1, df4_r.index)

        df4['rank2']=factorized
    
        
        plt.figure(figsize=(6,8))
        # plt.scatter(df4[thisParm][df4[thisParm]>0],df4['rank2'][df4[thisParm]>0],marker='|',s=5,c='k')
        # plt.scatter(df4[thisParm][df4[thisParm]<0],df4['rank2'][df4[thisParm]<0],marker='|',s=5,c='k')
        plt.plot([0,0],[0,max(f)],c='k',ls='--')
        plt.xlim([0, 1])
        plt.xlabel(f'{thisParm}_{prefix}')
        plt.ylim([0,max(df4['rank2'])+1])
        plt.ylabel('Ripple')
        plt.title(f'{Exper} dataset')
        # plt.scatter(df4['mRDI_C'], df4['rank'])
        plt.plot([0.05, 0.05],[0,max(f)],c='r',ls='--')
        # plt.plot([-1,-1],[0,max(f)],c='r',ls='--')
        
        # plt.scatter(df4[f'm{thisP}'], df4['rank2'],marker='o',s=5,c='b')
        plt.scatter(df4[f'{prefix}{thisP}'], df4['rank2'],marker='o',s=5,c='r')
        # plt.savefig(f'{ROOT_data}/plots/Reactivated Ensemble_raw/{thisP}_{Exper}-{index}.png')
        # plt.close()
        
plt.hist(df4['mRDI_C'])
plt.hist(df4['zRDI_L'])
plt.hist(df4['zRDI_R'])
