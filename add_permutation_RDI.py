# -*- coding: utf-8 -*-
"""
Created on Wed Dec 21 18:19:51 2022

@author: Jaemin
"""

import random
df_rip_valid[f'per{thisP}']=np.nan

for index,thisRip in df_rip_valid.iterrows():
    temp=np.empty((1000,1))
    
    thisUnits =df_unit_valid[ df_unit_valid['ID'].isin(df_act_valid[df_act_valid['RippleID']==thisRip.ID]['UnitID'])]
    thisPool = df_unit_valid[(df_unit_valid['rat']==thisRip['rat']) & (df_unit_valid['session']==thisRip['session'])]
    thisPool = thisPool[~np.isnan(thisPool[thisParm])]
    for i in range(1000):
        samp = random.sample(thisPool[thisParm].to_list(),len(thisUnits))
        temp[i,0]=np.mean(samp)
    
    m= np.nanmean(thisUnits[thisParm])
    df_rip_valid[f'per{thisP}'][index] = min(len(temp[temp>m]),len(temp[temp<m])) / len(temp)
        
    
