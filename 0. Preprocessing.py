# -*- coding: utf-8 -*-
"""
Created on Thu Dec 23 15:11:44 2021

0. Ripple & Active Unit Data Preprocessing
input: RipplesTable(session), ReactTable(session), UnitsTable(Session)
output: RipplesTable_all, ReactTable_all, UnitsTable_all
        RipplesTable_valid, ReactTable_valid, UnitsTable_valid
@author: JM_Seol
"""
import numpy as np
import pandas as pd
#%% import data
ROOT_data = 'D:/HPC-SWR project/Processed Data'
df_rip=pd.DataFrame()
df_unit=pd.DataFrame()
df_act=pd.DataFrame()

for thisSID in range(5):
    thisRID = 561
    thisSID = thisSID+1
    thisRegion = 'CA1'
    
    df_rip_temp = pd.read_excel(f'{ROOT_data}/RipplesTable_r{thisRID}_s{str(thisSID).zfill(2)}_{thisRegion}.xlsx')
    df_unit_temp = pd.read_excel(f'{ROOT_data}/UnitsTable_r{thisRID}_s{str(thisSID).zfill(2)}_{thisRegion}.xlsx')
    df_act_temp = pd.read_excel(f'{ROOT_data}/ActTable_r{thisRID}_s{str(thisSID).zfill(2)}_{thisRegion}.xlsx')
    
    df_rip=df_rip.append(df_rip_temp, ignore_index=True)
    df_unit=df_unit.append(df_unit_temp, ignore_index=True)
    df_act=df_act.append(df_act_temp, ignore_index=True)
#%% data processing
df_rip['RipID'] = df_rip['RipID'].astype(str).str.zfill(4)
df_rip['Session'] = df_rip['Session'].astype(str).str.zfill(1)
df_rip['RipID'] = df_rip['Session'] + '-' + df_rip['RipID']

df_unit['Session'] = df_unit['Session'].astype(str).str.zfill(1)
df_unit['TT'] = df_unit['TT'].astype(str).str.zfill(2)
df_unit['Unit'] = df_unit.Unit.astype(str).str.zfill(2)
df_unit['TT-Unit'] = df_unit['Session'] + '-' + df_unit['TT'] +  '-' + df_unit['Unit']

df_act['RipID'] = df_act['RipID'].astype(str).str.zfill(4)
df_act['SessionID'] = df_act['SessionID'].astype(str).str.zfill(1)
df_act['TT'] = df_act['TT'].astype(str).str.zfill(2)
df_act['UnitID'] = df_act.UnitID.astype(str).str.zfill(2)
df_act['TT-Unit'] = df_act['SessionID'] + '-' + df_act['TT'] +  '-' + df_act['UnitID']
df_act['RipID'] = df_act['SessionID'] + '-' + df_act['RipID']

valid_rips = df_rip.RipID[df_rip.Filter>0]
df_act_valid = df_act[df_act['RipID'].isin(valid_rips) & df_act['Type']>0]
df_unit_valid = df_unit[df_unit['TT-Unit'].isin(df_act_valid['TT-Unit'])]


df_rip['meanRDI_ZB'] = 0
df_rip['meanRDI_PM'] = 0
df_rip['meanRDI_LR'] = 0
for i in range(df_rip.shape[0]):
    thisRipID = df_rip.RipID[i]
    thisUnitID = df_act_valid['TT-Unit'][df_act_valid['RipID']==thisRipID]
    thisUnits = df_unit_valid[df_unit_valid['TT-Unit'].isin(thisUnitID)]
    df_rip['meanRDI_ZB'].iloc[i] = np.mean(thisUnits.RDI_ZB)
    df_rip['meanRDI_PM'].iloc[i] = np.mean(thisUnits.RDI_PM)
    df_rip['meanRDI_LR'].iloc[i] = np.mean(thisUnits.RDI_LR)
    
df_rip_valid = df_rip[df_rip['Filter']>0]


df_unit_valid=df_unit_valid.set_axis(df_unit_valid['TT-Unit'],axis=0)
df_rip_valid=df_rip_valid.set_axis(df_rip_valid['RipID'],axis=0)
df_act_valid=df_act_valid[df_act_valid['RipID'].isin(df_rip_valid['RipID'])]
#%% export data

df_unit.to_excel(f'{ROOT_data}/UnitTable_r{thisRID}_all_{thisRegion}.xlsx')
df_rip.to_excel(f'{ROOT_data}/RipplesTable_r{thisRID}_all_{thisRegion}.xlsx')
df_act.to_excel(f'{ROOT_data}/ActTable_r{thisRID}_all_{thisRegion}.xlsx')


df_unit_valid.to_excel(f'{ROOT_data}/UnitTable_r{thisRID}_all_{thisRegion}_v.xlsx')
df_rip_valid.to_excel(f'{ROOT_data}/RipplesTable_r{thisRID}_all_{thisRegion}_v.xlsx')
df_act_valid.to_excel(f'{ROOT_data}/ActTable_r{thisRID}_all_{thisRegion}_v.xlsx')
