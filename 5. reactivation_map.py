# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 17:58:22 2022

@author: Jaemin
"""



#%%
prefix='m'
thisParm='RDI_LScene'
thisP='RDI_L'
Rip_sort = ['RipNum']
Rip_sort = [f'{prefix}{thisP}','RipNum']
# Rip_sort = [f'StartTime_fromTrialEnd']
Unit_sort = [thisParm]

    
for index in range(len(df_session_list)):
    Session = df_session_list.iloc[index]
    thisSID = f'{Session.rat}-0{Session.session}'
    
    df4=df3.copy()
    df4=df4[~np.isnan(df4[thisParm])]
    df4=df4[(Session['rat']==df4['rat_x']) & (Session['session']==df4['session_x'])]
    df4=df4[((df4['nPCs'])>=3) & (df4['context']==3)]
    
    if len(df4)<=0:
        continue
    
    df4_r = df4.loc[:,Rip_sort].sort_values(by=Rip_sort,ascending=[True for i in range(len(Rip_sort))]).apply(tuple, axis=1)
    f, i = pd.factorize(df4_r)
    factorized = pd.Series(f + 1, df4_r.index)
    
    df4['rank_rip']=factorized
    
    df4_r = df4.loc[:,Unit_sort].sort_values(by=Unit_sort,ascending=[True for i in range(len(Unit_sort))]).apply(tuple, axis=1)
    f, i = pd.factorize(df4_r)
    factorized = pd.Series(f + 1, df4_r.index)
    
    df4['rank_unit']=factorized
    
    
    # dfx = pd.DataFrame([x,y])
    
    y = df4[thisParm]
    # y = df4['trial'].str[-3:].astype('int')
    x = df4['rank_rip']
    # x = df4['trial'].str[7:12].astype(int)
    
    bad_indices = np.isnan(x) | np.isnan(y)
    good_indices = ~bad_indices
    good_x = x[good_indices]
    good_y = y[good_indices]
    
    fig, ax=plt.subplots(figsize=[20, 10])
    h = plt.hist2d(good_x,good_y,bins=[x.max(),50],cmap=mpl.cm.get_cmap('gray_r'))
    plt.plot([0, max(good_x)],[0,0])
    
    plt.colorbar()
    plt.ylim([-1.8, 1.8])
    plt.clim([0, 3])
    plt.rc('font', size=20) 
    plt.title('Overall_Bamboo')
    plt.xlabel(thisParm)
    plt.ylabel('unit')

# plt.savefig(f'{ROOT_data}/plots/Reactivated Heatmap/Overall_Bamboo_all_RDI')
# plt.close('all')

    fig, ax=plt.subplots(figsize=[8, 8])
    plt.hist(good_x)
    plt.rc('font', size=12) 
    plt.title(thisSID)
    plt.xlabel('trial')
    plt.ylabel('# ripple')
    plt.savefig(f'{ROOT_data}/plots/RippleRate/{thisSID}_All')
    plt.close('all')        
#%%
fig, ax=plt.subplots()
plt.scatter(good_x,good_y)
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
thisParm='RDI_RScene'
thisP='RDI_R'

df3 = pd.merge(df,df2, how='left', left_on='UnitID', right_on='ID')

df3 = pd.merge(df3,df1, how='left', left_on='RippleID', right_on='ID')

df4=df3.copy()

df5=df3[(df3['DecodingP_all'])>=0.05]
    df6=df3[(df3['DecodingP_all'])<0.05]
    
    dat1 = df5[f'{thisParm}']
    dat2 = df6[f'{thisParm}']
    filtered_dat1 = dat1[~np.isnan(dat1)]
    filtered_dat2 = dat2[~np.isnan(dat2)]
    
    fig,ax=plt.subplots()
    ax.violinplot([filtered_dat1,filtered_dat2])
    ax.boxplot([filtered_dat1,filtered_dat2])
    ax.set_ylim(-2, 2)
    plt.show()
