# -*- coding: utf-8 -*-
"""
Created on Tue Nov 29 17:58:22 2022

@author: Jaemin
"""



#%%
thisParm='RDI_LScene'
thisP='RDI_L'
Rip_sort = 'RipNum'
Unit_sort = thisParm
for index,Session in df_session_list.iterrows():
    
    index=81
    Session = df_session_list.iloc[index]
    
    df4=df3.copy()
    df4=df4[(Session['rat']==df4['rat_x']) & (Session['session']==df4['session_x'])]
    df4=df4[(df3['nPCs'])>=3]
    
    
    
    df4_r = df4.loc[:,[Rip_sort]].sort_values(by=[Rip_sort],ascending=[True]).apply(tuple, axis=1)
    f, i = pd.factorize(df4_r)
    factorized = pd.Series(f + 1, df4_r.index)

    df4['rank_rip']=factorized
    
    df4_r = df4.loc[:,[Unit_sort]].sort_values(by=[Unit_sort],ascending=[True]).apply(tuple, axis=1)
    f, i = pd.factorize(df4_r)
    factorized = pd.Series(f + 1, df4_r.index)

    df4['rank_unit']=factorized


dfx = pd.DataFrame([x,y])

x = df4['rank_unit']
y = df4['trial'].str[-3:].astype('int')
# y = df4['rank_rip']

bad_indices = np.isnan(x) | np.isnan(y)
good_indices = ~bad_indices
good_x = x[good_indices]
good_y = y[good_indices]

# plt.hist2d(good_x,good_y)
plt.scatter(good_x,good_y)
plt.jet()
# plt.colorbar()
plt.show()

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

