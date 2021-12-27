# -*- coding: utf-8 -*-
"""
Created on Mon Dec 27 12:00:01 2021
Modules for SWR project

@author: JM_Seol
"""

def Exp_from_ClustSum(df_clust_summ,thisRID,df_unit, props):
    temp = str(thisRID) + '-0' + df_unit['TT-Unit']    
    s = df_clust_summ.loc[df_clust_summ['cellID'].isin(temp),props]    
    return s
    
def DrawDist_2samp(Signed,*dat,**opts):
    import seaborn as sns
    import matplotlib.pyplot as plt
    import scipy as sp
    f,axes=plt.subplots(1,2,figsize=(14,8))
    
    for i,d in enumerate(dat):
        axes[0] = sns.histplot(d,stat='probability',ax=axes[0],color=opts['cmap'][i],
                               bins=opts['bins'],binrange=opts['range'])
        
        axes[1] = sns.ecdfplot(d,ax=axes[1],color=opts['cmap'][i])

    
    axes[0].legend(opts['legend'],fontsize=10)
    axes[0].set_title(opts['title_hist'])
    plt.rc('font',size=15)
    
    if 'xlab' in opts:
        axes[0].set_xlabel(opts['xlab'])
        axes[1].set_xlabel(opts['xlab'])
    if len(dat)==2:
        p=sp.stats.ks_2samp(dat[0],dat[1])
        # plt.text(0.9,1.05, f'p={round(p[1],3)}',transform=axes[1].transAxes)
    axes[1].set_title(opts['title_dist'] + f'p={round(p[1],3)})')
    axes[1].legend(opts['legend'],fontsize=10)
    
    if Signed:
        for i,d in enumerate(dat):
            p=sp.stats.wilcoxon(d)
            axes[0].text(0.9,1.05+0.03*i,f'p={round(p[1],3)}',color=opts['cmap'][i],transform=axes[0].transAxes)
            
    return axes
