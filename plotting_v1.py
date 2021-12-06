import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df = pd.read_excel('D:/HPC-SWR project/Processed Data/UnitsTable_all.xlsx')
df[df.Type>0]
plt.scatter(x=df.PeakPos,
            y=df.Rat/(561*200),
            s=df.N_PeakPos*20,
            alpha=0.6,
            edgecolors='white')
plt.plot([31, 31], [0, 0.06], 'k--', lw=1)
sns.distplot(df.PeakPos,kde=False, norm_hist=True)


plt.ylim(0,0.06)
plt.show()
