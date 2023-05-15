# -*- coding: utf-8 -*-
"""
Created on Mon May 15 20:57:23 2023

@author: Jaemin
"""

import pandas as pd
from scipy.stats import chi2_contingency

swr = pd.DataFrame({'NS' : [91,19], 'N' : [203,43]})
swr = pd.DataFrame({'NS' : [177,79], 'N' : [592,134]})
swr = pd.DataFrame({'NS' : [177-91,79-19], 'N' : [592-203,134-43]})
swr.index = ['NR','R']

chiresult = chi2_contingency(swr, correction=False)
chiresult[1]


chi2, p, dof, expected = chi2_contingency(swr)

msg = 'Test Statistic: {}\np-value: {}\nDegree of Freedom: {}'
print(msg.format(chi2, p, dof))
print(expected)
