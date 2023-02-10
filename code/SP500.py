#%%
import os
os.chdir("/mnt/MyDoc/Dropbox/Research/MntStrBrk/code")
import numpy as np
import pandas as pd
import scipy.stats
import statsmodels.formula.api as sm
import Bsquid
import pyreadr
import matplotlib.pyplot as plt
pd.options.display.float_format = '{:.2f}'.format
'''
we use 2001-08-01 to 2001-09-11 as history period (0 to 28)
'''
sp2001 = pd.read_csv("../data/SP2001.csv")
#%%
pmDef = {'mu1_p_l':0.6,
        'mu1_p_u':2.0,
        'mu1_n_l':-2,
        'mu1_n_u':-0.6,
        'mu1_sp_size':50,
        'rho': 0.01,
        'smpl_sp_l': -4,
        'smpl_sp_u':4,
        'smpl_sp_size': 200
        }
pstmf = Bsquid.getPS_TMF(pmDef)
opt = pd.DataFrame(columns=['stkName','sbk'])
i = 0
for stk in sp2001.columns[1:501]:
    df = sp2001.loc[:,[stk]]
    hist = df.loc[0:28,:]
    hist = hist.dropna()
    mon = df.loc[29:,:]
    mon = mon.dropna()
    resid = (mon[stk]-hist[stk].mean())/hist[stk].std()
    y = resid.values
    b ={'runStr' : 5,
        'rho':0.01,
        'pi0':[0.5,0.5],
        'p21':0.00,
        'c':0.01,
        'pstmf':pstmf,
        'vec':y}
    r = Bsquid.Bsquid(b)
    opt.loc[i]=[stk,r[0]]
    i = i+1
#%%
