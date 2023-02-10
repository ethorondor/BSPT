import os
os.chdir("/mnt/MyDoc/Dropbox/Research/MntStrBrk/code")
import multiprocessing as mp
import numpy as np
import pandas as pd
import scipy.stats
from operator import mul
import Bsquid
import timeit
import matplotlib.pyplot as plt
import statsmodels.formula.api as sm
pd.options.display.float_format = '{:.2f}'.format
rho = 0.01
x1 = np.random.normal(0,1,600)
x2 = np.random.normal(0,1,600)
error = np.random.normal(0,1,600)
beta11 = 1
beta12 = 1.8
beta2 = 0.8
brkStr = 110
y = np.concatenate((beta11*x1[:brkStr],beta12*x1[brkStr:])) + beta2*x2 + error
df = pd.DataFrame()
df['x1'] = x1
df['x2'] = x2
df['y'] = y
reg = sm.ols(formula="y~x1+x2",data=df)
reslt = reg.fit()
print(reslt.summary())
pmDef = {'sigma1ul':1.3,
         'sigma1uu':3.0,
         'sigma1ll':0.1,
         'sigma1lu':0.5,
         'sigma1_sp_size':50,
         'rho':0.01,
         'smpl_sp_l': -4,
         'smpl_sp_u':4,
         'smpl_sp_size': 200
        }
pstmf = Bsquid.getPS_TMF_sigma(pmDef)
def getSimu(cn):
    np.random.seed(cn+123786)
    brkStr = 300
    brkSize = 1
    smpl = 900
    y = np.concatenate((np.random.normal(0,1,brkStr),np.random.normal(0,1+brkSize,smpl-brkStr)))
    b ={'runStr' : 99,
        'rho':rho,
        'pi0':[0.5,0.5],
        'p21':0.00,
        'c':0.003,
        'pstmf':pstmf,
        'vec':y}
    r = Bsquid.Bsquid_sigma(b)
    return(r)

##############################################################
#################  start simulation  #########################
##############################################################
start = timeit.default_timer()
nSim = 10000
pool = mp.Pool(10)
result = np.array(pool.map(getSimu,np.arange(nSim)))
pool.close()
pool.join()

df = pd.DataFrame(result,columns=['c','sigma1','pi','piStar','y'])
df.to_csv("../output/results.sigma.300.003.csv", sep=',',index=False)
stop = timeit.default_timer()
print('time ', stop - start)
