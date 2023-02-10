# %%
import timeit
import Bsquid
from operator import mul
import scipy.stats
import pandas as pd
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
from matplotlib import pyplot
import os
os.chdir("/home/elrond/Dropbox/Research/MntStrBrk/code")
pd.options.display.float_format = '{:.2f}'.format

rho = 0.01

pmDef = {'mu1_p_l': 0.5,
         'mu1_p_u': 2.0,
         'mu1_n_l': -2,
         'mu1_n_u': -0.5,
         'mu1_sp_size': 50,
         'rho': rho,
         'smpl_sp_l': -4,
         'smpl_sp_u': 4,
         'smpl_sp_size': 200
         }
pstmf = Bsquid.getPS_TMF(pmDef)

simParam = [[110, 0.8, 0],
            [150, 0.8, 0],
            [200, 0.8, 0],
            [300, 0.8, 0],
            [110, 0.3, 0.5],
            [150, 0.3, 0.5],
            [200, 0.3, 0.5],
            [300, 0.3, 0.5],
            [110, 0.1, 0.7],
            [150, 0.1, 0.7],
            [200, 0.1, 0.7],
            [300, 0.1, 0.7],
            [110, 0.0, 0.8],
            [150, 0.0, 0.8],
            [200, 0.0, 0.8],
            [300, 0.0, 0.8]]


def getSimu(cn, tau, alpha, gamma):
    smpl = 900
    np.random.seed(cn+3217865)
    y = np.empty(shape=[smpl, ])
    for i in range(smpl):
        if i < tau:
            y[i] = np.random.normal(0, 1, 1)
        if i >= tau:
            y[i] = alpha + gamma*y[i-1] + np.random.normal(0, 1, 1)
    b = {'runStr': 99,
         'rho': rho,
         'pi0': [0.5, 0.5],
         'p21': 0.00,
         'c': 0.008,
         'pstmf': pstmf,
         'vec': y}
    r = Bsquid.Bsquid(b)
    return(r)


def getBSPT(pm):
    tau = pm[0]
    alpha = pm[1]
    gamma = pm[2]
    nSim = 5000
    df = pd.DataFrame(columns=['c', 'mu1', 'pi', 'piStar', 'y'])
    for i in range(nSim):
        df.loc[i] = getSimu(i, tau, alpha, gamma)
    #df = pd.DataFrame(result, columns=['c', 'mu1', 'pi', 'piStar', 'y'])
    d = df['c'].dropna()
    typeII = (nSim-d.size)/nSim
    typeI = d[d < tau].size/nSim
    ed = np.mean(d[d >= tau]) - tau
    ed_sd = np.std(d[d >= tau])
    return([tau, alpha, gamma, ed, ed_sd, typeI, typeII])


result = pd.DataFrame(
    columns=['tau', 'alpha', 'gamma', 'expected depay', 'sd', 'type I', 'type II'])
result.loc[0] = getBSPT(simParam[13])
# for p in range(len(simParam)):
#    result.loc[p] = getBSPT(simParam[p])

print("done")

# %%
