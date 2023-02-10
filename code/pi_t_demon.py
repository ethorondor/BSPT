# %%
import timeit
import Bsquid
from operator import mul
import scipy.stats
import pandas as pd
import numpy as np
import multiprocessing as mp
import os
os.chdir("/home/elrond/Dropbox/Research/MntStrBrk/code")
pd.options.display.float_format = '{:.2f}'.format

y = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1]
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
b = {'runStr': 5,
     'rho': rho,
     'pi0': [0.5, 0.5],
     'p21': 0.00,
     'c': 0.01,
     'pstmf': pstmf,
     'vec': y}
r = Bsquid.Bsquid(b)
r_graph = Bsquid.Bsquid_crtGraph(b)
# %%
