# %%
import timeit
from operator import mul
import scipy.stats
import pandas as pd
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import os
os.chdir("/home/elrond/Dropbox/Research/MntStrBrk/code")
pd.options.display.float_format = '{:.2f}'.format

# %%
opt = pd.read_csv('../output/results.rho10.150.005.csv')
# %%
str = 150
fa = opt.loc[opt.c-str < 0, :]
delay = opt.loc[opt.c-str >= 0, :]
# %%
np.mean(delay.c-str)
# %%
np.std(delay.c-str)


# %%
fa


# %%
delay
# %%
opt

# %%
