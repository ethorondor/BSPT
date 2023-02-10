# %%
import pylab as p
import PIL
from matplotlib import pyplot
import matplotlib.pyplot as plt
import datetime
from statsmodels.tsa.arima_model import ARIMA
import statsmodels.formula.api as sm
import scipy.stats
import pandas as pd
import numpy as np
import matplotlib.image as mpimg
import matplotlib.cm as cm
from PIL import Image
import os
os.chdir("/home/elrond/Documents/Dropbox/Research/MntStrBrk")
pd.options.display.float_format = '{:.2f}'.format


def crt_df(hpi, state):
    df = hpi.loc[hpi.state == state]
    df.loc[:, 'date'] = pd.date_range('1991-01-01', periods=113, freq='Q')
    # annulized growth rate
    df.loc[:, 'annGrwRt'] = (1+df.index_sa.pct_change())**4-1
    df.set_index('date', inplace=True)
    return(df)


def crtDataFrame(ca, ma, oh, tx, fl):
    df = pd.DataFrame(index=ca.index)
    df['ca'] = ca['annGrwRt']
    df['oh'] = oh['annGrwRt']
    df['ma'] = ma['annGrwRt']
    df['fl'] = fl['annGrwRt']
    df['tx'] = tx['annGrwRt']
    return(df)


hpi = pd.read_csv("./data/HPI_PO_state.csv")
st_tab = pd.crosstab(index=hpi.state, columns='count')
ca = crt_df(hpi, 'CA')
ma = crt_df(hpi, 'MA')
oh = crt_df(hpi, 'OH')
tx = crt_df(hpi, 'TX')
fl = crt_df(hpi, 'FL')
df = crtDataFrame(ca, ma, oh, tx, fl)
fig, ax = plt.subplots(figsize=(10, 5.15))
plt.style.use('ggplot')
ax.plot(ca.index, ca.annGrwRt, linewidth=2, label='CA')
ax.plot(ma.index, ma.annGrwRt, '--', linewidth=2, label='MA')
ax.plot(oh.index, oh.annGrwRt, linewidth=2, label='OH')
ax.plot(fl.index, fl.annGrwRt, '-.', linewidth=2, label='FL')
ax.plot(tx.index, tx.annGrwRt, linewidth=2, label='TX')
ax.set_ylabel('Annual Growth Rate')
ax.legend(loc='best')
fig.savefig('HPI_5states.png')
# %%
figprops = dict(figsize=(10, 5.15))
fig1 = p.figure(**figprops)
adjustprops = dict()
image = Image.open('HPI_5states.png').convert("L")
arr = np.asarray(image)
p.figimage(arr, cmap=cm.Greys_r)
p.savefig('./document/graph/grayed_HPI_5States.png')
# %%
