# %%
from matplotlib import pyplot
import matplotlib.pyplot as plt
import datetime
import Bsquid
from statsmodels.tsa.arima_model import ARIMA
import statsmodels.formula.api as sm
import scipy.stats
import pandas as pd
import numpy as np
import datetime as dt
import os
os.chdir("/home/elrond/Documents/Dropbox/Research/MntStrBrk")
pd.options.display.float_format = '{:.2f}'.format
%matplotlib inline


def getResid(hpi):
    st_tab = pd.crosstab(index=hpi.state, columns='count')
    df = hpi.loc[hpi.state == 'TX']
    df.loc[:, 'date'] = pd.date_range('1991-01-01', periods=113, freq='Q')
    # annulized growth rate
    df.loc[:, 'annGrwRt'] = (1+df.index_sa.pct_change())**4-1
    df.set_index('date', inplace=True)
    df['lag1'] = df.annGrwRt.shift(1)
    df['lag2'] = df.annGrwRt.shift(2)
    df_his = df.loc['1991-06-30':'2000-12-31',
                    ['state', 'annGrwRt', 'lag1', 'lag2']]
    model = sm.ols(formula='annGrwRt~lag1+lag2', data=df_his).fit()
    df['pred'] = model.predict(df)
    df['fst_error'] = (df['annGrwRt'] - df['pred'])/np.std(model.resid)
    resid = df.loc['1999-12-31':, 'fst_error']
    return(resid)


def plotResult(r1graph):
    r1graph['date'] = resid.index[4:]
    r1graph.set_index('date', inplace=True)
    r1graph = r1graph.drop(r1graph.index[0])
    fig, (ax1, ax2) = plt.subplots(2, figsize=(10, 5.15))
    plt.style.use('ggplot')
    ax1.plot(r1graph.index, r1graph.y, 'b-',
             label="Forecast Errors", linewidth=2, color='k')
    ax1.legend(loc='best')
    ax1.set_ylabel('Forecast Error')
    ax1.axvline(dt.datetime(2010, 12, 31), linewidth=2, color='k')
    ax1.text(dt.datetime(2010, 12, 31), -2, "2010Q4")
    ax2.plot(r1graph.index, r1graph.pi, '--',
             label="Detecting Statistics:$\pi$", linewidth=2, color="k")
    ax2.plot(r1graph.index, r1graph.piStar, 'b-',
             label="Sequential Stopping Procedure:$\mathcal{T}^{*}$", linewidth=2, color='k')
    ax2.axvline(dt.datetime(2010, 12, 31), linewidth=2, color='black')
    ax2.legend(loc='best')
    ax2.set_ylabel('Probability')
    ax2.text(dt.datetime(2010, 12, 31), 0.1, "2010Q4")
    return(fig)


if __name__ == "__main__":
    hpi = pd.read_csv("./data/HPI_PO_state.csv")
    resid = getResid(hpi)
    pmDef = {'mu1_p_l': 0.6,
             'mu1_p_u': 2.0,
             'mu1_n_l': -2,
             'mu1_n_u': -0.6,
             'mu1_sp_size': 50,
             'rho': 0.01,
             'smpl_sp_l': -4,
             'smpl_sp_u': 4,
             'smpl_sp_size': 200
             }
    pstmf = Bsquid.getPS_TMF(pmDef)
    b1 = {'runStr': 4,
          'rho': 0.01,
          'pi0': [0.5, 0.5],
          'p21': 0.00,
          'c': 0.01,
          'pstmf': pstmf,
          'vec': resid.values}
    r1 = Bsquid.Bsquid(b1)
    r1graph = Bsquid.Bsquid_crtGraph(b1)
    fig = plotResult(r1graph)
    fig.savefig('./document/graph/Bsquid_TX.jpg')
