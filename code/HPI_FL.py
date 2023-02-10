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
os.chdir("/home/elrond/Dropbox/Research/MntStrBrk")
pd.options.display.float_format = '{:.2f}'.format
# %%


def getResid(hpi):
    st_tab = pd.crosstab(index=hpi.state, columns='count')
    df = hpi.loc[hpi.state == 'FL']
    df.loc[:, 'date'] = pd.date_range('1991-01-01', periods=113, freq='Q')
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


def getResidBust(hpi):
    df = hpi.loc[hpi.state == 'FL']
    df.loc[:, 'date'] = pd.date_range('1991-01-01', periods=113, freq='Q')
    df.loc[:, 'annGrwRt'] = (1+df.index_sa.pct_change())**4-1
    df.set_index('date', inplace=True)
    df['lag1'] = df.annGrwRt.shift(1)
    df['lag2'] = df.annGrwRt.shift(2)
    df_his = df.loc['1991-06-30':'2005-12-31',
                    ['state', 'annGrwRt', 'lag1', 'lag2']]
    model = sm.ols(formula='annGrwRt~lag1+lag2', data=df_his).fit()
    df['pred'] = model.predict(df)
    df['fst_error'] = (df['annGrwRt'] - df['pred'])/np.std(model.resid)
    resid = df.loc['2004-03-31':, 'fst_error']
    return(resid)


def plotResultBoom(r1graph):
    r1graph['date'] = resid.index[1:]
    r1graph.set_index('date', inplace=True)
    r1graph = r1graph.drop(r1graph.index[0])
    r1graph = r1graph.loc['2000-06-30':'2006-12-31', :]
    fig, (ax1, ax2) = plt.subplots(2, figsize=(10, 5.15))
    plt.style.use('ggplot')
    ax1.plot(r1graph.index, r1graph.y, 'b-',
             label="Forecast Errors", linewidth=2, color='black')
    ax1.legend(loc='best')
    ax1.set_ylabel('Forecast Error')
    ax1.axvline(dt.datetime(2004, 3, 31), linewidth=2, color='black')
    ax1.text(dt.datetime(2004, 3, 31), 0.6, "2004Q1")
    ax2.plot(r1graph.index, r1graph.pi, '--',
             label="Detecting Statistics:$\pi$", linewidth=2, color="black")
    ax2.plot(r1graph.index, r1graph.piStar, 'b-',
             label="Sequential Stopping Procedure:$\mathcal{T}^{*}$", linewidth=2, color='black')
    ax2.axvline(dt.datetime(2004, 3, 31), linewidth=2, color='black')
    ax2.legend(loc='best')
    ax2.set_ylim(0, 1.1)
    ax2.set_ylabel('Probability')
    ax2.text(dt.datetime(2004, 3, 31), 0.6, "2004Q1")
    return(fig)


def plotResultBust(r2graph, resid):
    r2graph['date'] = resid.index[1:]
    r2graph.set_index('date', inplace=True)
    r2graph = r2graph.drop(r2graph.index[0])
    r2graph = r2graph.loc['2004-03-31':'2009-12-31', :]
    fig, (ax1, ax2) = plt.subplots(2, figsize=(10, 5.15))
    plt.style.use('ggplot')
    ax1.plot(r2graph.index, r2graph.y, 'b-',
             label="Forecast Errors", linewidth=2, color='black')
    ax1.legend(loc='best')
    ax1.set_ylabel('Forecast Error')
    ax1.axvline(dt.datetime(2006, 3, 31), linewidth=2, color='black')
    ax1.text(dt.datetime(2006, 3, 31), 0.6, "2006Q1")
    ax2.plot(r2graph.index, r2graph.pi, '--',
             label="Detecting Statistics:$\pi$", linewidth=2, color="black")
    ax2.plot(r2graph.index, r2graph.piStar, 'b-',
             label="Sequential Stopping Procedure:$\mathcal{T}^{*}$", linewidth=2, color='black')
    ax2.axvline(dt.datetime(2006, 3, 31), linewidth=2, color='black')
    ax2.legend(loc='best')
    ax2.set_ylim(0, 1.1)
    ax2.set_ylabel('Probability')
    ax2.text(dt.datetime(2006, 3, 31), 0.6, "2006Q1")
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
    b = {'runStr': 1,
         'rho': 0.01,
         'pi0': [0.5, 0.5],
         'p21': 0.00,
         'c': 0.01,
         'pstmf': pstmf,
         'vec': resid.values}
    r1 = Bsquid.Bsquid(b)
    r1graph = Bsquid.Bsquid_crtGraph(b)
    fig1 = plotResultBoom(r1graph)
    fig1.savefig("./document/graph/Bsquid_FL_boom.jpg")
    resid_bust = getResidBust(hpi)
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
    b = {'runStr': 1,
         'rho': 0.01,
         'pi0': [0.5, 0.5],
         'p21': 0.00,
         'c': 0.01,
         'pstmf': pstmf,
         'vec': resid_bust.values}
    r2 = Bsquid.Bsquid(b)
    r2graph = Bsquid.Bsquid_crtGraph(b)
    fig2 = plotResultBust(r2graph, resid_bust)
    fig2.savefig('./document/graph/Bsquid_FL_bust.jpg')
    print("the main is done")


# %%
