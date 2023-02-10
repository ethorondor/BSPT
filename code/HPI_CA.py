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
os.chdir("/home/thorondor/Dropbox/Research/MntStrBrk")
pd.options.display.float_format = '{:.2f}'.format
#%%
'''
We invesgate the boom and bust of CA HPI. We study the forecast error, did not find the evidence of boom.
but we did find evidence of bust.
'''
def getResidual(hpi):
    st_tab = pd.crosstab(index=hpi.state, columns='count')
    df = hpi.loc[hpi.state == 'CA']
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


def plotResult(r_graph):
    r_graph['date'] = resid.index[1:]
    r_graph.set_index('date', inplace=True)
    r_graph = r_graph.drop(r_graph.index[0])
    r_graph = r_graph.loc['2000-06-30':'2009-12-31', :]
    fig, (ax1, ax2) = plt.subplots(2, figsize=(10, 5.15))
    plt.style.use('ggplot')
    ax1.plot(r_graph.index, r_graph.y, 'b-',
             label="Forecast Errors", linewidth=2, color='black')
    ax1.legend(loc='best')
    ax1.set_ylabel('Forecast Error')
    ax1.axvline(dt.datetime(2005, 9, 30), linewidth=2, color='black')
    ax1.text(dt.datetime(2005, 9, 30), 0.6, "2005Q3")
    ax2.plot(r_graph.index, r_graph.pi, '--',
             label="Detecting Statistics:$\pi$", linewidth=2, color="black")
    ax2.plot(r_graph.index, r_graph.piStar, 'b-',
             label="Sequential Stopping Procedure:$\mathcal{T}^{*}$", linewidth=2, color='black')
    ax2.axvline(dt.datetime(2005, 9, 30), linewidth=2, color='black')
    ax2.legend(loc='best')
    ax2.set_ylabel('Probability')
    ax2.text(dt.datetime(2005, 9, 30), 0.8, "2005Q3")
    return(fig)


if __name__ == "__main__":
    hpi = pd.read_csv("./data/HPI_PO_state.csv")
    resid = getResidual(hpi)
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
    r = Bsquid.Bsquid(b)
    r_graph = Bsquid.Bsquid_crtGraph(b)
    resid.index[23]
    fig = plotResult(r_graph)
    fig.savefig('./document/graph/Bsquid_CA_bust.jpg')
    print("the main is done")


# %%
