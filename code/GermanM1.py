# %%
import matplotlib.pyplot as plt
import pyreadr
import Bsquid
import statsmodels.formula.api as sm
import scipy.stats
import pandas as pd
import numpy as np
import datetime as dt
import os
os.chdir("/home/elrond/Dropbox/Research/MntStrBrk")
pd.options.display.float_format = '{:.2f}'.format

# %%


def getResidule(df):
    GermanM1 = df['GermanM1']
    historyM1 = df['historyM1']
    monitorM1 = df['monitorM1']
    historyM1['ecm_res'] = historyM1['ecm.res']
    GermanM1['ecm_res'] = GermanM1['ecm.res']
    GermanM1['dm_actual'] = GermanM1['dm']
    historyM1['t'] = pd.date_range('1961-3', '1990-8', freq='Q')
    GermanM1['t'] = pd.date_range('1961-3', '1996-1', freq='Q')
    reg = sm.ols(formula="dm~dy2+dR+dR1+dp+ecm_res+season",
                 data=historyM1).fit()
    pred = reg.predict(GermanM1)
    sd = np.std(reg.resid)
    return((GermanM1['dm'] - pred)/sd)


def plotResult(rGraph, df):
    GermanM1 = df['GermanM1']
    graph = rGraph.drop(rGraph.index[0])
    mon = GermanM1.loc[100:, :]
    mon['pi'] = graph.pi
    mon['piStar'] = graph.piStar
    mon['resid'] = graph.y
    mon.set_index('t', inplace=True)
    fig, (ax1, ax2) = plt.subplots(2, figsize=(10, 5.15))
    plt.style.use('ggplot')
    ax1.plot(mon.index, mon.resid, 'b-',
             label="Forecast Errors", linewidth=2, color='black')
    ax1.legend(loc='best')
    ax1.set_ylabel('Forecast Error')
    ax1.axvline(dt.datetime(1990, 12, 31), linewidth=2, color="black")
    ax2.plot(mon.index, mon.pi, 'b-',
             label="Detecting Statistics:$\pi$", linewidth=2, color="black")
    ax2.plot(mon.index, mon.piStar, '--',
             label="Sequential Stopping Procedure:$\mathcal{T}^{*}$", linewidth=2, color='black')
    ax2.axvline(dt.datetime(1990, 12, 31), linewidth=2, color="black")
    ax2.set_ylim(0, 1.1)
    ax2.legend(loc='best')
    ax2.set_ylabel('Probability')
    ax2.text(dt.datetime(1990, 12, 31), 0.6, "Q4 1990")
    return(fig)


if __name__ == "__main__":
    df = pyreadr.read_r('./data/GermanM1.rda')
    pred_resid = getResidule(df)
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
    b = {'runStr': 99,
         'rho': 0.01,
         'pi0': [0.5, 0.5],
         'p21': 0.00,
         'c': 0.01,
         'pstmf': pstmf,
         'vec': pred_resid}
    r = Bsquid.Bsquid(b)
    rGraph = Bsquid.Bsquid_crtGraph(b)
    g = plotResult(rGraph, df)
    g.savefig('./document/graph/German_M1.jpg')
    print("main is done")


# %%
