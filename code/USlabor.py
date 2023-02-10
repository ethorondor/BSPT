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
os.chdir("/home/elrond/Documents/Dropbox/Research/MntStrBrk")
pd.options.display.float_format = '{:.2f}'.format


def getResidual(durab):
    durab['t'] = pd.date_range('1947-3', '2001-5', freq='M')
    durab.set_index('t', inplace=True)
    durab_hist = durab.loc['1964-01-31':'1980-12-31', :]
    durab_mon = durab.loc['1980-01-31':'1997-12-31', ]
    reg = sm.ols(formula="y~lag", data=durab_hist).fit()
    # print(reg.summary())
    pred = reg.predict(durab_mon)
    sd = np.std(reg.resid)
    durab_mon['resid'] = (durab_mon['y'] - pred)/sd
    return(durab_mon)


def plotResidual(rgraph, durab_mon):
    graph = rgraph.drop(rgraph.index[0])
    mon = durab_mon.loc['1980-06-30':, :]
    mon['pi'] = graph.pi.values
    mon['piStar'] = graph.piStar.values
    mon['resid_d'] = graph.y.values
    fig, (ax1, ax2) = plt.subplots(2, figsize=(10, 5.15))
    plt.style.use('ggplot')
    ax1.plot(mon.index, mon.resid, 'b-', label="Forecast Errors",
             linewidth=2, color='black')
    ax1.legend(loc='best')
    ax1.set_ylabel('Forecast Error')
    ax1.axvline(dt.datetime(1997, 8, 31), linewidth=2, color='black')
    ax2.text(dt.datetime(1995, 9, 30), 0.4, "Aug 1997")
    ax1.hlines(y=0, xmin=dt.datetime(1980, 6, 30),
               xmax=dt.datetime(1997, 12, 31), colors='black')
    ax2.plot(mon.index, mon.pi, '--', label="Detecting Statistics:$\pi$",
             linewidth=2, color='black')
    ax2.plot(mon.index, mon.piStar, 'b-',
             label="Sequential Stopping Procedure:$\mathcal{T}^{*}$", linewidth=2, color='black')
    ax2.axvline(dt.datetime(1997, 8, 31), linewidth=2, color='black')
    ax2.legend(loc='best')
    ax2.set_ylabel('Probability')
    ax2.text(dt.datetime(1995, 9, 30), 0.4, "Aug 1997")
    return(fig)


if __name__ == "__main__":
    durab = pd.read_csv('./data/durab.csv')
    durab_mon = getResidual(durab)
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
    b = {'runStr': 4,
         'rho': 0.01,
         'pi0': [0.5, 0.5],
         'p21': 0.00,
         'c': 0.01,
         'pstmf': pstmf,
         'vec': durab_mon.resid.values}
    r = Bsquid.Bsquid(b)
    rgraph = Bsquid.Bsquid_crtGraph(b)
    fig = plotResidual(rgraph, durab_mon)
    fig.savefig('./document/graph/US_labor_prod.jpg')
    print("the main is done")
