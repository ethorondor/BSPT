# %%
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import pyplot
import matplotlib.pyplot as plt
import timeit
import Bsquid
from operator import mul
import scipy.stats
import pandas as pd
import numpy as np
import multiprocessing as mp
import os
os.chdir("/home/elrond/Dropbox/Research/MntStrBrk")
pd.options.display.float_format = '{:.2f}'.format


def getTransitionMatrix():
    pmDef = {'mu1_p_l': 0.6,
             'mu1_p_u': 2.0,
             'mu1_n_l': -2,
             'mu1_n_u': -0.6,
             'mu1_sp_size': 50,
             'rho': rho,
             'smpl_sp_l': -4,
             'smpl_sp_u': 4,
             'smpl_sp_size': 200
             }
    return(Bsquid.getPS_TMF(pmDef))


def crtBsquid(pstmf, rho):
    np.random.seed(3218)
    brkStr = 150
    brkSize = 0.8
    smpl = 900
    y = np.random.normal(0, 1, smpl)
    y[brkStr:] = y[brkStr:] + brkSize
    b = {'runStr': 99,
         'rho': rho,
         'pi0': [0.5, 0.5],
         'p21': 0.00,
         'c': 0.05,
         'pstmf': pstmf,
         'vec': y}
    return(b)


def crt3DPlot():
    upper = postDist[0:50]
    lower = postDist[50:100]
    mid = pd.DataFrame(data=np.linspace(-0.576, 0.6, 50,
                                        endpoint=False), columns=['mu1'])
    df = pd.DataFrame()
    df = upper.append(mid)
    df = df.append(lower)
    df = df.fillna(0)
    data = df[df.columns[6:100]]
    x = df.mu1.values
    fig = plt.figure(figsize=(10, 5.15))
    # plt.style.use('dark_background')
    plt.style.use('seaborn-paper')
    ax = fig.add_subplot(111, projection='3d')
    y1 = np.full(150, 99)
    z1 = data[99].values
    y2 = np.full(150, 110)
    z2 = data[110].values
    y3 = np.full(150, 120)
    z3 = data[120].values
    y4 = np.full(150, 130)
    z4 = data[130].values
    y5 = np.full(150, 140)
    z5 = data[140].values
    y6 = np.full(150, 150)
    z6 = data[150].values
    y7 = np.full(150, 160)
    z7 = data[160].values
    y8 = np.full(150, 170)
    z8 = data[170].values
    y9 = np.full(150, 180)
    z9 = data[180].values
    y10 = np.full(150, 190)
    z10 = data[190].values
    ax.plot(x, y1, z1, 'b-')
    ax.plot(x, y2, z2, 'b-')
    ax.plot(x, y3, z3, 'b-')
    ax.plot(x, y4, z4, 'b-')
    ax.plot(x, y5, z5, 'b-')
    ax.plot(x, y6, z6, 'b-')
    ax.plot(x, y7, z7, 'b-')
    ax.plot(x, y8, z8, 'b-')
    ax.plot(x, y9, z9, 'b-')
    ax.plot(x, y10, z10, 'b-')
    ax.set_ylabel('Number of Obs')
    ax.set_xlabel("$\\theta$")
    plt.show()
    fig.savefig('./document/graph/threeDplot.png')


def PlotValueFunction(v):
    v.set_index(0, inplace=True)
    v['pi'] = 1 - v.index.values
    fig, ax = plt.subplots(figsize=(10, 5.15))
    plt.style.use('grayscale')
    ax.plot(v.index, v['pi'],'--', label="$1-\pi$", color='firebrick')
    ax.plot(v.index, v[20], 'b-', color='blue',
            label='value function iteration')
    ax.plot(v.index, v[50], 'b-', color='blue')
    ax.plot(v.index, v[70], 'b-', color='blue')
    ax.plot(v.index, v[100], 'b-', color='blue')
    ax.plot(v.index, v[110], 'b-', color='blue')
    ax.plot(v.index, v[120], 'b-', color='blue')
    ax.plot(v.index, v[127], 'b-', color='blue')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.legend(loc='best')
    ax.axvline(x=0.8, linewidth=2, color='black')
    ax.text(0.8, 0.8, "$\pi^{*}:0.8$")
    ax.xaxis.set_label_text('$\pi$')
    ax.yaxis.set_label_text('$v(\pi)$')
    fig.savefig('./document/graph/value_function.jpg')


def plotPostDist(postDist):
    upper = postDist[0:50]
    lower = postDist[50:100]
    mid = pd.DataFrame(data=np.linspace(-0.576, 0.6, 50,
                                        endpoint=False), columns=['mu1'])
    df = pd.DataFrame()
    df = upper.append(mid)
    df = df.append(lower)
    df = df.fillna(0)
    postDist.set_index('mu1', inplace=True, append=False, drop=True)
    fig0, ax = plt.subplots(figsize=(10, 5.15))
    plt.style.use('ggplot')
    ax.plot(df.mu1, df[99], label='100', linewidth=2)
    ax.plot(df.mu1, df[129], label='130', linewidth=2)
    ax.plot(df.mu1, df[159], label='160', linewidth=2)
    ax.plot(df.mu1, df[169], label='170', linewidth=2)
    ax.plot(df.mu1, df[189], label='190', linewidth=2)
    ax.legend(loc='best')
    ax.xaxis.set_label_text('Parameter Space')
    fig0.savefig('./document/graph/distribution.jpg')


def plotSimulation(r_graph):
    r_graph1 = r_graph.loc[100:200, :]
    fig, (ax1, ax2) = plt.subplots(2, figsize=(10, 5.15))
    plt.style.use('ggplot')
    ax1.plot(r_graph1.index, r_graph1.y, 'b-',
             label="Simulation", linewidth=2, color='black')
    ax1.axvline(x=150, linewidth=2, color='black')
    ax1.text(140, -1, "structural break at 150")
    ax1.legend(loc='best')
    ax1.set_ylabel('Simulation')
    ax2.plot(r_graph1.index, r_graph1.pi, '--', color='black',
             label="Detecting Statistics:$\pi$", linewidth=2)
    ax2.plot(r_graph1.index, r_graph1.piStar, 'b-',
             label="Sequential Stopping Procedure:$\mathcal{T}^{*}$", linewidth=2, color='black')
    ax2.axvline(x=160, linewidth=2, color='black')
    ax2.text(150, 0.6, "structural break announced at 160")
    ax2.legend(loc='best')
    ax2.set_ylabel('Probability')
    ax2.set_ylim(0, 1.1)
    fig.savefig('./document/graph/simu_mean_15.jpg')


if __name__ == "__main__":
    rho = 0.01
    pstmf = getTransitionMatrix()
    b = crtBsquid(pstmf, rho)
    v = Bsquid.Bsquid_getVF(b)
    r = Bsquid.Bsquid(b)
    r_graph = Bsquid.Bsquid_crtGraph(b)
    postDist = Bsquid.Bsquid_getPostDist(b)
    #crt3DPlot()
    PlotValueFunction(v)
    #plotPostDist(postDist)
    #plotSimulation(r_graph)
    print("this is the main function")


  # %%
