#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug  3 23:04:50 2019

@author: elrond
"""
import numpy as np
import scipy.stats
import os
os.chdir("/mnt/MyDoc/Dropbox/Research/MntStrBrk/code")
import numpy as np
import pandas as pd
import scipy.stats
import statsmodels.formula.api as sm
import pyreadr
import matplotlib.pyplot as plt
pd.options.display.float_format = '{:.2f}'.format
durab = pd.read_csv('../data/durab.csv')
durab['t'] = pd.date_range('1947-3','2001-5',freq='M')
durab_hist = durab.loc[202:393,]
durab_mon = durab.loc[202:,]
reg = sm.ols(formula="y~lag",data=durab_hist).fit()
reg_mon = sm.ols(formula="y~lag",data=durab_mon).fit()
plt.plot(reg_mon.resid)
print(reg.summary())
pred = reg.predict(durab_mon)
sd = np.std(reg.resid) 
pred_resid = (durab_mon['y'] - pred)/sd
y = pred_resid.values
plt.plot(y)
pi0 = [0.5,0.5]
p = np.array([[1-0.01, 0],
              [0.01, 1-0]])
p21 = 0
sigma1 = 0.5
brkStr = 300
brkSize = 0.8
smpl = 900
mu1 = 0.8
y = np.concatenate((np.random.normal(0,1,brkStr),np.random.normal(0,1+brkSize,smpl-brkStr)))

def getLikelihood(pi0,p,p21,mu1,y):
    l0 = scipy.stats.norm(0,1).pdf(y)
    l1 = scipy.stats.norm(mu1,1).pdf(y)
    pi = []
    f = []
    for i in range(len(y)):
        if i == 0:
                f.append(np.dot(np.dot(p,pi0),[l0[i],l1[i]]))
                pi.append(np.dot(p,pi0)*[l0[i],l1[i]]/f[i])
        else:
                f.append(np.dot(np.dot(p,pi[i-1]),[l0[i],l1[i]]))
                pi.append(np.dot(p,pi[i-1])*[l0[i],l1[i]]/f[i])
    return([mu1,sum(np.log(f)),pi[len(pi)-1][1]])
    
lk = getLikelihood(pi0,p,p21,mu1,y)

def getLikelihood_sigma(pi0,p,p21,sigma1,y):
    l0 = scipy.stats.norm(0,1).pdf(y)
    l1 = scipy.stats.norm(0,sigma1).pdf(y)
    pi = []
    f = []
    for i in range(len(y)):
        if i == 0:
                f.append(np.dot(np.dot(p,pi0),[l0[i],l1[i]]))
                pi.append(np.dot(p,pi0)*[l0[i],l1[i]]/f[i])
        else:
                f.append(np.dot(np.dot(p,pi[i-1]),[l0[i],l1[i]]))
                pi.append(np.dot(p,pi[i-1])*[l0[i],l1[i]]/f[i])
    return([sigma1,sum(np.log(f)),pi[len(pi)-1][1]])
    
lk = getLikelihood_sigma(pi0,p,p21,sigma1,y)