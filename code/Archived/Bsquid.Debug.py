import numpy as np
import pandas as pd
import scipy.stats
from operator import mul
np.set_printoptions(suppress=True)
######################################################
######### define function ############################
######################################################
def Bsquid(theObject):
    pmSp = getPmSp(theObject)
    T = getTbar(pmSp,theObject)
    piStar = getPistar(theObject,T)
    mu1_bar = sum(pmSp['w']*pmSp['mu1'])
    pi_bar = sum(pmSp['w']*pmSp['pi']) 
    return(mu1_bar,pi_bar,piStar)   

def getPmSp(theObject):
    #rho = theObject['rho']
    #pi0 = theObject['pi0']
    #p21 = theObject['p21']
    #y = theObject['vec']
    mu1Sp_n = list(np.linspace(theObject['mu1_n_l'],
                    theObject['mu1_n_u'],theObject['mu1_sp_size']))
    mu1Sp_p = list(np.linspace(theObject['mu1_p_l'],
                    theObject['mu1_p_u'],theObject['mu1_sp_size']))

    lk_n = pd.DataFrame([getLikelihood(mu1,theObject=theObject) for mu1 in mu1Sp_n],
            columns=['mu1','l','pi'])
    lk_p = pd.DataFrame([getLikelihood(mu1,theObject=theObject) for mu1 in mu1Sp_p],
            columns=['mu1','l','pi'])
    if sum(lk_n['l']) > sum(lk_p['l']):
        pmSp = lk_n
    else:
        pmSp = lk_p
    pmSp['w'] = np.exp(pmSp['l'] - pmSp['l'].mean(axis=0))/sum(np.exp(pmSp['l'] - pmSp['l'].mean(axis=0)))
    return(pmSp)

def getLikelihood(mu1,theObject):
    rho = theObject['rho']
    pi0 = theObject['pi0']
    p21 = theObject['p21']
    y = theObject['vec']
    p = [[1-rho, p21],
        [rho, 1-p21]]
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

def getTbar(pmSp,theObject):
    #rho = theObject['rho']
    TrMx = [getTM(mu1,theObject=theObject) for mu1 in pmSp['mu1']]
    TT = list(map(mul, TrMx,pmSp['w']))
    return(np.array(sum(TT)))

def getTM(mu1,theObject):
    #rho = theObject['rho']
    smplSp = list(np.linspace(-2,2,200))
    #piSp = np.array(range(0,101))/100
    piSp = np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
    df = pd.DataFrame(np.zeros((len(smplSp),6)),columns=['l0','l1','p0','p1','lr','pi_prime'])
    df['l0'] = scipy.stats.norm(0,1).pdf(smplSp)
    df['l1'] = scipy.stats.norm(mu1,1).pdf(smplSp)
    df['p0'] = df['l0']/sum(df['l0'])
    df['p1'] = df['l1']/sum(df['l1'])
    df['lr'] = df['l1']/df['l0']
    TM = np.empty([len(piSp),len(piSp)])
    for i in range(len(piSp)): 
        TM[i,] = getTV(theObject,piSp[i],df)
    return(TM)

def getTV(theObject,pi,df):
    rho = theObject['rho']
    #_T = pd.DataFrame(np.arange(0,101)/100,columns=['pi_prime'])
    _T = pd.DataFrame([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0],columns=['pi_prime'])
    
    df['pi_prime'] = np.round((df['lr']*(pi+rho*(1-pi)))/(df['lr']*(pi+rho*(1-pi))+(1-rho)*(1-pi)),1)
    #df['pi_prime'] = np.round(df.pi_prime[df.pi_prime<0.6],1)
    d = df.groupby('pi_prime',as_index=False).sum()
    d['p'] = (1-pi)*(1-rho)*d['p0']+(pi+(1-pi)*rho)*d['p1']
    _T = pd.merge(_T,d[['pi_prime','p']],on='pi_prime',how='left')
    _T = np.array(_T.iloc[:,1].fillna(0))
    return(_T.reshape(1,-1))

def getPistar(theObject,T):
    c = theObject['c']
    #piSp = np.array(range(0,101))/100
    piSp = np.array([0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0])
    h = c*np.matmul(T,piSp)+1+(c-1)*piSp
    g = 1-piSp
    Q = np.amin([g,h],axis=0)
    ep = 1
    cnt = 0
    while ep > 0.001 and cnt < 500:
        cnt + 1
        Q1 = Q
        h = np.matmul(T,Q) + c*piSp
        Q = np.amin([g,h],axis=0)
        ep = max(abs(Q1-Q))
        diff = abs(g-h)
        if cnt < 500:
            piStar = piSp[np.argmin(diff)]
        else:
            piStar = 'NaN'
    return(piStar)
import os
os.chdir("/mnt/MyDoc/Dropbox/Research/MntStrBrk/code")
import numpy as np

brkStr = 109
brkSize = 0.8
smpl = 600
y = np.random.normal(0,1,smpl)
y[brkStr:] = y[brkStr:] + brkSize
param = {'mu1_p_l':0.5,
          'mu1_p_u':1.5,
          'mu1_n_l':-1.5,
          'mu1_n_u':-0.5,
          'mu1_sp_size':100,
          'rho':0.01,
          'pi0':[0.5,0.5],
          'p21':0,
          'c':0.01,
          'vec':y}
b = Bsquid(param)

"""
def getSimu(c):
    np.random.seed(c)
    brkStr = 109
    brkSize = 0.8
    smpl = 600
    y = np.random.normal(0,1,smpl)
    y[brkStr:] = y[brkStr:] + brkSize
    for i in range(98,len(y)):
        param = {'mu1_p_l':0.5,
                'mu1_p_u':1.5,
                'mu1_n_l':-1.5,
                'mu1_n_u':-0.5,
                'mu1_sp_size':100,
                'rho':0.01,
                'pi0':[0.5,0.5],
                'p21':0,
                'c':0.01,
                'vec':y[0:i]}
        b = Bsquid(param)
        if b[1] > b[2]:
            return([b[0],b[1],b[2],i])
            break
    return(bs)

    
##############################################################
#################  start simulation  #########################
##############################################################
getSimu(1)

nSim = 10
pool = mp.Pool()
result = pool.map(getSimu,np.arange(nSim))
#result = [pool.apply(getSimu,args=(simPam,c)) for c in list(np.arange(nSim))]
pool.close()
pool.join()

df = pd.DataFrame(result,columns=['mu1','pi','piStar','os'])
df.to_csv("./results.csv", sep=',',index=False)
"""