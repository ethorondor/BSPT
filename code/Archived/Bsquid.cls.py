import numpy as np
import pandas as pd
import scipy.stats
from operator import mul
import multiprocessing as mp
pool = mp.Pool(processes=16)
np.set_printoptions(suppress=True)

class Bsquid(object):
     # initiate a Bsquid class
     # The param is a dictionary contains all the information regarding the Bsquid
     # param = {'mu1_p_l':0.5,
     #          'mu1_p_u':1.5,
     #          'mu1_n_l':-1.5,
     #          'mu1_n_u':-0.5,
     #          'mu1_sp_size':100,
     #          'rho':0.01,
     #          'pi0':[0.5,0.5]
     #          'p21':0,
     #          'c':0.01}
     # bs = Bsquid(param, y)
    def __init__(self, param, vec):
        self.param = param
        self.vec = vec


    def main(self):
        pmSp = getPmSp(self)
        T = getTbar(self,pmSp)
        piStar = getPistar(self,T)
        pool.close()
        pool.join()
        mu1_bar = sum(pmSp['w']*pmSp['mu1'])
        pi_bar = sum(pmSp['w']*pmSp['pi']) 
        return(mu1_bar,pi_bar,piStar)   

    def getPmSp(self):
        rho = self.param['rho']
        pi0 = self.param['pi0']
        p21 = self.param['p21']
        y = self.vec
        mu1Sp_n = list(np.linspace(self.param['mu1_n_l'],
                        self.param['mu1_n_u'],self.param['mu1_sp_size']))
        mu1Sp_p = list(np.linspace(self.param['mu1_p_l'],
                        self.param['mu1_p_u'],self.param['mu1_sp_size']))
        lk_n = pd.DataFrame([pool.apply(getLikelihood,args=(self,mu1)) for mu1 in mu1Sp_n],
                columns=['mu1','l','pi'])
        lk_p = pd.DataFrame([pool.apply(getLikelihood,args=(self,mu1)) for mu1 in mu1Sp_p],
                columns=['mu1','l','pi'])
        if sum(lk_n['l']) > sum(lk_p['l']):
            pmSp = lk_n
        else:
            pmSp = lk_p
        pmSp['w'] = np.exp(pmSp['l'] - pmSp['l'].mean(axis=0))/sum(np.exp(pmSp['l'] - pmSp['l'].mean(axis=0)))
        return(pmSp)

    def getLikelihood(self,mu1):
        rho = self.param['rho']
        pi0 = self.param['pi0']
        p21 = self.param['p21']
        y = self.vec
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

    def getTbar(self,pmSp):
        rho = self.param['rho']
        TrMx = [pool.apply(getTM,args=(self,mu1)) for mu1 in pmSp['mu1']]
        TT = list(map(mul, TrMx,pmSp['w']))
        return(np.array(sum(TT)))

    def getTM(self,mu1):
        rho = self.param['rho']
        smplSp = list(np.linspace(-3,3,300))
        piSp = np.array(range(0,101))/100
        df = pd.DataFrame(np.zeros((len(smplSp),6)),columns=['l0','l1','p0','p1','lr','pi_prime'])
        df['l0'] = scipy.stats.norm(0,1).pdf(smplSp)
        df['l1'] = scipy.stats.norm(mu1,1).pdf(smplSp)
        df['p0'] = df['l0']/sum(df['l0'])
        df['p1'] = df['l1']/sum(df['l1'])
        df['lr'] = df['l1']/df['l0']
        TM = np.empty([len(piSp),len(piSp)])
        for i in range(len(piSp)): 
            TM[i,] = getTV(self,piSp[i],df)
        return(TM)

    def getTV(self,pi,df):
        rho = self.param['rho']
        _T = pd.DataFrame(np.arange(0,101)/100,columns=['pi_prime'])
        df['pi_prime'] = np.round((df['lr']*(pi+rho*(1-pi)))/(df['lr']*(pi+rho*(1-pi))+(1-rho)*(1-pi)),2)
        d = df.groupby('pi_prime',as_index=False).sum()
        d['p'] = (1-pi)*(1-rho)*d['p0']+(pi+(1-pi)*rho)*d['p1']
        _T = pd.merge(_T,d[['pi_prime','p']],on='pi_prime',how='left')
        _T = np.array(_T.iloc[:,1].fillna(0))
        return(_T.reshape(1,-1))

    def getPistar(self,T):
        c = self.param['c']
        piSp = np.array(range(0,101))/100
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
     



 



