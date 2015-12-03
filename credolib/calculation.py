__all__=['correlmatrix']

import numpy as np

def correlmatrix(data1d, std_multiplier=None, logarithmic=True):
    cm=np.zeros((len(data1d),len(data1d)),np.double)
    for i in range(len(data1d)):
        for j in range(i+1,len(data1d)):
            I1=data1d[i].Intensity
            E1=data1d[i].Error
            I2=data1d[j].Intensity
            E2=data1d[j].Error
            idx=(E1>0)&(E2>0)
            if logarithmic:
                idx&=(I1>0)&(I2>0)
                E1=E1[idx]/I1[idx]
                E2=E2[idx]/I2[idx]
                I1=np.log(I1[idx])
                I2=np.log(I2[idx])
            w=(E1**2+E2**2)
            cm[i,j]=cm[j,i]=((I1-I2)**2/w).sum()/(1/w).sum()

    rowavg=cm.sum(axis=0)/(len(data1d)-1)
    cm=cm+np.diagflat(rowavg)
    if std_multiplier is None:
        return cm

    return cm, np.abs(rowavg-rowavg.mean())>std_multiplier*rowavg.std(), rowavg
