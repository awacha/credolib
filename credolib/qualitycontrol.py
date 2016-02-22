__all__=['assess_flux_stability','sum_measurement_times','assess_sample_stability','assess_instrumental_background','assess_transmission']
from IPython.core.getipython import get_ipython
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
import os
import ipy_table
from IPython.display import display
from sastool.misc.easylsq import nonlinear_leastsquares
from sastool.classes import SASHeader, SASMask
from sastool.misc.errorvalue import ErrorValue

def assess_flux_stability(samplename='Glassy_Carbon'):
    ip = get_ipython()
    f=plt.figure()
    ax1=f.add_subplot(1,1,1)
    plt.xlabel('Date of exposure')
    plt.ylabel('Beam flux (photon/sec), continuous lines')
    ax2=plt.twinx()
    plt.ylabel('Vacuum pressure (mbar), dotted lines')
    plt.title('Beam flux stability')
    lines=[]
    samplenames=sorted([sn_ for sn_ in ip.user_ns['_headers_tosave'] if samplename in sn_])
    linestyles=['b','g','r','c','m','y','k']
    lines=[]
    for sn,ls in zip(samplenames,linestyles):
        print(sn)
        heds=ip.user_ns['_headers_tosave'][sn]
        allheds=[]
        for k in heds.keys():
            allheds.extend(heds[k])
        allheds=sorted(allheds,key=lambda x:x['FSN'])
        fsns=np.array([h['FSN'] for h in allheds])
        cctstyleheaders=['devices.pilatus.exptime' in h for h in allheds]
        divisor=[[h['XPixel']*h['YPixel'],1]['devices.pilatus.exptime' in h] for h in allheds]
        factors=np.array([1/h['NormFactor']/d/0.96 for h,d in zip(allheds,divisor)])

        dates=[h['Date'] for h in allheds]
        lines.extend(ax1.plot(dates,factors,ls+'o',label='Flux (%s)'%sn))
        vacuums=np.array([h['Vacuum'] for h in allheds])
        lines.extend(ax2.plot(dates,vacuums,ls+'s',label='Vacuum (%s)'%sn,lw=2))
        print('  Measurement duration: %.2f h'%((dates[-1]-dates[0]).total_seconds()/3600.))
        print('  Mean flux: ',factors.mean(),'+/-',factors.std(),'photons/sec')
        print('  RMS variation of flux: ',factors.std()/factors.mean()*100,'%')
        print('  P-P variation of flux: ',(factors.max()-factors.min())/factors.mean()*100,'%')
    ax1.legend(lines,[l.get_label() for l in lines],loc='best')
    plt.show()

def sum_measurement_times():
    ip=get_ipython()
    _headers_tosave=ip.user_ns['_headers_tosave']
    tab=[['Sample name','Distance (mm)','Measurement time (h)']]
    for samplename in sorted(_headers_tosave):
        for dist in _headers_tosave[samplename]:
            tab.append([samplename, "%.f"%dist, sum([h['MeasTime'] for h in _headers_tosave[samplename][dist]])/3600])
    tab=ipy_table.IpyTable(tab)
    tab.apply_theme('basic')
    display(tab)

def assess_sample_stability(end_cutoff=3):
    ip=get_ipython()
    rowavg=ip.user_ns['_rowavg']
    tab=[['Sample name','Distance','Slope of autocorrelation function','Stability']]
    plt.figure()
    for sn in sorted(rowavg):
        for dist in sorted(rowavg[sn]):
            rowavg_rescaled=rowavg[sn][dist]/rowavg[sn][dist].mean()
            rowavg_std=rowavg_rescaled.std()*np.ones_like(rowavg_rescaled)
            try:
                A,B,stat=nonlinear_leastsquares(np.arange(len(rowavg_rescaled)),rowavg_rescaled,rowavg_std,lambda x,a,b:a*x+b,[0,0])
                problematic=(A.val>A.err*3)
            except TypeError:
                A='N/A'
                problematic=2
            tab.append([sn,'%.2f'%dist,A,["\u2713","\u2718\u2718\u2718\u2718\u2718",'\u274e'][problematic]])
            plt.errorbar(np.arange(len(rowavg_rescaled)),rowavg_rescaled,label=sn+' %.2f mm'%dist)#,diags_std[1:])
    plt.xlabel('Separation in time (FSN units)')
    plt.ylabel('Average discrepancy between curves')
    plt.legend(loc='best')
    tab=ipy_table.IpyTable(tab)
    tab.apply_theme('basic')
    display(tab)
    plt.show()

def assess_instrumental_background(Wx=20, Wy=20, emptyname='Empty_Beam', maskname='mask.mat', raw=True):
    ip = get_ipython()
    for l in os.listdir(ip.user_ns['saveto_dir']):
        if not (l.startswith(emptyname) and l.endswith('.npz') and
                ((not raw) or ('raw' in l))):
            continue
        print(l)
        I=np.load(os.path.join(ip.user_ns['saveto_dir'],l))['Intensity']
        header=SASHeader(os.path.join(ip.user_ns['saveto_dir'],l[:-4]+'.log'),plugin='CREDO Reduced')
        I=I/header['MeasTime']
        mask=SASMask(maskname,dirs=ip.user_ns['datadirs']).mask
        m=np.zeros_like(I)
        print('   Mean intensity per pixel:',I[mask==1].mean(),'cps')
        print('   STD intensity per pixel:',I[mask==1].std(),'cps')
        print('   Total intensity:',I[mask==1].sum(),'cps')
        for row in range(m.shape[0]):
            for col in range(m.shape[1]):
                m[row,col]=I[max(row-Wy,0):min(row+Wy,m.shape[0]-1),
                             max(col-Wx,0):min(col+Wx,m.shape[1]-1)][
                                 mask[max(row-Wy,0):min(row+Wy,m.shape[0]-1),
                                      max(col-Wx,0):
                                      min(col+Wx,m.shape[1]-1)]==1].mean()
        plt.figure()
        plt.subplot(1,2,1)
        plt.imshow(I,norm=LogNorm())
        plt.subplot(1,2,2)
        plt.imshow(m)
        plt.tight_layout()
        plt.suptitle(l)

def assess_transmission():
    ip = get_ipython()
    tab=[['Sample name','Distance','Transmission','Linear absorption coefficient (1/cm)','Absorption length (cm)']]
    for sample in sorted(ip.user_ns['_headers_tosave']):
        for dist in sorted(ip.user_ns['_headers_tosave'][sample]):
            transms_seen=[]
            for h in ip.user_ns['_headers_tosave'][sample][dist]:
                if h['Transm'] not in transms_seen:
                    transms_seen.append(h['Transm'])
                    transm=ErrorValue(h['Transm'],h['TransmError'])
                    thickness=ErrorValue(h['Thickness'],h['ThicknessError'])
                    mu=(-transm.log())/thickness
                    try:
                        invmu=(1/mu).tostring(extra_digits=2)
                    except ZeroDivisionError:
                        invmu='Infinite'
                    tab.append([sample,dist,transm.tostring(extra_digits=2),mu.tostring(extra_digits=2),invmu])
    tab=ipy_table.IpyTable(tab)
    tab.apply_theme('basic')
    display(tab)
