__all__=['guinieranalysis']

from .io import getsascurve
from .atsas import autorg, datgnom, shanum
from .utils import writemarkdown
import matplotlib.pyplot as plt
import numpy as np
from IPython.display import display
import ipy_table

def guinieranalysis(samplenames, qranges=None):
    figpr=plt.figure()
    axpr=figpr.add_subplot(1,1,1)
    if qranges is None:
        qranges={}
    table_autorg=[['Name', 'Rg (nm)', 'I$_0$ (cm$^{-1}$ sr$^{-1}$)', 'q$_{min}$ (nm$^{-1}$)', 'q$_{max}$ (nm$^{-1}$)', 'quality', 'aggregation', 'Dmax (nm)', 'q$_{shanum}$ (nm$^{-1}$)']]
    table_gnom=[['Name', 'Rg (nm)', 'I$_0$ (cm$^{-1}$ sr$^{-1}$)', 'qmin (nm$^{-1}$)', 'qmax (nm$^{-1}$)', 'Dmin (nm)', 'Dmax (nm)']]
    for sn in samplenames:
        if sn not in qranges:
            qrange=(0,1000000)
        else:
            qrange=qranges[sn]
        curve=getsascurve(sn)[0].trim(*qrange)
        curve.save(sn+'.dat')
        Rg, I0, qmin, qmax, quality, aggregation=autorg(sn+'.dat')
        dmax, nsh, nopt, qmaxopt=shanum(sn+'.dat')
        curve.trim(qmin,qmaxopt).save(sn+'_optrange.dat')
        gnompr,metadata=datgnom(sn+'_optrange.dat', Rg=Rg.val,noprint=True)
        axpr.errorbar(gnompr[:,0],gnompr[:,1],gnompr[:,2],None,label=sn)
        figsample=plt.figure()
        axgnomfit=figsample.add_subplot(1,2,1)
        curve.errorbar('b.',axes=axgnomfit,label='measured')
        axgnomfit.errorbar(metadata['qj'],metadata['jexp'], metadata['jerror'],None,'g.',label='gnom input')
        axgnomfit.loglog(metadata['qj'], metadata['jreg'], 'r-',label='regularized by GNOM')
        figsample.suptitle(sn)
        axgnomfit.set_xlabel('q (nm$^{-1}$)')
        axgnomfit.set_ylabel('$d\Sigma/d\Omega$ (cm$^{-1}$ sr$^{-1}$)')
        axgnomfit.axvline(qmaxopt,0,1,linestyle='dashed', color='black',lw=2)
        axgnomfit.grid(True, which='both')
        axgnomfit.axis('tight')
        axgnomfit.legend(loc='best')
        axguinier=figsample.add_subplot(1,2,2)
        axguinier.errorbar(curve.q,curve.Intensity, curve.Error, curve.qError,'.',label='Measured')
        q=np.linspace(qmin,qmax,100)
        axguinier.plot(q,I0.val*np.exp(-q**2*Rg.val**2/3),label='AutoRg')
        axguinier.plot(q,metadata['I0_gnom'].val*np.exp(-q**2*metadata['Rg_gnom'].val**2/3),label='Gnom')
        idxmin=np.arange(len(curve))[curve.q<=qmin].max()
        idxmax=np.arange(len(curve))[curve.q>=qmax].min()
        idxmin=max(0,idxmin-5)
        idxmax=min(len(curve)-1,idxmax+5)
        axguinier.set_xscale('power',exponent=2)
        axguinier.set_yscale('log')
        axguinier.set_xlabel('q (nm$^{-1}$)')
        axguinier.set_ylabel('$d\Sigma/d\Omega$ (cm$^{-1}$ sr$^{-1}$)')
        axguinier.legend(loc='best')
        curveguinier=curve.trim(curve.q[idxmin],curve.q[idxmax])
        axguinier.axis(xmax=curve.q[idxmax],xmin=curve.q[idxmin], ymin=curveguinier.Intensity.min(),ymax=curveguinier.Intensity.max())
        axguinier.grid(True, which='both')
        table_gnom.append([sn, metadata['Rg_gnom'].tostring(extra_digits=2), metadata['I0_gnom'].tostring(extra_digits=2), metadata['qmin'], metadata['qmax'], metadata['dmin'], metadata['dmax']])
        table_autorg.append([sn,Rg.tostring(extra_digits=2), I0,'%.3f'%qmin,'%.3f'%qmax,'%.1f %%'%(quality*100), aggregation, '%.3f'%dmax, '%.3f'%qmaxopt])
    axpr.set_xlabel('r (nm)')
    axpr.set_ylabel('P(r)')
    axpr.legend(loc='best')
    axpr.grid(True,which='both')
    writemarkdown('## Results from autorg and shanum')
    tab=ipy_table.IpyTable(table_autorg)
    tab.apply_theme('basic')
    display(tab)
    writemarkdown('## Results from gnom')
    tab=ipy_table.IpyTable(table_gnom)
    tab.apply_theme('basic')
    display(tab)

