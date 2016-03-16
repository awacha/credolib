__all__=['guinieranalysis']

from .io import getsascurve
from .atsas import autorg, datgnom, shanum
from .utils import writemarkdown
import matplotlib.pyplot as plt
import numpy as np
from IPython.display import display
from IPython.core.getipython import get_ipython
import ipy_table
import os

def guinieranalysis(samplenames, qranges=None,qmax_from_shanum=True, prfunctions_postfix='', dist=None, plotguinier=True, graph_extension='.png'):
    figpr=plt.figure()
    ip=get_ipython()
    axpr=figpr.add_subplot(1,1,1)
    if qranges is None:
        qranges={}
    table_autorg=[['Name', 'Rg (nm)', 'I$_0$ (cm$^{-1}$ sr$^{-1}$)', 'q$_{min}$ (nm$^{-1}$)', 'q$_{max}$ (nm$^{-1}$)', 'quality', 'aggregation', 'Dmax (nm)', 'q$_{shanum}$ (nm$^{-1}$)']]
    table_gnom=[['Name', 'Rg (nm)', 'I$_0$ (cm$^{-1}$ sr$^{-1}$)', 'qmin (nm$^{-1}$)', 'qmax (nm$^{-1}$)', 'Dmin (nm)', 'Dmax (nm)']]
    results={}
    for sn in samplenames:
        if sn not in qranges:
            qrange=(0,1000000)
        else:
            qrange=qranges[sn]
        curve=getsascurve(sn,dist)[0].trim(*qrange)
        curve.save(sn+'.dat')
        try:
            Rg, I0, qmin, qmax, quality, aggregation=autorg(sn+'.dat')
        except ValueError:
            print('Error running autorg on %s'%sn)
            continue
        dmax, nsh, nopt, qmaxopt=shanum(sn+'.dat')
        if qmax_from_shanum:
            curve.trim(qmin,qmaxopt).save(sn+'_optrange.dat')
        else:
            curve.trim(qmin,qrange[1]).save(sn+'_optrange.dat')
        gnompr,metadata=datgnom(sn+'_optrange.dat', Rg=Rg.val,noprint=True)
        axpr.errorbar(gnompr[:,0],gnompr[:,1],gnompr[:,2],None,label=sn)
        if plotguinier:
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
            axguinier.set_xscale('power',exponent=2)
            axguinier.set_yscale('log')
            axguinier.set_xlabel('q (nm$^{-1}$)')
            axguinier.set_ylabel('$d\Sigma/d\Omega$ (cm$^{-1}$ sr$^{-1}$)')
            axguinier.legend(loc='best')
        idxmin=np.arange(len(curve))[curve.q<=qmin].max()
        idxmax=np.arange(len(curve))[curve.q>=qmax].min()
        idxmin=max(0,idxmin-5)
        idxmax=min(len(curve)-1,idxmax+5)
        if plotguinier:
            curveguinier=curve.trim(curve.q[idxmin],curve.q[idxmax])
            axguinier.axis(xmax=curve.q[idxmax],xmin=curve.q[idxmin], ymin=curveguinier.Intensity.min(),ymax=curveguinier.Intensity.max())
            axguinier.grid(True, which='both')
        table_gnom.append([sn, metadata['Rg_gnom'].tostring(extra_digits=2), metadata['I0_gnom'].tostring(extra_digits=2), metadata['qmin'], metadata['qmax'], metadata['dmin'], metadata['dmax']])
        table_autorg.append([sn,Rg.tostring(extra_digits=2), I0,'%.3f'%qmin,'%.3f'%qmax,'%.1f %%'%(quality*100), aggregation, '%.3f'%dmax, '%.3f'%qmaxopt])
        if plotguinier:
            figsample.tight_layout()
            figsample.savefig(os.path.join(ip.user_ns['auximages_dir'],'guinier_%s%s'%(sn,graph_extension)),dpi=600)
        results[sn]={'Rg_autorg':Rg,'I0_autorg':I0,
                     'qmin_autorg':qmin,'qmax_autorg':qmax,
                     'quality':quality,'aggregation':aggregation,
                     'dmax_autorg':dmax,'qmax_shanum':qmaxopt,
                     'Rg_gnom':metadata['Rg_gnom'],
                     'I0_gnom':metadata['I0_gnom'],
                     'qmin_gnom':metadata['qmin'],
                     'qmax_gnom':metadata['qmax'],
                     'dmin_gnom':metadata['dmin'],
                     'dmax_gnom':metadata['dmax']}
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
    if prfunctions_postfix and prfunctions_postfix[0]!='_':
        prfunctions_postfix='_'+prfunctions_postfix
    figpr.tight_layout()
    figpr.savefig(os.path.join(ip.user_ns['auximages_dir'],'prfunctions%s%s'%(prfunctions_postfix,graph_extension)),dpi=600)
    display(tab)
    return results
