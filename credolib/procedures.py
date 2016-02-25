__all__= ['summarize','unite','subtract_bg']

from IPython.core.getipython import get_ipython
from IPython.display import display
from .utils import writemarkdown, print_abscissavalue, putlogo
from .calculation import correlmatrix
from .plotting import plotsascurve
from .io import get_different_distances
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid import make_axes_locatable
import numpy as np
import os
from sastool.classes import SASExposure, SASCurve
from sastool.misc.errorvalue import ErrorValue
from sastool.misc.easylsq import nonlinear_odr, FixedParameter
from sastool.libconfig import qunit
import ipy_table


def _collect_data_for_summarization(headers,raw,reintegrate,qrange):
    ip=get_ipython()
    data1d=[]
    data2d=0
    headersout=[]
    if not headers:
        return
    for head in headers:
        mo = ip.user_ns['mask_override'](head)
        if raw:
            fileformat = ip.user_ns['crd_prefix']+'_%05d.cbf'
            dirs = ip.user_ns['datadirs']
        else:
            fileformat = ip.user_ns['crd_prefix']+'_%05d.npz'
            dirs = ip.user_ns['evaldirs']
        try:
            if mo is not None:
                ex = SASExposure(
                    fileformat, head['FSN'], dirs=dirs, maskfile=mo)
            else:
                ex = SASExposure(
                    fileformat, head['FSN'], dirs=dirs)
        except IOError:
            print('Could not load 2D file: %s'%(fileformat%head['FSN']))
            ip.user_ns['badfsns'].append(head['FSN'])
            continue
        ex.header = head
#                 for keytodelete in ['Wavelength', 'WavelengthCalibrated']:
#                     try:
#                         del ex.header[keytodelete]
#                     except KeyError:
#                         pass
        if not reintegrate:
            data1d.append(SASCurve(
                os.path.join(ip.user_ns['onedim_folder'], (ip.user_ns['crd_prefix']+'_%d.txt') % head['FSN'])))
        else:
            data1d.append(
                ex.radial_average(qrange, errorpropagation=3,
                                  abscissa_errorpropagation=3))
        data1d[-1].save(os.path.join(ip.user_ns['saveto_dir'],'curve_%05d.txt'%head['FSN']))
        data2d=data2d+ex
        headersout.append(ex.header)
    data2d /= len(data1d)
    data2d['MeasTime'] = sum(
        [h['MeasTime'] for h in headersout])
    return data1d, data2d, headersout

def summarize(reintegrate=True, dist_tolerance=3, qranges=None,
              samples=None, raw=False, late_radavg=True, graph_ncols=3,
              std_multiplier=3, graph_extension='png',
              graph_dpi=80, correlmatrix_colormap='coolwarm',
              image_colormap='viridis'):
    """Summarize scattering patterns and curves for all samples defined 
    by the global `allsamplenames`.
    
    Inputs:
        reintegrate (bool, default=True): if the curves are to be obained
            by reintegrating the patterns. Otherwise 1D curves are loaded.
        dist_tolerance (float, default=3): sample-to-detector distances
            nearer than this are considered the same
        qranges (dict): a dictionary mapping approximate sample-to-detector
            distances (within dist_tolerance) to one-dimensional np.ndarrays
            of the desired q-range of the reintegration.
        samples (list or None): the names of the samples to summarize. If
            None, all samples defined by ``allsamplenames`` are used.
        raw (bool, default=False): if raw images are to be treated instead
            the evaluated ones (default).
        late_radavg (bool, default=True): if the scattering curves are to
            be calculated from the summarized scattering pattern. If False,
            scattering curves are calculated from each pattern and will be
            averaged.
        graph_ncols: the number of columns in graphs (2D patterns, 
            correlation matrices)
        std_multiplier: if the absolute value of the relative discrepancy 
            is larger than this limit, the exposure is deemed an outlier.
        graph_extension: the extension of the produced hardcopy files.
        graph_dpi: resolution of the graphs
        correlmatrix_colormap: name of the colormap to be used for the
            correlation matrices (resolved by matplotlib.cm.get_cmap())
        image_colormap: name of the colormap to be used for the scattering
            patterns (resolved by matplotlib.cm.get_cmap())
    """
    if qranges is None:
        qranges={}
    ip = get_ipython()
    data2d = {}
    data1d = {}
    headers_tosave = {}
    rowavg={}
    if raw:
        writemarkdown('# Summarizing RAW images.')
        headers = ip.user_ns['_rawheaders']
        rawpart = '_raw' # this will be added in the filenames saved
    else:
        writemarkdown('# Summarizing CORRECTED images.')
        headers = ip.user_ns['_evalheaders']
        rawpart = '' # nothing will be added in the filenames saved
    
    if samples is None:
        samples = sorted(ip.user_ns['allsamplenames'])
    for samplename in samples:
        writemarkdown('## '+samplename)
        headers_sample=[h for h in headers if h['Title']==samplename]
        data2d[samplename] = {}
        data1d[samplename] = {}
        rowavg[samplename] = {}
        headers_tosave[samplename] = {}
        dists=get_different_distances([h for h in headers if h['Title']==samplename], dist_tolerance)
        if not dists:
            writemarkdown('No measurements from sample, skipping.')
            continue
        fig_2d = plt.figure()
        fig_curves=plt.figure()
        fig_correlmatrices = plt.figure()
        distaxes = {}
        correlmatrixaxes={}
        ncols=min(len(dists),graph_ncols)
        nrows=int(np.ceil(len(dists)/ncols))
        onedimaxes=fig_curves.add_axes((0.1,0.3,0.8,0.5))
        onedimstdaxes=fig_curves.add_axes((0.1,0.1,0.8,0.2))
        for distidx, dist in enumerate(dists):
            writemarkdown("### Distance "+str(dist)+" mm")
            headers_narrowed=[h for h in headers_sample if abs(h['DistCalibrated']-dist)<dist_tolerance]
            distaxes[dist] = fig_2d.add_subplot(
                nrows, ncols, distidx + 1)
            correlmatrixaxes[dist]=fig_correlmatrices.add_subplot(
                nrows, ncols, distidx+1)
            # determine the q-range to be used from the qranges argument.
            try:
                distkey_min = min([np.abs(k - dist)
                                  for k in qranges if np.abs(k - dist) < dist_tolerance])
            except ValueError:
                # no matching key in qranges dict
                qrange = None  # request auto-determination of q-range
            else:
                distkey = [
                    k for k in qranges if np.abs(k - dist) == distkey_min][0]
                qrange = qranges[distkey]

            (data1d[samplename][dist], data2d[samplename][dist],headers_tosave[samplename][dist]) = \
                _collect_data_for_summarization(headers_narrowed,raw, reintegrate, qrange)

            # calculate and plot correlation matrix
            cmatrix,badidx, rowavg[samplename][dist]=correlmatrix(data1d[samplename][dist],std_multiplier)
            rowavgmean=rowavg[samplename][dist].mean()
            rowavgstd=rowavg[samplename][dist].std()
            writemarkdown('#### Assessing sample stability')
            writemarkdown("- Mean of row averages: "+str(rowavgmean))
            writemarkdown("- Std of row averages: "+str(rowavgstd)+' (%.2f %%)'% (rowavgstd/rowavgmean*100))

            img=correlmatrixaxes[dist].imshow(cmatrix,interpolation='nearest',cmap=matplotlib.cm.get_cmap(correlmatrix_colormap))
            cax=make_axes_locatable(correlmatrixaxes[dist]).append_axes('right',size="5%",pad=0.05)
            fig_correlmatrices.colorbar(img,cax=cax)
            correlmatrixaxes[dist].set_title('%.2f mm'%dist)
            fsnaxes=np.array([h['FSN'] for h in headers_tosave[samplename][dist]])
            correlmatrixaxes[dist].set_xticks(list(range(len(data1d[samplename][dist]))))
            correlmatrixaxes[dist].set_xticklabels([str(f) for f in fsnaxes],rotation='vertical')
            correlmatrixaxes[dist].set_yticks(list(range(len(data1d[samplename][dist]))))
            correlmatrixaxes[dist].set_yticklabels([str(f) for f in fsnaxes])
            np.savez_compressed(os.path.join(ip.user_ns['saveto_dir'],
                'correlmatrix_%s_%s' % (
                    samplename,
                    ('%.2f' % dist).replace('.' , '_')) + rawpart + '.npz'),
                correlmatrix=cmatrix, fsns=fsnaxes)

            # Plot the image
            try:
                data2d[samplename][dist].plot2d(
                   zscale='log10', axes=distaxes[dist], crosshair=False,
                   cmap=matplotlib.cm.get_cmap(image_colormap))
            except ValueError:
                print('Error plotting 2D image for sample %s, distance %.2f'%(samplename, dist))
            distaxes[dist].set_xlabel('q (' + qunit() + ')')
            distaxes[dist].set_ylabel('q (' + qunit() + ')')
            distaxes[dist].set_title(
                '%.2f mm (%d curve%s)' % (dist, len(headers_tosave[samplename][dist]),
                    ['', 's'][len(headers_tosave[samplename][dist]) > 1]))

            # Plot the curves
            Istd=np.zeros((len(data1d[samplename][dist][0]),len(data1d[samplename][dist])),np.float)
            for i,c in enumerate(data1d[samplename][dist]):
                c.loglog(axes=onedimaxes)
                Istd[:,i]=c.Intensity
            if Istd.shape[1]>1:
                onedimstdaxes.loglog(data1d[samplename][dist][0].q,Istd.std(axis=1)/Istd.mean(axis=1)*100,'b-')
            if not late_radavg:
                data1d[samplename][dist] = SASCurve.average(
                    *data1d[samplename][dist])
            else:
                data1d[samplename][dist] = (
                    data2d[samplename][dist].radial_average(
                        qrange,
                        errorpropagation=3,
                        abscissa_errorpropagation=3))
            data1d[samplename][dist].loglog(
                label='Average', lw=2, color='k', axes=onedimaxes)

            #Saving image, headers, mask and curve
            data2d[samplename][dist].write(
                os.path.join(ip.user_ns['saveto_dir'],
                             samplename + '_'+(
                                 '%.2f' % dist).replace('.', '_') +
                             rawpart + '.npz'), plugin='CREDO Reduced')
            data2d[samplename][dist].header.write(
                os.path.join(ip.user_ns['saveto_dir'],
                             samplename + '_'+(
                                 '%.2f' % dist).replace('.', '_') +
                             rawpart +'.log'), plugin='CREDO Reduced')
            data2d[samplename][dist].mask.write_to_mat(
                os.path.join(ip.user_ns['saveto_dir'],
                             data2d[samplename][dist].mask.maskid+'.mat'))
            data1d[samplename][dist].save(os.path.join(ip.user_ns['saveto_dir'], samplename + '_'+('%.2f' % dist).replace('.', '_') + rawpart + '.txt'))

            #Report table on sample stability
            tab=[['FSN','Date','Discrepancy','Relative discrepancy ((x-mean(x))/std(x))','Quality']]
            discrmean=rowavg[samplename][dist].mean()
            discrstd=rowavg[samplename][dist].std()
            for h,bad,discr in zip(headers_tosave[samplename][dist],
                                   badidx,rowavg[samplename][dist]):
                tab.append([h['FSN'],h['Date'].isoformat(),discr,(discr-discrmean)/discrstd,["\u2713","\u2718\u2718\u2718\u2718\u2718"][bad]])
                if bad:
                    ip.user_ns['badfsns'].append(h['FSN'])
            tab=ipy_table.make_table(tab)
            ipy_table.apply_theme('basic')
            display(tab)

            # Report on qrange and flux
            q_ = data1d[samplename][dist].q
            qmin = q_[q_ > 0].min()
            writemarkdown('#### Q-range & flux')
            writemarkdown('- $q_{min}$: '+ print_abscissavalue(qmin, data2d[samplename][dist]['Wavelength'], dist))
            writemarkdown('- $q_{max}$: '+ print_abscissavalue(data1d[samplename][dist].q.max(), data2d[samplename][dist]['Wavelength'], dist))
            writemarkdown('- Number of $q$ points: '+str(len(data1d[samplename][dist])))
            try:
                normfac = [h['NormFactor']
                           for h in headers_tosave[samplename][dist]]
                cct_type_headers=['devices.pilatus.exptime' in h for h in headers_tosave[samplename][dist]]
                if not any (cct_type_headers):
                    flux=(1 / ErrorValue(np.mean(normfac), np.std(normfac))/
                          headers_tosave[samplename][dist][0]['XPixel']/
                          headers_tosave[samplename][dist][0]['YPixel']/
                          0.96).tostring()
                elif all(cct_type_headers):
                    flux=(1/ErrorValue(np.mean(normfac),np.std(normfac))/0.96).tostring()
                else:
                    flux='Mixed SAXSCtrl and CCT headers, cannot estimate normfactor correctly.'
                writemarkdown("- beam flux (photon/sec): %s"%flux)
                writemarkdown("- from %d exposures, total exposure time %.0f sec <=> %.2f hr" % (len(headers_tosave[samplename][dist]), data2d[samplename][dist]['MeasTime'], data2d[samplename][dist]['MeasTime']/3600.))
            except KeyError:
                writemarkdown("- *No information on beam flux: dealing with raw data.*")
        onedimaxes.set_xlabel('')
        onedimaxes.set_ylabel('$d\\Sigma/d\\Omega$ (cm$^{-1}$ sr$^{-1}$)')
        # plt.legend(loc='best')
        onedimaxes.grid(True,which='both')
        onedimaxes.axis('tight')
        onedimaxes.set_title(samplename)
        onedimstdaxes.set_xlabel('q (' + qunit() + ')')
        onedimstdaxes.set_ylabel('Rel.std.dev. of intensity (%)')
        onedimstdaxes.grid(True, which='both')
        onedimstdaxes.set_xlim(*onedimaxes.get_xlim())
        putlogo(fig_curves)
        putlogo(fig_2d)
        fig_2d.tight_layout()
        fig_correlmatrices.suptitle(samplename)
        fig_correlmatrices.tight_layout()
        #fig_curves.tight_layout()
        fig_2d.savefig(
            os.path.join(ip.user_ns['auximages_dir'],
                         'averaging2D_' +
                         samplename + rawpart + '.'+graph_extension),
            dpi=graph_dpi)
        fig_curves.savefig(
            os.path.join(ip.user_ns['auximages_dir'],
                         'averaging1D_' +
                         samplename + rawpart + '.'+graph_extension),
            dpi=graph_dpi)
        putlogo(fig_correlmatrices)
        fig_correlmatrices.savefig(
            os.path.join(ip.user_ns['auximages_dir'],
                         'correlation_' +
                         samplename + rawpart + '.'+graph_extension),
            dpi=graph_dpi)
        writemarkdown("### Collected images from all distances")
        plt.show()
    writemarkdown("Updated badfsns list:")
    writemarkdown('['+', '.join(str(f) for f in ip.user_ns['badfsns'])+']')
    ip.user_ns['_data1d'] = data1d
    ip.user_ns['_data2d'] = data2d
    ip.user_ns['_headers_tosave'] = headers_tosave
    ip.user_ns['_rowavg']=rowavg

def _merge_two_curves(curve1, curve2, qmin, qmax, qsep, use_additive_constant=False):
    """Merge two scattering curves

    :param curve1: the first curve (longer distance)
    :type curve1: sastool.classes.curve.GeneralCurve
    :param curve2: the second curve (shorter distance)
    :type curve2: sastool.classes.curve.GeneralCurve
    :param qmin: lower bound of the interval for determining the scaling factor
    :type qmin: float
    :param qmax: upper bound of the interval for determining the scaling factor
    :type qmax: float
    :param qsep: separating (tailoring) point for the merge
    :type qsep: float
    :return: merged_curve, factor, background, stat
    :rtype tuple of a sastool.classes.curve.SASCurve and a float
    """
    if len(curve1.trim(qmin,qmax))>len(curve2.trim(qmin,qmax)):
        curve2_interp=curve2.trim(qmin,qmax)
        curve1_interp=curve1.interpolate(curve2_interp.q)
        print('Interpolated curve1 to curve2')
    else:
        curve1_interp=curve1.trim(qmin,qmax)
        curve2_interp=curve2.interpolate(curve1_interp.q)
        print('Interpolated curve2 to curve1')
    if use_additive_constant:
        bg_init=0
    else:
        bg_init=FixedParameter(0)
    factor, bg, stat=nonlinear_odr(curve2_interp.y, curve1_interp.y, curve2_interp.dy, curve1_interp.dy, lambda x,factor, bg:x*factor+bg,[1,0])
    return SASCurve.merge(curve1-bg, factor*curve2, qsep), factor, bg, stat


def unite(samplename, uniqmin=[], uniqmax=[], uniqsep=[], graph_ncols=2, graph_subplotpars={'hspace':0.3}, graph_extension='png', graph_dpi=80, additive_constant=False):
    ip = get_ipython()
    data1d = ip.user_ns['_data1d'][samplename]
    print("Uniting measurements of sample %s at different s-d distances" % samplename)
    uniparams = {'qmin': uniqmin, 'qmax': uniqmax, 'qsep': uniqsep}
    for p in uniparams:
        uniparams[p] = uniparams[p] + [None] * \
            max(0, len(data1d) - 1 - len(uniparams[p]))
    dists = list(reversed(sorted(data1d.keys())))
    if len(dists) < 2:
        print("Less than two distances found for sample %s; no point of uniting." % samplename)
        return
    united = None
    graph_nrows = int(
        np.ceil((len(dists)) / (graph_ncols * 1.0)))
    fig = plt.figure()
    unitedaxis = fig.add_subplot(graph_nrows, graph_ncols, 1)
    factor = 1
    for idx, dist1, dist2, qmin, qmax, qsep in zip(list(range(len(dists) - 1)),
                                                   dists[:-1], dists[1:],
                                                   uniparams['qmin'],
                                                   uniparams['qmax'],
                                                   uniparams['qsep']):
        print("    Scaling together distances %f and %f mm" % (dist1, dist2),flush=True)
        if united is None:
            united = data1d[dist1]
        if qmin is None:
            qmin = data1d[dist2].sanitize(minval=0,fieldname='q').q.min()
            print("        Auto-detected qmin:", qmin,flush=True)
        if qmax is None:
            qmax = data1d[dist1].sanitize(minval=0,fieldname='q').q.max()
            print("        Auto-detected qmax:", qmax,flush=True)
        if qsep is None:
            qsep = 0.5 * (qmin + qmax)
            print("        Auto-detected qsep:", qsep,flush=True)
        ax = fig.add_subplot(graph_nrows, graph_ncols, 2 + idx)
        (factor * data1d[dist1]).loglog(axes=ax, label='%.2f mm' % dist1)
        united, factor1, bg, stat = _merge_two_curves(united,
            data1d[dist2], qmin, qmax, qsep, use_additive_constant=additive_constant)
        factor=factor*factor1
        uniparams['qmin'][idx] = qmin
        uniparams['qmax'][idx] = qmax
        uniparams['qsep'][idx] = qsep
        print("        Scaling factor is", factor.tostring(), flush=True)
        if not additive_constant:
            print("        Additive constant has not been used.",flush=True)
        else:
            print("        Additive constant is:",bg.tostring(), flush=True)
        print("        Reduced Chi^2 of the ODR fit:",stat['Chi2_reduced'], flush=True)
        print("        DoF of the ODR fit:",stat['DoF'],flush=True)
        (factor * data1d[dist2]+bg).loglog(axes=ax, label='%.2f mm' % dist2)
        ax.set_xlabel('q (' + qunit() + ')')
        ax.set_ylabel('$d\\Sigma/d\\Omega$ (cm$^{-1}$ sr$^{-1}$)')
        ax.legend(loc='best')
        # ax.grid(which='both')
        ax.axis('tight')
        ax.set_title('Factor: ' + str(factor))
        lims = ax.axis()
        ax.plot([qmin, qmin], lims[2:], '--r', lw=2)
        ax.plot([qmax, qmax], lims[2:], '--r', lw=2)
        ax.plot([qsep, qsep], lims[2:], '--k')
    if '_data1dunited' not in ip.user_ns:
        ip.user_ns['_data1dunited'] = {}
    united.loglog(axes=unitedaxis)
    unitedaxis.set_xlabel('q (' + qunit() + ')')
    unitedaxis.set_ylabel('$d\\Sigma/d\\Omega$ (cm$^{-1}$ sr$^{-1}$)')
    unitedaxis.legend(loc='best')
    unitedaxis.set_title('United scattering of %s' % samplename)
    # unitedaxis.grid(which='both')
    unitedaxis.axis('tight')
    lims = unitedaxis.axis()
    for qs in uniparams['qsep']:
        unitedaxis.plot([qs] * 2, lims[2:], '--r')
    ip.user_ns['_data1dunited'][samplename] = united
    putlogo()
    fig.subplots_adjust(**graph_subplotpars)
    plt.savefig(
        os.path.join(ip.user_ns['auximages_dir'], 'uniting_' + samplename + '.'+ graph_extension), dpi=graph_dpi)
    print("    United curve spans the following ranges:")
    print("        q_min: ", print_abscissavalue(united.q.min(), ip.user_ns['_data2d'][samplename][dists[0]]['Wavelength']))
    print("        q_max: ", print_abscissavalue(united.q.max(), ip.user_ns['_data2d'][samplename][dists[0]]['Wavelength']))
    print("        q_max/q_min:", united.q.max() / united.q.min())
    print("        I_min: ", united.Intensity.min(), "cm^{-1}")
    print("        I_max: ", united.Intensity.max(), "cm^{-1}")
    print("        I_max/I_min:", united.Intensity.max() / united.Intensity.min())
    print("        # of points: ", len(united))
    united.save(os.path.join(ip.user_ns['saveto_dir'], 'united_' + samplename + '.txt'))
    plt.show()

def subtract_bg(samplename, bgname, factor=1, distance=None, disttolerance=2,
                subname=None, qrange=(), graph_extension='png', graph_dpi=80):
    """Subtract background from measurements.

    Inputs:
        samplename: the name of the sample
        bgname: the name of the background measurements. Alternatively, it can
            be a numeric value (float or ErrorValue), which will be subtracted.
            If None, this constant will be determined by integrating the
            scattering curve in the range given by qrange.
        factor: the background curve will be multiplied by this
        distance: if None, do the subtraction for all sample-to-detector distances.
            Otherwise give here the value of the sample-to-detector distance.
        qrange: a tuple (qmin, qmax)
        disttolerance: the tolerance in which two distances are considered
            equal.
        subname: the sample name of the background-corrected curve. The default
            is samplename + '-' + bgname
    """
    ip = get_ipython()
    data1d = ip.user_ns['_data1d']
    data2d = ip.user_ns['_data2d']
    if 'subtractedsamplenames' not in ip.user_ns:
        ip.user_ns['subtractedsamplenames']=set()
    subtractedsamplenames=ip.user_ns['subtractedsamplenames']
    if subname is None:
        if isinstance(bgname,str):
            subname=samplename+'-'+bgname
        else:
            subname=samplename+'-const'
    if distance is None:
        dists=data1d[samplename]
    else:
        dists=[d for d in data1d[samplename] if abs(d-distance)<disttolerance]
    for dist in dists:
        if isinstance(bgname,str):
            if not disttolerance:
                if dist not in data1d[bgname]:
                    print('Warning: Missing distance %g for background measurement (samplename: %s, background samplename: %s)'%(dist,samplename,bgname))
                    continue
                else:
                    bgdist=dist
            else:
                bgdist=sorted([(d,r) for (d,r) in [(d,np.abs(d-dist)) for d in list(data1d[bgname].keys())] if r<=disttolerance],key=lambda x:x[1])[0][0]
        if subname not in data1d:
            data1d[subname]={}
        if subname not in data2d:
            data2d[subname]={}
        data1_s=data1d[samplename][dist]
        data2_s=data2d[samplename][dist]
        if isinstance(bgname, str):
            data1_bg=data1d[bgname][bgdist]
            data2_bg=data2d[bgname][bgdist]
            if factor is None:
                factor=data1_s.trim(*qrange).momentum(0)/data1_bg.trim(*qrange).momentum(0)
        elif bgname is None:
            data1_bg=data1_s.trim(*qrange).momentum(0)
            data2_bg=data1_bg
        else:
            data1_bg=bgname
            data2_bg=bgname
        if factor is None:
            factor=1
        data1d[subname][dist]=data1_s-factor*data1_bg
        data2d[subname][dist]=data2_s-factor*data2_bg
        data1d[subname][dist].save(os.path.join(ip.user_ns['saveto_dir'], subname+'_' + ('%.2f'%dist).replace('.','_')+ '.txt'))
        plt.figure()
        plotsascurve(samplename,dist=dist)
        if isinstance(bgname, str):
            plotsascurve(bgname,dist=dist,factor=factor)
        plotsascurve(subname,dist=dist)
        plt.savefig(os.path.join(ip.user_ns['auximages_dir'],
                                 'uniting_' + samplename + '.'+graph_extension),
                    dpi=graph_dpi)

        subtractedsamplenames.add(subname)
