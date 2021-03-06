__all__ = ['assess_flux_stability', 'sum_measurement_times', 'assess_sample_stability',
           'assess_instrumental_background', 'assess_transmission', 'assess_gc_fit', 'assess_fitting_results']

import ipy_table
import matplotlib.pyplot as plt
import numpy as np
from IPython.core.getipython import get_ipython
from IPython.display import display
from matplotlib.colors import LogNorm
from sastool.classes2 import Curve, Exposure
from sastool.misc.cormap import cormaptest
from sastool.misc.easylsq import nonlinear_leastsquares

from .io import load_exposure, load_mask


def assess_flux_stability(samplename='Glassy_Carbon'):
    ip = get_ipython()
    f = plt.figure()
    ax1 = f.add_subplot(1, 1, 1)
    plt.xlabel('Date of exposure')
    plt.ylabel('Beam flux (photon/sec), continuous lines')
    ax2 = plt.twinx()
    plt.ylabel('Vacuum pressure (mbar), dotted lines')
    plt.title('Beam flux stability')
    samplenames = sorted([sn_ for sn_ in ip.user_ns['_headers_sample'] if samplename in sn_])
    linestyles = ['b', 'g', 'r', 'c', 'm', 'y', 'k']
    lines = []
    for sn, ls in zip(samplenames, linestyles):
        print(sn)
        heds = ip.user_ns['_headers_sample'][sn]
        allheds = []
        for k in heds.keys():
            allheds.extend(heds[k])
        allheds = sorted(allheds, key=lambda x: x.fsn)
        flux = np.array([float(h.flux) for h in allheds])
        dates = [h.date for h in allheds]
        lines.extend(ax1.plot(dates, flux, ls + 'o', label='Flux (%s)' % sn))
        vacuums = np.array([float(h.vacuum) for h in allheds])
        lines.extend(ax2.plot(dates, vacuums, ls + 's', label='Vacuum (%s)' % sn, lw=2))
        print('  Measurement duration: %.2f h' % ((dates[-1] - dates[0]).total_seconds() / 3600.))
        print('  Mean flux: ', flux.mean(), '+/-', flux.std(), 'photons/sec')
        print('  RMS variation of flux: ', flux.std() / flux.mean() * 100, '%')
        print('  P-P variation of flux: ', flux.ptp() / flux.mean() * 100, '%')
    ax1.legend(lines, [l.get_label() for l in lines], loc='best')
    plt.show()


def sum_measurement_times():
    ip = get_ipython()
    headers = ip.user_ns['_headers_sample']
    tab = [['Sample name', 'Distance (mm)', 'Measurement time (h)']]
    for samplename in sorted(headers):
        for dist in sorted(headers[samplename]):
            tab.append([samplename, "%.2f" % dist,
                        '%.2f' % (sum([h.exposuretime for h in headers[samplename][dist]]) / 3600)])
    tab = ipy_table.IpyTable(tab)
    tab.apply_theme('basic')
    display(tab)


def assess_sample_stability(end_cutoff=3):
    ip = get_ipython()
    rowavg = ip.user_ns['_rowavg']
    tab = [['Sample name', 'Distance', 'Slope of autocorrelation function', 'Stability']]
    plt.figure()
    for sn in sorted(rowavg):
        for dist in sorted(rowavg[sn]):
            rowavg_rescaled = rowavg[sn][dist] / rowavg[sn][dist].mean()
            rowavg_std = rowavg_rescaled.std() * np.ones_like(rowavg_rescaled)
            try:
                A, B, stat = nonlinear_leastsquares(np.arange(len(rowavg_rescaled)), rowavg_rescaled, rowavg_std,
                                                    lambda x, a, b: a * x + b, [0, 0])
                problematic = (A.val > A.err * 3)
            except TypeError:
                A = 'N/A'
                problematic = 2
            tab.append([sn, '%.2f' % dist, A, ["\u2713", "\u2718\u2718\u2718\u2718\u2718", '\u274e'][problematic]])
            plt.errorbar(np.arange(len(rowavg_rescaled)), rowavg_rescaled,
                         label=sn + ' %.2f mm' % dist)  # ,diags_std[1:])
    plt.xlabel('Separation in time (FSN units)')
    plt.ylabel('Average discrepancy between curves')
    plt.legend(loc='best')
    tab = ipy_table.IpyTable(tab)
    tab.apply_theme('basic')
    display(tab)
    plt.show()


def assess_instrumental_background(Wx=20, Wy=20, emptyname='Empty_Beam', maskname='mask.mat'):
    ip = get_ipython()
    for dist in ip.user_ns['_headers_sample'][emptyname]:
        data = ip.user_ns['_data2d'][emptyname][dist]
        assert isinstance(data, Exposure)
        intensity = data.intensity / data.header.exposuretime
        mask = load_mask(maskname).astype(np.bool)
        m = np.zeros_like(intensity)
        print('   Mean intensity per pixel:', intensity[mask].mean(), 'cps')
        print('   STD intensity per pixel:', intensity[mask == 1].std(), 'cps')
        print('   Total intensity:', intensity[mask == 1].sum(), 'cps')
        for row in range(m.shape[0]):
            for col in range(m.shape[1]):
                m[row, col] = intensity[max(row - Wy, 0):min(row + Wy, m.shape[0] - 1),
                              max(col - Wx, 0):min(col + Wx, m.shape[1] - 1)][
                    mask[max(row - Wy, 0):min(row + Wy, m.shape[0] - 1),
                    max(col - Wx, 0):
                    min(col + Wx, m.shape[1] - 1)] == 1].mean()
        plt.figure()
        plt.subplot(1, 2, 1)
        plt.imshow(intensity, norm=LogNorm())
        plt.subplot(1, 2, 2)
        plt.imshow(m)
        plt.tight_layout()
        plt.suptitle('Empty beam, {} mm'.format(dist))


def assess_transmission():
    ip = get_ipython()
    tab = [
        ['Sample name', 'Distance', 'Transmission', 'Linear absorption coefficient (1/cm)', 'Absorption length (cm)']]
    for sample in sorted(ip.user_ns['_headers_sample']):
        for dist in sorted(ip.user_ns['_headers_sample'][sample]):
            transms_seen = []
            for h in ip.user_ns['_headers_sample'][sample][dist]:
                if float(h.transmission) not in transms_seen:
                    transms_seen.append(float(h.transmission))
                    transm = h.transmission
                    thickness = h.thickness
                    try:
                        mu = (-transm.log()) / thickness
                    except ZeroDivisionError:
                        mu = 'Infinite'
                        invmu = 0
                    else:
                        try:
                            invmu = (1 / mu).tostring(extra_digits=2)
                        except ZeroDivisionError:
                            invmu = 'Infinite'
                        mu = mu.tostring(extra_digits=2)
                    tab.append([sample, dist, transm.tostring(extra_digits=2), mu, invmu])
    tab = ipy_table.IpyTable(tab)
    tab.apply_theme('basic')
    display(tab)


def assess_gc_fit(reffile=None, gcname='Glassy_Carbon'):
    ip = get_ipython()
    if reffile is None:
        reffile = ip.user_ns['_loaders'][0].get_subpath('config/GC_data_nm.dat')
    refcurve = Curve.new_from_file(reffile)
    f = plt.figure()
    f.add_subplot(1, 1, 1)
    rads = {}
    for fsn in sorted([h.fsn for h in ip.user_ns['_headers']['processed'] if h.title == gcname]):
        try:
            ex = load_exposure(fsn, raw=False, processed=True)
        except:
            continue
        rads[ex.header.fsn] = ex.radial_average(refcurve.q)
        del ex
    qmin = max([r.sanitize().q.min() for r in rads.values()])
    qmax = min([r.sanitize().q.max() for r in rads.values()])
    refcurve.trim(qmin, qmax).loglog('o', mfc='none', ms=10)
    for r in sorted(rads):
        rads[r].loglog('.', label='#{:d}'.format(r))
    plt.axis('tight')
    plt.legend(loc='best', numpoints=1)
    plt.xlabel('q (nm$^{-1}$)')
    plt.ylabel('$d\Sigma/d\Omega$ (cm$^{-1}$ sr$^{-1}$)')
    plt.grid(True, which='both')
    plt.draw()


def calc_chi2(y, dy, fittedy):
    return (((y - fittedy) / dy) ** 2).sum() / (len(y) - 1)


def calc_R2(y, fittedy):
    SStot = ((y - np.mean(y)) ** 2).sum()
    SSres = ((fittedy - y) ** 2).sum()
    return 1 - SSres / SStot


def assess_fitting_results(basename, cormap_alpha=0.01):
    """Assess the results of a fit based on the .fit and .fir files created by
    various programs from the ATSAS suite."""
    plt.figure(figsize=(12, 4))
    plt.subplot2grid((1, 4), (0, 0), colspan=2)
    fir = np.loadtxt(basename + '.fir', skiprows=1)  # q, Iexp, Errexp, Ifitted
    # do a cormap test to compare the raw data and the model.
    pvalf, Cf, cormapf = cormaptest(fir[:, 1], fir[:, 3])
    cormapstatusf = ['Reject', 'Accept'][pvalf >= cormap_alpha]

    plt.errorbar(fir[:, 0], fir[:, 1], fir[:, 2], None, 'bo-', label='Raw data')
    plt.plot(fir[:, 0], fir[:, 3], 'r-', label='Fitted')
    chi2 = calc_chi2(fir[:, 1], fir[:, 2], fir[:, 3])
    R2 = calc_R2(fir[:, 1], fir[:, 3])
    try:
        skiprows = 0
        while True:
            try:
                fit = np.loadtxt(basename + '.fit', skiprows=skiprows)  # q, Ismoothed, Ifitted
                break
            except ValueError as ve:
                if ve.args[0].startswith('could not convert string to float'):
                    skiprows += 1
                    continue
                else:
                    raise
        # do a cormap test to compare the raw data to the smoothed data
        smoothed = fit[(fit[:, 0] >= fir[:, 0].min()) & (fit[:, 0] <= fir[:, 0].max()), 1]
        pvals, Cs, cormaps = cormaptest(fir[:, 1], smoothed)
        cormapstatuss = ['Reject', 'Accept'][pvals >= cormap_alpha]
        plt.plot(fit[:, 0], fit[:, 1], 'g.-', label='Smoothed, extrapolated')
        plt.plot(fit[:, 0], fit[:, 2], 'm-', label='Fitted to smoothed, extrapolated')
    except ValueError as ve:
        print('Error while loading file: {}.fit: {}'.format(basename, ve))
    except FileNotFoundError:
        fit = None
        cormaps = cormapstatuss = pvals = Cs = None
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='best')
    plt.grid(which='both')
    if fit is not None:
        plt.subplot2grid((1, 4), (0, 2))
        plt.imshow(cormaps, cmap='gray', interpolation='nearest')
        plt.title('CorMap of the smoothing')
    plt.subplot2grid((1, 4), (0, 3))
    plt.imshow(cormapf, cmap='gray', interpolation='nearest')
    plt.title('CorMap of the fitting')
    print('R2: ', R2)
    print('Chi2: ', chi2)
    if fit is not None:
        print('Cormap test of the smoothing: {} (p={}, C={}, N={})'.format(cormapstatuss, pvals, Cs, cormaps.shape[0]))
    print('Cormap test of fit: {} (p={}, C={}, N={})'.format(cormapstatusf, pvalf, Cf, cormapf.shape[0]))
