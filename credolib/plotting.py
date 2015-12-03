__all__=['plotsascurve','guinierplot','kratkyplot']
from .io import getsascurve
import matplotlib.pyplot as plt
from sastool.libconfig import qunit, dunit

def plotsascurve(samplename, *args, **kwargs):
    if 'dist' not in kwargs:
        kwargs['dist'] = None
    data1d, dist = getsascurve(samplename, kwargs['dist'])
    del kwargs['dist']
    if 'factor' in kwargs:
        factor=kwargs['factor']
        del kwargs['factor']
    else:
        factor=1
    if 'label' not in kwargs:
        if isinstance(dist, str):
            kwargs['label'] = samplename + ' ' + dist
        else:
            kwargs['label'] = samplename + ' %g mm' % dist
    if 'errorbar' in kwargs:
        errorbars = bool(kwargs['errorbar'])
        del kwargs['errorbar']
    else:
        errorbars = False
    if errorbars:
        ret = (data1d*factor).errorbar(*args, **kwargs)
        plt.xscale('log')
        plt.yscale('log')
    else:
        ret = (data1d*factor).loglog(*args, **kwargs)
    plt.xlabel('q (' + qunit() + ')')
    plt.ylabel('$d\\Sigma/d\\Omega$ (cm$^{-1}$ sr$^{-1}$)')
    plt.legend(loc='best')
    plt.grid(True, which='both')
    plt.axis('tight')
    return ret


def guinierplot(*args, **kwargs):
    """Make a Guinier plot. This is simply a wrapper around plotsascurve()."""
    ret=plotsascurve(*args, **kwargs)
    plt.xscale('power',exponent=2)
    plt.yscale('log')
    return ret


def kratkyplot(samplename, *args, **kwargs):
    if 'dist' not in kwargs:
        kwargs['dist'] = None
    data1d, dist = getsascurve(samplename, kwargs['dist'])
    del kwargs['dist']
    if 'factor' in kwargs:
        factor=kwargs['factor']
        del kwargs['factor']
    else:
        factor=1
    if 'label' not in kwargs:
        if isinstance(dist, str):
            kwargs['label'] = samplename + ' ' + dist
        else:
            kwargs['label'] = samplename + ' %g mm' % dist
    if 'errorbar' in kwargs:
        errorbars = bool(kwargs['errorbar'])
        del kwargs['errorbar']
    else:
        errorbars = False
    data1dscaled=data1d*factor
    if errorbars:
        if hasattr(data1dscaled, 'dx'):
            dx=data1dscaled.qError
            dy=(data1dscaled.Error ** 2 * data1dscaled.q ** 4 +
                data1dscaled.Intensity ** 2 * data1dscaled.qError ** 2
                * data1dscaled.q ** 2 * 4) ** 0.5
        else:
            dx=None
            dy=data1dscaled.Error
        ret = plt.errorbar(data1dscaled.q,
                           data1dscaled.q ** 2 * data1dscaled.Intensity,
                           dy, dx, *args, **kwargs)
    else:
        ret = plt.plot(data1dscaled.q,
                       data1dscaled.Intensity * data1dscaled.q ** 2,
                       *args, **kwargs)
    plt.xlabel('q (' + dunit() + ')')
    plt.ylabel('$q^2 d\\Sigma/d\\Omega$ (' +
               qunit() +
               '$^{-2}$ cm$^{-1}$ sr$^{-1}$)')
    plt.legend(loc='best')
    plt.grid(True, which='both')
    plt.axis('tight')
    return ret
