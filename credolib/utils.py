__all__=['writemarkdown','putlogo','print_abscissavalue']
from IPython.display import display,Markdown
from IPython.core.getipython import get_ipython
import matplotlib.pyplot as plt
import sastool
import numpy as np
import pkg_resources
from scipy.misc import imread

credo_logo = imread(pkg_resources.resource_filename('credolib','resource/credo_logo.png'))


def writemarkdown(*args):
    display(Markdown(' '.join(str(a) for a in args)))

def putlogo(figure=None):
    """Puts the CREDO logo at the bottom right of the current figure (or
    the figure given by the ``figure`` argument if supplied).
    """
    ip = get_ipython()
    if figure is None:
        figure=plt.gcf()
    curraxis= figure.gca()
    logoaxis = figure.add_axes([0.89, 0.01, 0.1, 0.1], anchor='NW')
    logoaxis.set_axis_off()
    logoaxis.xaxis.set_visible(False)
    logoaxis.yaxis.set_visible(False)
    logoaxis.imshow(credo_logo)
    figure.subplots_adjust(right=0.98)
    figure.sca(curraxis)

def print_abscissavalue(q, wavelength=None, distance=None, digits=10):
    qunit = sastool.libconfig.qunit()
    dunit = sastool.libconfig.dunit()
    formatstring='%%.%df'%digits
    retval = str(q) + ' ' + qunit
    retval = retval + "("
    retval = retval + " <=> " + formatstring %(2 * np.pi / q) + " " + dunit + "(d)"
    retval = retval + " <=> " + formatstring %(1 / q) + " " + dunit + "(Rg)"
    if wavelength is not None:
        tth_rad = 2 * np.arcsin((q * wavelength) / 4 / np.pi)
        tth_deg = tth_rad * 180.0 / np.pi
        retval = retval + " <=> " + formatstring %(tth_deg) + "\xb0"
        if distance is not None:
            radius = np.tan(tth_rad) * distance
            retval = retval + " <=> " + formatstring % (radius) + " mm(r)"
    retval = retval + ")"
    return retval

class figsize(object):
    def __init__(self, sizex, sizey):
        self._originalsize=plt.rcParams['figure.figsize']
        plt.rcParams['figure.figsize']=(sizex, sizey)
    def __enter__(self):
        pass
    def __exit__(self, exc_type, exc_val, exc_tb):
        plt.rcParams['figure.figsize']=self._originalsize
        return False # we don't want to suppress the exception, if any
