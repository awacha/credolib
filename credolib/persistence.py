__all__ = ['restoredata', 'storedata']
import pickle
import warnings

from IPython.core.getipython import get_ipython


def storedata(filename=None):
    """Store the state of the current credolib workspace in a pickle file."""
    if filename is None:
        filename = 'credolib_state.pickle'
    ns = get_ipython().user_ns
    with open(filename, 'wb') as f:
        d = {}
        for var in ['_headers', '_loaders', '_data1d', '_data2d', '_data1dunited',
                    'allsamplenames', '_headers_sample', 'badfsns', '_rowavg',
                    'saveto_dir', 'mask_override',
                    'badfsns_datcmp', 'auximages_dir', 'subtractedsamplenames', 'outputpath', 'saveto_dir_rel',
                    'auximages_dir_rel', 'crd_prefix']:
            try:
                d[var] = ns[var]
            except KeyError:
                warnings.warn('Skipping storage of unavailable variable "%s"' % var)
        pickle.dump(d, f)


def restoredata(filename=None):
    """Restore the state of the credolib workspace from a pickle file."""
    if filename is None:
        filename = 'credolib_state.pickle'
    ns = get_ipython().user_ns
    with open(filename, 'rb') as f:
        d = pickle.load(f)
    for k in d.keys():
        ns[k] = d[k]
