__all__=['set_length_units','init_dirs']

import os

import sastool.libconfig as libconfig
from IPython.core.getipython import get_ipython
from sastool.classes2 import Loader
from sastool.io import credo_saxsctrl, credo_cct

def set_length_units(units):
    """Set the length units: either 'nm' or 'A'.
    """
    libconfig.LENGTH_UNIT = units
    print("Length units have been set to:", units,flush=True)


def init_dirs(rootdir_or_loader, outputpath, saveto_dir='data',
              auximages_dir='auximages', prefix='crd'):
    """Initialize the directiories.

    Inputs:
        rootdir_or_loader: depends on the type:
            str:
                the root directory of the SAXSCtrl/CCT
                software, i.e. where the subfolders ``eval2d``, ``param``,
                ``images``, ``mask`` etc. reside.
            sastool.classes2.Loader instance:
                a fully initialized loader, which will be used to acquire
                headers and exposures.
            list:
                a list of sastool.classes2.Loader instances, which will
                be used to open headers and exposures. When opening something,
                always the first item will be tried first, and if it fails with
                FileNotFoundError, the second, third, etc. will be tried until
                either the file can be opened or the last one fails.

        outputpath: the directory where the produced files are
            written. This is usually the working directory of
            the IPython notebook.

        saveto_dir: the subdirectory where averaged, united,
            subtracted etc. datasets are written.

        auximages_dir: the subdirectory where automatically produced
            images reside.

    Remarks:
        If a single root directory is given, a list of four loaders will be
        constructed in this order: CCT (processed), CCT (raw), SAXSCtrl (processed),
        SAXSCtrl (raw). Raw and processed loaders are handled separately.
    """
    ip = get_ipython()
    if isinstance(rootdir_or_loader, str):
        ip.user_ns['_loaders'] = [
            credo_cct.Loader(rootdir_or_loader, prefix),
            credo_saxsctrl.Loader(rootdir_or_loader, prefix),
        ]
    elif isinstance(rootdir_or_loader, Loader):
        ip.user_ns['_loaders'] = [rootdir_or_loader]
    elif isinstance(rootdir_or_loader, list) and all([isinstance(l, Loader) for l in rootdir_or_loader]):
        ip.user_ns['_loaders'] = rootdir_or_loader[:]
    else:
        raise TypeError(rootdir_or_loader)
    if not os.path.isdir(outputpath):
        os.makedirs(outputpath)
    print("Output files will be written to:", outputpath)
    os.chdir(outputpath)
    ip.user_ns['outputpath'] = outputpath
    if not os.path.isdir(os.path.join(ip.user_ns['outputpath'], saveto_dir)):
        os.mkdir(os.path.join(ip.user_ns['outputpath'], saveto_dir))
    if not os.path.isdir(os.path.join(ip.user_ns['outputpath'], auximages_dir)):
        os.mkdir(os.path.join(ip.user_ns['outputpath'], auximages_dir))
    ip.user_ns['auximages_dir'] = os.path.join(outputpath, auximages_dir)
    ip.user_ns['saveto_dir'] = os.path.join(outputpath, saveto_dir)
    ip.user_ns['saveto_dir_rel'] = saveto_dir
    ip.user_ns['auximages_dir_rel'] = auximages_dir
    ip.user_ns['crd_prefix']=prefix
    set_length_units('nm')