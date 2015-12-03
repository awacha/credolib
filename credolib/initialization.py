__all__=['set_length_units','init_dirs']

from sastool.misc import find_subdirs
import sastool.libconfig as libconfig
import os
from IPython.core.getipython import get_ipython

def set_length_units(units):
    """Set the length units: either 'nm' or 'A'.
    """
    libconfig.LENGTH_UNIT = units
    print("Length units have been set to:", units,flush=True)

def init_dirs(rootdir, outputpath, saveto_dir='data',
              auximages_dir='auximages', prefix='crd'):
    """Initialize the directiories.

    Inputs:
        rootdir: the root directory of the SAXSCtrl software,
            i.e. where the subfolders ``eval2d``, ``param``,
            ``images``, ``mask`` etc. reside.
        outputpath: the directory where the produced files are
            written. This is usually the working directory of
            the IPython notebook.
        saveto_dir: the subdirectory where averaged, united,
            subtracted etc. datasets are written.
        auximages_dir: the subdirectory where automatically produced
            images reside.
    """
    ip = get_ipython()
    ip.user_ns['credo_root'] = rootdir
    ip.user_ns['evaldirs'] =  [os.path.join(rootdir, 'eval2d'),
                               os.path.join(rootdir, 'eval2d', 'crd'),
                               os.path.join(rootdir, 'eval1d'),
                               os.path.join(rootdir, 'eval1d', 'crd')] + \
        find_subdirs(os.path.join(rootdir, 'mask'))
    ip.user_ns['datadirs'] = [os.path.join(rootdir, 'param_override'),
                              os.path.join(rootdir, 'images'),
                              os.path.join(rootdir, 'images', prefix),
                              os.path.join(rootdir, 'param'),
                              os.path.join(rootdir, 'param', prefix)] + \
        find_subdirs(os.path.join(rootdir, 'mask'))
    ip.user_ns['onedim_folder'] = os.path.join(rootdir, 'eval1d')
    if not os.path.isabs(outputpath):
        outputpath = os.path.join(rootdir, 'processing', outputpath)
    if not os.path.isdir(outputpath):
        os.makedirs(outputpath)
    if not os.path.isdir(os.path.join(ip.user_ns['outputpath'], saveto_dir)):
        os.mkdir(os.path.join(ip.user_ns['outputpath'], saveto_dir))
    if not os.path.isdir(os.path.join(ip.user_ns['outputpath'], auximages_dir)):
        os.mkdir(os.path.join(ip.user_ns['outputpath'], auximages_dir))
    print("Output files will be written to:", outputpath)
    os.chdir(outputpath)
    ip.user_ns['outputpath'] = outputpath
    ip.user_ns['auximages_dir'] = os.path.join(outputpath, auximages_dir)
    ip.user_ns['saveto_dir'] = os.path.join(outputpath, saveto_dir)
    ip.user_ns['saveto_dir_rel'] = saveto_dir
    ip.user_ns['auximages_dir_rel'] = auximages_dir
    ip.user_ns['crd_prefix']=prefix
    set_length_units('nm')