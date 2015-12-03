__all__=['restoredata']
from IPython.core.getipython import get_ipython
import os
import re
import numpy as np
from sastool.classes import SASHeader, SASCurve, SASExposure, SASMask

def restoredata(rootpath=None):
    ip = get_ipython()
    if rootpath is None:
        saveto_dir=ip.user_ns['saveto_dir']
    else:
        saveto_dir=os.path.join(rootpath,ip.user_ns['saveto_dir_rel'])

    if 'allsamplenames' not in ip.user_ns:
        ip.user_ns['allsamplenames']=set()
    allsamplenames=ip.user_ns['allsamplenames']
    if '_data2d' not in ip.user_ns:
        ip.user_ns['_data2d']={}
    if '_data1d' not in ip.user_ns:
        ip.user_ns['_data1d']={}
    if '_headers_tosave' not in ip.user_ns:
        ip.user_ns['_headers_tosave']={}
    if '_data1dunited' not in ip.user_ns:
        ip.user_ns['_data1dunited']={}
    _data2d=ip.user_ns['_data2d']
    _data1d=ip.user_ns['_data1d']
    _headers_tosave=ip.user_ns['_headers_tosave']
    _data1dunited=ip.user_ns['_data1dunited']
    dirlist= os.listdir(saveto_dir)
    basenames={os.path.splitext(f)[0] for f in dirlist}
    extns=['.npz','.txt','.log']
    datasetnames=[f for f in basenames if all([f+extn in dirlist for extn in extns])]
    remaining=[f for f in dirlist if os.path.splitext(f)[0] not in datasetnames]
    print('Loading complete datasets')
    for dsn in sorted(datasetnames):
        dsn_base,mangleddist1,mangleddist2=dsn.rsplit('_',2)
        if not (re.match(r'\d+',mangleddist1) and re.match(r'\d+',mangleddist2)):
            continue
        print('  ',dsn_base)
        if dsn_base not in _data2d:
            _data2d[dsn_base]={}
        if dsn_base not in _data1d:
            _data1d[dsn_base]={}
        if dsn_base not in _headers_tosave:
            _headers_tosave[dsn_base]={}
        h=SASHeader(os.path.join(saveto_dir, dsn+'.log'), plugin='B1 log')
        d=np.load(os.path.join(saveto_dir,dsn+'.npz'))
        c=SASCurve(os.path.join(saveto_dir, dsn+'.txt'))
        ex=SASExposure(d['Intensity'],d['Error'])
        ex.header=h
        m=SASMask(os.path.join(saveto_dir, ex.header['maskid']+'.mat'))
        ex.set_mask(m)
        _data2d[dsn_base][ex.header['Dist']]=ex
        _data1d[dsn_base][ex.header['Dist']]=c
        _headers_tosave[dsn_base][ex.header['Dist']]=h
        allsamplenames.add(dsn_base)
    print('Loading united curves')
    for fn in [fn_ for fn_ in remaining if fn_.startswith('united_') and fn_.endswith('.txt')]:
        _data1dunited[fn[7:-4]]=SASCurve(os.path.join(saveto_dir,fn))
        print('  ',fn[7:-4])
        remaining.remove(fn)
    print('Loading remaining .txt files')
    for fn in [fn_ for fn_ in remaining if fn_.endswith('.txt')]:
        if re.match(r'.*_\d+_\d\d\.txt$',fn) is not None:
            basename,distmangled1,distmangled2=fn.rsplit('.',1)[0].rsplit('_',2)
            dist=float(distmangled1)+float('0.'+distmangled2)
            print('  ',basename)
            if basename not in _data1d:
                _data1d[basename]={}
            _data1d[basename][dist]=SASCurve(os.path.join(saveto_dir,fn))
            remaining.remove(fn)
#    print('Remaining (not loaded): ',', '.join(remaining))
