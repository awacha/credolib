__all__=['load_headers','getsascurve','getsasexposure','getheaders', 'getdists']
from IPython.core.getipython import get_ipython
from sastool.classes import SASHeader

def load_headers(fsns, skip_errorflags=True):
    """Load header files
    """
    ip = get_ipython()
    for type_, dirs, varname in [('eval', ip.user_ns['evaldirs'], '_evalheaders'), ('raw', ip.user_ns['datadirs'], '_rawheaders')]:
        print("Loading %d headers (%s)" % (len(fsns), type_), flush=True)
        headers = SASHeader(
            ip.user_ns['crd_prefix']+'_%05d.param', fsns, dirs=dirs, error_on_not_found=False)
        if skip_errorflags:
            headers = [h for h in headers if not h['ErrorFlags']]
        print("%d headers (%s) out of %d have been loaded successfully." % (len(headers), type_, len(fsns)))
        print('Read FSN range:', min([h['FSN'] for h in headers]), 'to', max([h['FSN'] for h in headers]))
        allsamplenames = {h['Title'] for h in headers}
        print("Samples covered by these headers:")
        print("    " + "\n    ".join(sorted(allsamplenames)),flush=True)
        if varname == '_evalheaders':
            ip.user_ns['allsamplenames'] = allsamplenames
        ip.user_ns[varname] = headers

def getsascurve(samplename, dist=None):
    ip = get_ipython()
    if dist == 'united':
        data1d = ip.user_ns['_data1dunited'][samplename]
    elif dist is None:
        try:
            data1d = ip.user_ns['_data1dunited'][samplename]
            dist = 'united'
        except KeyError:
            data1d = ip.user_ns['_data1d'][samplename]
            dist = sorted(data1d.keys())[0]
            data1d = data1d[dist]
    else:
        data1d = ip.user_ns['_data1d'][samplename]
        dist = sorted(list(data1d.keys()), key=lambda k:abs(float(dist) - k))[0]
        data1d = data1d[dist]
    return data1d, dist

def getsasexposure(samplename, dist=None):
    ip = get_ipython()
    if dist is None:
        data2d = ip.user_ns['_data2d'][samplename]
        dist = sorted(data2d.keys())[0]
        data2d = data2d[dist]
    else:
        data2d = ip.user_ns['_data2d'][samplename]
        dist = sorted(list(data2d.keys()), key=lambda k:abs(float(dist) - k))[0]
        data2d = data2d[dist]
    return data2d, dist

def getheaders(eval=True):
    ip = get_ipython()
    if eval:
        return ip.user_ns['_evalheaders']
    else:
        return ip.user_ns['_rawheaders']

def getdists(samplename):
    ip = get_ipython()
    return sorted([d for d in ip.user_ns['_data2d'][samplename]])

def get_different_distances(headers, tolerance=2):
    alldists = {h['DistCalibrated'] for h in headers}
    dists = []
    for d in alldists:
        if [d_ for d_ in dists if abs(d - d_) < tolerance]:
            continue
        dists.append(d)
    return sorted(dists)
