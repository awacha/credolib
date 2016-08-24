__all__ = ['load_headers', 'getsascurve', 'getsasexposure', 'getheaders', 'getdists', 'filter_headers']
from IPython.core.getipython import get_ipython


def filter_headers(criterion):
    """Filter already loaded headers against some criterion.

    The criterion function must accept a single argument, which is an instance
    of sastool.classes2.header.Header, or one of its subclasses. The function
    must return True if the header is to be kept or False if it needs to be
    discarded. All manipulations on the header (including sample name changes,
    etc.) carried out by this function are preserved.
    """
    ip = get_ipython()
    for headerkind in ['processed', 'raw']:
        for h in ip.user_ns['_headers'][headerkind][:]:
            if not criterion(h):
                ip.user_ns['_headers'][headerkind].remove(h)

def load_headers(fsns):
    """Load header files
    """
    ip = get_ipython()
    ip.user_ns['_headers'] = {}
    for type_ in ['raw', 'processed']:
        print("Loading %d headers (%s)" % (len(fsns), type_), flush=True)
        processed = type_ == 'processed'
        headers = []
        for f in fsns:
            for l in [l_ for l_ in ip.user_ns['_loaders'] if l_.processed == processed]:
                try:
                    headers.append(l.loadheader(f))
                    break
                except FileNotFoundError:
                    continue
        allsamplenames = {h.title for h in headers}
        if not headers:
            print('NO HEADERS READ FOR TYPE "%s"' % type_)
        else:
            print("%d headers (%s) out of %d have been loaded successfully." % (len(headers), type_, len(fsns)))
            print('Read FSN range:', min([h.fsn for h in headers]), 'to', max([h.fsn for h in headers]))
            print("Samples covered by these headers:")
            print("    " + "\n    ".join(sorted(allsamplenames)), flush=True)
        if processed:
            ip.user_ns['allsamplenames'] = allsamplenames
        ip.user_ns['_headers'][type_] = headers

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


def getheaders(processed=True):
    ip = get_ipython()
    if processed:
        return ip.user_ns['_headers']['processed']
    else:
        return ip.user_ns['_headers']['raw']

def getdists(samplename):
    ip = get_ipython()
    return sorted([d for d in ip.user_ns['_headers_sample'][samplename]])

def get_different_distances(headers, tolerance=2):
    alldists = {float(h.distance) for h in headers}
    dists = []
    for d in alldists:
        if [d_ for d_ in dists if abs(d - d_) < tolerance]:
            continue
        dists.append(d)
    return sorted(dists)
