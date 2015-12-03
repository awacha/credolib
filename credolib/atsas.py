__all__=['read_gnom_pr','execute_command','autorg', 'shanum','datgnom','dammif']
import numpy as np
import subprocess
import itertools
import tempfile
import os
from sastool.classes import GeneralCurve, SASCurve
from sastool.misc.errorvalue import ErrorValue

def read_gnom_pr(filename):
    with open(filename, encoding='utf-8') as f:
        l=f.readline()
        while 'Distance distribution' not in l:
            l=f.readline()
        data=[]
        while True:
            l=f.readline()
            if not l:
                break
            if not l.strip():
                continue
            try:
                data.append([float(f_) for f_ in l.strip().split()])
            except ValueError:
               pass
            except:
                raise

        return np.array(data)

def execute_command(cmd, input_to_command=None,eat_output=False, noprint=False):
    if isinstance(input_to_command,str):
        stdin=subprocess.PIPE
    else:
        stdin=input_to_command
    popen=subprocess.Popen(cmd,stdout=subprocess.PIPE,stderr=subprocess.PIPE,stdin=stdin)
    if (isinstance(input_to_command,str)):
        input_to_command=input_to_command.encode('utf-8')
    if isinstance(input_to_command, bytes):
        popen.stdin.write(input_to_command)
    lines_iterator=itertools.chain(popen.stdout,popen.stderr)
    resultinglines=[]
    for line in lines_iterator:
        if not noprint:
            if not eat_output:
                print(str(line[:-1], encoding='utf-8'), flush=True)
            else:
                print(".", end='', flush=True)
        resultinglines.append(str(line[:-1], encoding='utf-8'))
    return resultinglines

def autorg(filename, mininterval=None, qminrg=None, qmaxrg=None, noprint=True):
    """Execute autorg.

    Inputs:
        filename: either a name of an ascii file, or an instance of SASCurve.
        mininterval: the minimum number of points in the Guinier range
        qminrg: the maximum value of qmin*Rg. Default of autorg is 1.0
        qmaxrg: the maximum value of qmax*Rg. Default of autorg is 1.3
        noprint: if the output of autorg should be redirected to the null 
            device.

    Outputs:
        Rg as an ErrorValue
        I0 as an ErrorValue
        qmin: the lower end of the chosen Guinier range
        qmax: the upper end of the chosen Guinier range
        quality: the quality parameter, between 0 and 1
        aggregation: float, the extent of aggregation
    """
    if isinstance(filename, GeneralCurve):
        curve=filename
        with tempfile.NamedTemporaryFile('w+b',
                                         delete=False) as f:
            curve.save(f)
            filename=f.name
    cmdline=['autorg', filename, '-f', 'ssv']
    if mininterval is not None:
        cmdline.extend(['--mininterval',str(mininterval)])
    if qminrg is not None:
        cmdline.extend(['--sminrg',str(qminrg)])
    if qmaxrg is not None:
        cmdline.extend(['--smaxrg',str(qmaxrg)])
    result = execute_command(cmdline,noprint=noprint)
    Rg, dRg, I0, dI0, idxfirst, idxlast, quality, aggregation, filename = result[0].split(None, 8)
    try:
        curve
    except NameError:
        curve=SASCurve(filename)
    else:
        os.unlink(filename)    
    return ErrorValue(float(Rg),float(dRg)),ErrorValue(float(I0),float(dI0)),curve.q[int(idxfirst)-1],curve.q[int(idxlast)-1],float(quality),float(aggregation)

def datgnom(filename, Rg=None, noprint=True):
    if Rg is None:
        Rg, I0, idxfirst, idxlast, quality, aggregation = autorg(filename)
    execute_command(['datgnom', filename, '-r', '%f' % float(Rg)],
                    noprint=noprint)
    gnomoutputfilename = filename.rsplit('.',1)[0] + '.out'
    gnomdata = read_gnom_pr(gnomoutputfilename)
    return gnomdata

def dammif(gnomoutputfilename, prefix=None, mode='fast', symmetry='P1', N=None,
           noprint=True):
    if prefix is None:
        prefix='dammif_'+gnomoutputfilename.rsplit('.',1)[0]
    if N is None:
        execute_command(['dammif','--prefix=%s'%prefix, '--omit-solvent',
                         '--mode=%s'%mode, '--symmetry=%s'%symmetry,
                         '--unit=NANOMETER', gnomoutputfilename],
                        noprint=noprint)
        return prefix+'-1.pdb'
    else:
        ret=[]
        for i in range(N):
            execute_command(['dammif','--prefix=%s_%03d'%(prefix,i), '--omit-solvent',
                             '--mode=%s'%mode, '--symmetry=%s'%symmetry,
                             '--unit=NANOMETER', gnomoutputfilename],
                            noprint=noprint)
            ret.append('%s_%03d-1.pdb'%(prefix,i))
        return ret

def shanum(filename,dmax=None, noprint=True):
    """Execute the shanum program to determine the optimum qmax
    according to an estimation of the optimum number of Shannon
    channels.
    
    Inputs:
        filename: either a name of an ascii file, or an instance
            of SASCurve
        dmax: the cut-off of the P(r) function, if known. If None,
            this will be determined by the shanum program
        noprint: if the printout of the program is to be suppressed.
    
    Outputs: dmax, nsh, nopt, qmaxopt
        dmax: the cut-off of the P(r) function.
        nsh: the estimated number of Shannon channels
        nopt: the optimum number of Shannon channels
        qmaxopt: the optimum value of the high-q cutoff
    """
    if isinstance(filename, GeneralCurve):
        curve=filename
        with tempfile.NamedTemporaryFile('w+b', delete=False) as f:
            curve.save(f)
            filename=f.name
    cmdline=['shanum', filename]
    if dmax is not None:
        cmdline.append(str(float(dmax)))
    result=execute_command(cmdline, noprint=noprint)
    for l in result:
        l=l.strip()
        if l.startswith('Dmax='):
            dmax=float(l.split('=')[1])
        elif l.startswith('Smax='):
            qmax=float(l.split('=')[1])
        elif l.startswith('Nsh='):
            nsh=float(l.split('=')[1])
        elif l.startswith('Nopt='):
            nopt=float(l.split('=')[1])
        elif l.startswith('Sopt='):
            qmaxopt=float(l.split('=')[1])
            
    return dmax, nsh, nopt, qmaxopt