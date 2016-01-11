__all__=['read_gnom_pr','execute_command','autorg', 'shanum','datgnom','dammif','bodies']
import numpy as np
import subprocess
import itertools
import ipy_table
from IPython.display import display
import tempfile
import os
import re
from sastool.classes import GeneralCurve, SASCurve
from sastool.misc.errorvalue import ErrorValue

def read_gnom_pr(filename, get_metadata=False):
    metadata={}
    with open(filename, 'rt', encoding='utf-8') as f:
        l=f.readline()
        while 'Final results' not in l:
            l=f.readline()
        assert(not f.readline().strip()) # skip empty line
        assert(f.readline().strip()=='Parameter    DISCRP    OSCILL    STABIL    SYSDEV    POSITV    VALCEN')
        parameters={'DISCRP':{},'OSCILL':{}, 'STABIL':{}, 'SYSDEV':{}, 'POSITV':{}, 'VALCEN':{}}
        for i in range(6):
            line=f.readline().strip().split()
            if i==4:
                # this line contains only a dashed line: "- - - - - - etc."
                assert(all([l=='-' for l in line]))
                continue
            what=line[0]
            (parameters['DISCRP'][what], parameters['OSCILL'][what],
             parameters['STABIL'][what], parameters['SYSDEV'][what],
             parameters['POSITV'][what], parameters['VALCEN'][what]) = tuple([
                float(x) for x in line[1:]])
        assert(not f.readline().strip()) # skip empty line
        match=re.match(r'Angular\s+range\s+:\s+from\s+(?P<qmin>\d+\.\d+)\s+to\s+(?P<qmax>\d+\.\d+)', f.readline().strip())
        assert(match is not None)
        metadata['qmin']=float(match.groupdict()['qmin'])
        metadata['qmax']=float(match.groupdict()['qmax'])
        match=re.match(r'Real\s+space\s+range\s+:\s+from\s+(?P<dmin>\d+\.\d+)\s+to\s+(?P<dmax>\d+\.\d+)', f.readline().strip())
        assert(match is not None)
        metadata['dmin']=float(match.groupdict()['dmin'])
        metadata['dmax']=float(match.groupdict()['dmax'])
        assert(not f.readline().strip())
        match=re.match(r'Highest ALPHA \(theor\) :\s+(?P<highestalpha>\d+\.\d+E[+-]?\d+)',f.readline().strip())
        assert(match is not None)
        metadata['highestalpha']=float(match.groupdict()['highestalpha'])
        match=re.match(r'Current ALPHA\s+:\s+(?P<currentalpha>\d+\.\d+E[+-]\d+)\s+Rg :  (?P<Rg>\d+\.\d+E[+-]\d+)\s+I\(0\) :\s+(?P<I0>\d+\.\d+E[+-]\d+)', f.readline().strip())
        assert(match is not None)
        metadata['currentalpha']=float(match.groupdict()['currentalpha'])
        metadata['Rg_guinier']=float(match.groupdict()['Rg'])
        metadata['I0_guinier']=float(match.groupdict()['I0'])
        assert(not f.readline().strip()) # skip empty line
        match=re.match(r'Total  estimate : (?P<totalestimate>\d+\.\d+)\s+ which is \s+(?P<qualitystring>.*)\s+solution', f.readline().strip())
        assert(match is not None)
        metadata['totalestimate']=float(match.groupdict()['totalestimate'])
        metadata['qualitystring']=match.groupdict()['qualitystring']
        assert(not f.readline().strip()) # skip empty line
        assert(f.readline().strip().split()==['S','J','EXP','ERROR','J','REG','I','REG'])
        assert(not f.readline().strip()) # skip empty line
        s=[]
        sj=[]
        jexp=[]
        jerror=[]
        jreg=[]
        ireg=[]
        l=f.readline()
        while l.strip():
            terms=[float(x) for x in l.strip().split()]
            s.append(terms[0])
            ireg.append(terms[-1])
            if len(terms)>2:
                sj.append(terms[0])
                jexp.append(terms[1])
                jerror.append(terms[2])
                jreg.append(terms[3])
            l=f.readline()
        metadata['q']=np.array(s)
        metadata['qj']=np.array(sj)
        metadata['jexp']=np.array(jexp)
        metadata['jerror']=np.array(jerror)
        metadata['jreg']=np.array(jreg)
        metadata['ireg']=np.array(ireg)
        assert('Distance distribution  function of particle' == f.readline().strip())
        assert(not f.readline().strip()) # skip empty line
        assert(not f.readline().strip()) # skip empty line
        assert(f.readline().strip().split()==['R','P(R)','ERROR'])
        assert(not f.readline().strip()) # skip empty line

        data=[]
        while True:
            l=f.readline()
            if not l.strip():
                break
            if not l.strip():
                continue
            try:
                data.append([float(f_) for f_ in l.strip().split()])
            except ValueError:
                if 'Reciprocal space' in l:
                    break
            except:
                raise
        l=f.readline()
        match=re.match(r'Real space: Rg =\s+(?P<Rg>\d+\.\d+(E[+-]?\d+)?) \+- (?P<dRg>\d+\.\d+(E[+-]?\d+)?)\s+I\(0\) =\s+(?P<I0>\d+\.\d+(E[+-]?\d+)?) \+-\s+(?P<dI0>\d+\.\d+(E[+-]?\d+)?)',l.strip())
        assert(match is not None)
        metadata['Rg_gnom']=ErrorValue(float(match.groupdict()['Rg']), float(match.groupdict()['dRg']))
        metadata['I0_gnom']=ErrorValue(float(match.groupdict()['I0']), float(match.groupdict()['dI0']))
        if get_metadata:
            return (np.array(data), metadata)
        else:
            return (np.array(data),)

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
    gnomdata, metadata = read_gnom_pr(gnomoutputfilename, get_metadata=True)
    return gnomdata, metadata

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

def bodies(filename, bodytypes=None, prefix=None, fit_timeout=10,Ndummyatoms=2000,noprint=True):
    BODIES=['ellipsoid','rotation-ellipsoid','cylinder','elliptic-cylinder','hollow-cylinder','parallelepiped','hollow-sphere','dumbbell']
    if bodytypes is None:
        bodytypes=BODIES
    unknownbodies=[b for b in bodytypes if b not in BODIES]
    if unknownbodies:
        raise ValueError('Unknown body type(s): '+', '.join(unknownbodies))
    if isinstance(filename, GeneralCurve):
        curve=filename
        with tempfile.NamedTemporaryFile('w+b', delete=False) as f:
            curve.save(f)
            filename=f.name
        assert(prefix is not None)
    else:
        if prefix is None:
            prefix=filename.rsplit('.',1)[0]
    fittingresults={}
    for b in bodytypes:
        print('Fitting geometrical body %s'%b,flush=True)
        p=subprocess.Popen(['bodies'],stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            stdout,stderr=p.communicate(input=b'f\n%s\n%d\n\n\n\n\n\n\n%s\n'%(filename.encode('utf-8'),BODIES.index(b)+1,prefix.encode('utf-8')),timeout=fit_timeout)
        except subprocess.TimeoutExpired:
            print('Fitting timed out.')
            continue
        stdout=stdout.decode('utf-8')
        stderr=stderr.decode('utf-8')
        if stderr:
            print('Error: ',stderr,flush=True)
        printing_on=False
        parameter_recording_on=False
        bodyparameters=[]
        bodyparameternames=[]
        fittingresults[b]={}
        for s in stdout.split('\n'):
            if s.startswith(' Input file name'):
                printing_on=True
            if printing_on and not noprint:
                print(s,flush=True)
            if s.startswith(' Body type'):
                parameter_recording_on=True
            if s.startswith(' Parameter \'scale\''):
                parameter_recording_on=False
            if parameter_recording_on and s.startswith(' Parameter \''):
                bodyparameters.append(float(s.split(':')[1].strip()))
                bodyparameternames.append(s[s.index("'")+1:(s.index("'")+s[s.index("'")+1:].index("'")+1)])
            if s.startswith(' Expected Radius of Gyration'):
                fittingresults[b]['Rgexp']=float(s.split(':')[1].strip())
            elif s.startswith(' Expected I0'):
                fittingresults[b]['I0exp']=float(s.split(':')[1].strip())
            elif s.startswith(' Expected Volume'):
                fittingresults[b]['Volexp']=float(s.split(':')[1].strip())
            elif s.startswith(' Fit Radius of Gyration'):
                fittingresults[b]['Rgfit']=float(s.split(':')[1].strip())
            elif s.startswith(' Fit I0'):
                fittingresults[b]['I0fit']=float(s.split(':')[1].strip())
            elif s.startswith(' Fit Volume'):
                fittingresults[b]['Volfit']=float(s.split(':')[1].strip())
            elif s.startswith(' Goodness of Fit (chi-square)'):
                fittingresults[b]['Chi2']=float(s.split(':')[1].strip())
        fittingresults[b]['type']=b
        fittingresults[b]['bodyparameters']=bodyparameters
        fittingresults[b]['bodyparameternames']=bodyparameternames
        print('Creating DAM model')
        damoutputfile=prefix+'-'+b+'.pdb'
        p=subprocess.Popen(['bodies'],stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            stdout,stderr=p.communicate(input=b'd\n%d\n'%(BODIES.index(b)+1)+b'\n'.join([b'%6f'%(10*v) for v in bodyparameters])+b'\n1\n%d\n%s\n'%(Ndummyatoms,damoutputfile.encode('utf-8')),timeout=fit_timeout)
        except subprocess.TimeoutExpired:
            print('Error creating DAM model.')
        if stderr:
            print(stderr)
    tab=[['Body','Goodness of Fit ($\chi^2$)','Rg mismatch','I0 mismatch','Volume mismatch']]
    for b in sorted(fittingresults):
        tab.append([
            fittingresults[b]['type']+' ('+', '.join(['%s=%.3f nm'%(var,val) for var,val in zip(fittingresults[b]['bodyparameternames'],
                                                                                                fittingresults[b]['bodyparameters'])])+')',
            fittingresults[b]['Chi2'],
            '%.2f nm'%(fittingresults[b]['Rgfit']-fittingresults[b]['Rgexp']),
            '%5g cm$^{-1}$ sr$^{-1}$'%(fittingresults[b]['I0fit']-fittingresults[b]['I0exp']),
            '%.2f nm^3'%(fittingresults[b]['Volfit']-fittingresults[b]['Volexp']),
        ])
    tab=ipy_table.IpyTable(tab)
    tab.apply_theme('basic')
    display(tab)
    return fittingresults

