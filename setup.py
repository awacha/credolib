#!/usb/bin/env python
import os

from setuptools import setup, find_packages

def getresourcefiles():
    print('Generating resource list',flush=True)
    reslist=[]
    for directory, subdirs, files in os.walk('credolib/resource'):
        reslist.extend([os.path.join(directory,f).split('/',1)[1] for f in files])
    print('Generated resource list:\n  '+'\n  '.join(x for x in reslist)+'\n',flush=True)
    return reslist


setup(name='credolib', author='Andras Wacha',
      author_email='awacha@gmail.com', url='http://github.com/awacha/credolib',
      description='CREDO Data processing library',
      packages=find_packages(),
      install_requires=['numpy>=1.0.0', 'scipy>=0.7.0', 'matplotlib', 'sastool', 'IPython', 'ipy_table', 'ipywidgets'],
      keywords="saxs sans sas small-angle scattering x-ray instrument control",
      setup_requires=['setuptools_scm'],
      use_scm_version=True,
      license="BSD 3-clause",
      package_data={'': getresourcefiles()},
      zip_safe=False,
      )
