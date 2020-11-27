import setuptools
from setuptools import find_packages
#from setuptools import setup, find_packages
from numpy.distutils.core import Extension, setup

ext1 = Extension(name='fmsutils.fort_interp',
                 sources=['fmsutils/fort_interp/lib_array2.f90',
                          'fmsutils/fort_interp/fort_interp.f90'],
                 )


setup(name='fmsutils',
      version='1.0',
      description='Some interpolation and plotting tools for FMS output',
      author='Hamish Innes',
      author_email='hamish.innes@physics.ox.ac.uk',
      url='https://github.com/hinnes97/fmsutils',
      packages=find_packages(),
      ext_modules=[ext1]
)

