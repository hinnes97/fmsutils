import setuptools
from setuptools import find_packages
#from setuptools import setup, find_packages
from numpy.distutils.core import Extension, setup

ext1 = Extension(name='fmsplot.fort_interp',
                 sources=['fmsplot/fort_interp/lib_array2.f90','fmsplot/fort_interp/fort_interp.f90'],
                 )


setup(name='fmsplot',
      version='1.0',
      description='Plot GCM output',
      author='Hamish Innes',
      author_email='hamish.innes@physics.ox.ac.uk',
      url='https://github.com/hinnes97/fmsplot',
      packages=find_packages(),
      #package_data={'fmsplot': ['fort_interp.cpython-37m-x86_64-linux-gnu.so']},)
      ext_modules=[ext1]
)

