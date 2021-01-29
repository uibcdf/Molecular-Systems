from setuptools import setup, find_packages
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension
#import distutils.extension

extensions_list=[]

setup(
    name='OpenMolecularSystems',
    version='0.0.1',
    author='UIBCDF Lab',
    author_email='uibcdf@gmail.com',
    package_dir={'openmolecularsystems': 'openmolecularsystems'},
    packages=find_packages(),
    ext_modules=extensions_list,
    package_data={'openmolecularsystems': []},
    scripts=[],
    url='http://uibcdf.org',
    download_url ='https://github.com/uibcdf/OpenMolecularSystems',
    license='MIT',
    description="---",
    long_description="---",
)

