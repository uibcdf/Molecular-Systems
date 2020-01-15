from setuptools import setup, find_packages
from numpy.distutils.core import setup
from numpy.distutils.extension import Extension
#import distutils.extension

extensions_list=[]

setup(
    name='UIBCDF_test_systems',
    version='0.0.1',
    author='UIBCDF Lab',
    author_email='uibcdf@gmail.com',
    package_dir={'uibcdf_test_systems': 'uibcdf_test_systems'},
    packages=find_packages(),
    ext_modules=extensions_list,
    package_data={'uibcdf_test_systems': []},
    scripts=[],
    url='http://uibcdf.org',
    download_url ='https://github.com/uibcdf/UIBCDF_test_systems',
    license='MIT',
    description="---",
    long_description="---",
)

