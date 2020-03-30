# -*- coding: utf-8 -*-

from setuptools import setup, find_packages
from distutils.core import setup, Extension

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()
    
setup(
    name='doris',
    version='1.0.0',
    description='Supporting code for DORIS.',
    long_description=readme,
    author='James Tuck',
    author_email='jtuck@ncsu.edu',
    url='https://github.com/jamesmtuck/DORIS',
    license=license,
    packages=find_packages(exclude=('scripts','data','docs')),
    ext_modules= []
)
