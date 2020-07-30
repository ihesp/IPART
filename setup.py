#!/usr/bin/python

# Image-Processing based Atmospheric River Tracking (IPART) algorithms

from setuptools import find_packages
from distutils.core import setup

with open('README.md', 'r') as fin:
    long_description=fin.read()

setup(name='ipart',
        version='2.0.2',
        description='IPART is a Python package for the detection and tracking of atmospheric rivers from gridded IVT data using image-processing techniques.',
        long_description=long_description,
        long_description_content_type='text/markdown',
        author='Guangzhi XU',
        author_email='xugzhi1987@gmail.com',
        url='https://github.com/ihesp/IPART',
        classifiers=[
            'Development Status :: 3 - Alpha',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
            'Natural Language :: English',
            'Operating System :: POSIX :: Linux',
            'Operating System :: MacOS',
            'Operating System :: Microsoft :: Windows',
            'Programming Language :: Python :: 3',
            'Topic :: Scientific/Engineering :: Atmospheric Science'
            ],
        install_requires=[
            "netcdf4",
            "matplotlib==2.2.3",
            "scipy",
            "scikit-image",
            "pandas",
            "basemap==1.2.0",
            "networkx"
        ],
        packages=find_packages(include=['ipart', 'ipart.*']),
        license='GPL-3'
        )
