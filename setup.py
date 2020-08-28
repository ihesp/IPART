#!/usr/bin/python

# Image-Processing based Atmospheric River Tracking (IPART) algorithms

from setuptools import find_packages, setup
#from distutils.core import setup

setup(name='ipart',
        version='3.0.4',
        description='IPART is a Python package for the detection and tracking of atmospheric rivers from gridded IVT data using image-processing techniques.',
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
            "numpy",
            "scipy",
            "netcdf4",
            "pandas",
            "scikit-image",
            "networkx",
            "matplotlib",
            "cartopy"
        ],
        packages=find_packages(include=['ipart', 'ipart.*']),
        license='GPL-3.0-or-later'
        )
