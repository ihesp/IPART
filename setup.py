#!/usr/bin/python

# Image-Processing based Atmospheric River Tracking (IPART) algorithms

from distutils.core import setup
from setuptools import find_packages

setup(name='ipart',
        version='2.0',
        description='IPART is a Python package for the detection and tracking of atmospheric rivers from gridded IVT data using image-processing techniques.',
        author='Guangzhi XU',
        author_email='xugzhi1987@gmail.com',
        url='https://github.com/ihesp/IPART',
        install_requires=[
            "cdms2",
            "matplotlib",
            "scipy",
            "scikit-image",
            "pandas",
            "cdutil",
            "basemap",
            "networkx",
        ],
        packages=find_packages(include=['ipart', 'ipart.*']),
        license='GPL-3'
        )
