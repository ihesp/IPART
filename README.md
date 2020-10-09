# Image-Processing based Atmospheric River Tracking (IPART) algorithms


## Introduction

IPART (Image-Processing based Atmospheric River Tracking) is a Python package
for automated Atmospheric River (AR) detection, axis finding and AR tracking
from gridded Integrated Vapor Transport (IVT) data, for instance Reanalysis
datasets, or model simulations.

IPART is intended for researchers and students who are interested in the
field of atmospheric river studies in the present day climate or future
projections. Unlike the convectional detection methods that rely on magnitude
thresholding on the intensities of atmospheric vapor fluxes, IPART tackles the
detection task from a spatio-temporal scale perspective and is thus
free from magnitude thresholds.

## Documentation

Further documentation can be found at [https://ipart.readthedocs.io/en/latest/](https://ipart.readthedocs.io/en/latest/).


## Example use case


| ![fig3](joss/fig3.png) |
| :--: |
|*(a) The IVT field in kg/m/s at 1984-01-26 00:00 UTC over the North Hemisphere. (b) the IVT reconstruction field (IVT_rec) at the same time point. (c) the IVT anomaly field (IVT_ano) from the THR process at the same time point.*|

| ![](joss/ar_track_198424.png) |
| :--: |
|*Locations of a track labelled "198424" found in year 1984. Black to yellow color scheme indicates the evolution.*|



## Dependencies

* Python2.7 or Python3.7.
* netCDF4 (tested 1.4.2, 1.5.3 in py2, tested 1.5.3 in py3)
* numpy (developed in 1.16.5 in py2, tested 1.18.1, 1.19.0 in py3)
* scipy (developed in 1.2.1 in py2, tested 1.4.1, 1.5.1 in py3)
* matplotlib (tested 2.2.5 in py2, tested 3.3.1 in py3)
* pandas (developed in 0.23.4, 0.24.2 in py2, tested 1.0.3, 1.0.5 in py3)
* networkx (developed in 1.11 and 2.2 in py2, tested 2.4 in py3)
* scikit-image (developed in 0.14.2, 0.14.3 in py2, tested 0.16.2, 0.17.2 in py3)
* cartopy (optional, only used for plotting. Tested 0.17.0 in py2, tested 1.18.0 in py3)
* OS: Linux or Mac, may work in Windows.

## Installation

Recommend building the Python environment using [Anaconda](https://www.anaconda.com/distribution/).


### Install from conda-forge

In your working Python environment:

```
conda install -c conda-forge ipart
```

will install `ipart` and its dependencies for Python 3.


### Create conda environment using environment file

This way will install the optional `cartopy` package and allow you to run
the notebook examples.

After Anaconda installation, git clone this repository:

```
git clone https://github.com/ihesp/IPART
```

Then build a new conda environment using the environment file provided. For example:

```
cd IPART
conda env create -f environment_py3.yml
```

This creates a new environment named `ipartpy3`. Activate the environment using

```
conda activate ipartpy3
```

After that, you can check the list of packages installed by

```
conda list
```

Similarly for Python 2.7, use

```
conda env create -f environment_py2.yml
```

Finally install IPART using:

```
pip install -e .
```


## tests

To validate installation, issue a new Python session and run

```
import ipart
```

If nothing prints out, installation is successful.

The `tests` folder also contains a number of `unittest`s, to run them (only if you have done a source code install):

```
python -m unittest discover -s tests
```



## Inventory

* docs: readthedocs documentation.
* ipart: core module functions.
* notebooks: a series of jupyter notebooks illustrating the major functionalities of the package.
* scripts: example computation scripts. Can be used as templates to quickly develop your own working scripts.


## Changelog

### v3.0

Make algorithms zonally cyclic.

### v2.0

* restructure into a module `ipart`, separate module from scripts.
* add a `findARsGen()` generator function to yield results at each time point separately.

### v1.0

* initial upload. Can perform AR detection and tracing through time.



## Contributing

Following the guidelines by the [Neurohackademy 2020 curriculum](https://github.com/neurohackademy/nh2020-curriculum), we welcome
contributions from the community. Please create a fork of the project on GitHub
and use a pull request to propose your changes. We strongly encourage creating
an issue before starting to work on major changes, to discuss these changes
first.

## Getting help

Please post issues on the project GitHub page.
