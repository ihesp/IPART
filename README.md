# Atmospheric River (AR) detection and tracking algorithms



## Dependencies

* Python2.7: may not be fully compatible with Python3 yet.
* CDAT (Climatic Data Analysis Tool): (https://github.com/CDAT/cdat). Use a version compatible with Py2. cdat-lite will do
* numpy (developed in 1.16.5)
* scipy (developed in 1.2.1)
* matplotlib (developed in 2.2.3)
* basemap (developed in 1.2.0)
* pandas (developed in 0.23.4)
* networkx (tested 1.11 and 2.2)
* scikit-image (developed in 0.14.2)



## Inventory

* utils/func.py: general purpose functions.
* utils/plot.py: plotting functions.
* utils/rdp.py: Ramer-Douglas-Peucker algorithm.
* utils/peak_prominence2d.py: Separate local peaks using peak prominence.
* compute_thr_singlefile.py: compute THR on IVT data. Process a single file.
* compute_thr_multifile.py: compute THR on IVT data. Process multiple files.
* river_tracker1.py: detect ARs on individual time steps. Requires outputs from `compute_thr.py`.
* river_tracker2.py: track individual AR records across time. Requires outputs from `river_tracker1.py`.
* river_tracker1_funcs.py: utility functions used in `river_tracker1.py`.


## Contribution

If you encounter problems or would like to help improve the code, please don't
hesitate to fire up an issue or pull request.



