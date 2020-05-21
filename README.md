[![Build Status](https://travis-ci.org/telegraphic/pygdsm.svg?branch=master)](https://travis-ci.org/telegraphic/pygdsm)
[![codecov](https://codecov.io/gh/telegraphic/pygdsm/branch/master/graph/badge.svg)](https://codecov.io/gh/telegraphic/pygdsm)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/) 
 
PyGDSM
=====

**Note: LFSS methods are in beta!!**

`PyGDSM` is a Python interface for global diffuse sky models: all-sky maps in Healpix format of diffuse Galactic radio emission.

This package includes interfaces to:
 * **GSM2008:** A model of diffuse Galactic radio emission from 10 MHz to 100 GHz, [Oliveira-Costa et. al., (2008)](https://ui.adsabs.harvard.edu/abs/2008MNRAS.388..247D/abstract).
 * **GSM2016:** An improved model of diffuse galactic radio emission from 10 MHz to 5 THz, [Zheng et. al., (2016)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.464.3486Z/abstract).
 * **LFSS:** The LWA1 Low Frequency Sky Survey (10-408 MHz) [Dowell et. al. (2017)](https://ui.adsabs.harvard.edu/abs/2008MNRAS.388..247D/abstracthttp://arxiv.org/abs/1605.04920).

This is *not* a wrapper of the original code, it is a python-based equivalent that provides a uniform API with some additional features and advantages, such as healpy integration for imaging, and sky rotation for observed skies. 


Quickstart
----------

The first thing to do will be to make sure you've got the dependencies: 

* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org/install.html)
* [healpy](http://healpy.readthedocs.org/en/latest/)
* [h5py](http://www.h5py.org/)
* [astropy](http://www.astropy.org/)

Then you should be able to install with:

        pip install git+https://github.com/telegraphic/pygdsm

Alternatively, clone the directory:

        git clone https://github.com/telegraphic/PyGDSM
       
An run `pip install .`. On first run, the sky model data will be downloaded from Zenodo. These are about 500 MB total, and will be downloaded into your astropy cache (`~/.astropy/`). The data are hosted on [Zenodo](https://zenodo.org/record/3479985#.XaASx79S-AY).

Examples
---------

To get a quick feel of what `PyGDSM` does, have a look at the 
[GSM2008 quickstart guide](http://nbviewer.ipython.org/github/telegraphic/PyGSM/blob/master/docs/pygsm_quickstart.ipynb), and the new
[GSM2016 quickstart guide](http://nbviewer.ipython.org/github/telegraphic/PyGSM/blob/master/docs/pygsm2016_quickstart.ipynb).

Q & A
-----

**Q. What's the difference between this and the `gsm.f` from the main GSM2008 website?**
     The `gsm.f` is a very basic Fortran code, which reads and writes values to and from
     ASCII files, and uses a command line interface for input. If you want to run this code
     on an ancient computer with nothing by Fortran installed, then `gsm.f` is the way to go. 
     In contrast, `PyGSM` is a Python code that leverages a lot of other Packages so that you 
     can do more stuff more efficiently. For example: you can view a sky model in a healpy 
     image; you can write a sky model to a Healpix FITS file; and believe it or not, the 
     Python implementation is *much faster*. Have a look at the 
     [quickstart guide](http://nbviewer.ipython.org/github/telegraphic/PyGSM/blob/master/docs/pygsm_quickstart.ipynb)
     to get a feel for what `PyGSM` does.

**Q. Are the outputs of `gsm.f` and `pygsm` identical?** At the moment: **no**. The cubic
     spline interpolation implementation differs, so values will differ by as much as 
     a few percent. The interpolation code used in `gsm.f` does not have an open-source
     license (it's from [Numerical Recipes](http://www.nr.com/licenses/) ), so we haven't 
     implemented it (one could probably come up with an equivalent that didn't infringe).
     Nevertheless, the underlying PCA data are identical, and I've run tests to check that
     the two outputs are indeed comparable. 

**Q. What's the difference between this and the Zheng et. al. github repo?**
     `pygsm` provides two classes: `GlobalSkyModel2016()` and `GSMObserver2016()`, which once instantiated
     provide methods for programatically generating sky models. The Zheng et. al. github repo is a 
     simple, low-dependency, command line tool. Have a look at the 
     [GSM2016 quickstart guide](http://nbviewer.ipython.org/github/telegraphic/PyGSM/blob/master/docs/pygsm2016_quickstart.ipynb)
     to get a feel for what `PyGSM` does.

**Q. Why is this package so large?**
     The package size is dominated by the PCA healpix maps, which have about 3 million points each.
     They're compressed using HDF5 LZF, so are actually about 3x smaller than the `*.dat`
     files that come in the original `gsm.tar.gz` file. The next biggest thing is test data,
     so that the output can be compared against precomputed output from `gsm.f`. The package now also includes
     the Zheng et. al. data, which is another ~300 MB.

**Q. Why do I need h5py?**
     `h5py` is required to read the PCA data, which are stored in a HDF5 file. Reading from
     HDF5 into Python is incredibly efficient, and the compression is transparent to the end user.
     This means that you can't eyeball the data using `vim` or `less` or a text editor, but if
     you're trying to do that on a file with millions of data points you're doing science wrong anyway.
   

References
----------

The PCA data contained here is from:
* GSM2008 http://space.mit.edu/~angelica/gsm/index.html 
* GSM2016 https://github.com/jeffzhen/gsm2016
* LFSS https://lda10g.alliance.unm.edu/LWA1LowFrequencySkySurvey/


```
A model of diffuse Galactic radio emission from 10 MHz to 100 GHz
A. de Oliveira-Costa, M. Tegmark, B.M. Gaensler, J. Jonas, T.L. Landecker and P. Reich
MNRAS 388, 247-260 (2008)
https://ui.adsabs.harvard.edu/abs/2008MNRAS.388..247D/abstract

An Improved Model of Diffuse Galactic Radio Emission from 10 MHz to 5 THz
H. Zheng, M. Tegmark, J. Dillon, A. Liu, A. Neben, J. Jonas, P. Reich, W.Reich
MNRAS, 464, 3, 3486-3497 (2017)
https://ui.adsabs.harvard.edu/abs/2017MNRAS.464.3486Z/abstract

The LWA1 Low Frequency Sky Survey
J. Dowell, G. B. Taylor, F. Schinzel, N. E. Kassim, K. Stovall
MNRAS, 469, 4, 4537-4550 (2017)
```

PyGSDM has an [ascl.net entry](https://ascl.net/1603.013):

```
D. C. Price, 2016, 2.0.0, Astrophysics Source Code Library, 1603.013
```

License
-------

All *code* in PyGSM is licensed under the MIT license (not the underlying *data*). 
The PCA data, by Zheng et. al. is licensed under MIT also (see https://github.com/jeffzhen/gsm2016).

