[![ascl:1603.013](https://img.shields.io/badge/ascl-1603.013-blue.svg)](https://ascl.net/1603.013)
[![codecov](https://codecov.io/gh/telegraphic/pygdsm/branch/master/graph/badge.svg)](https://codecov.io/gh/telegraphic/pygdsm)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/) 
 
PyGDSM
=====

![skymodels.jpg](https://github.com/telegraphic/pygdsm/raw/master/docs/skymodels.jpg)



`PyGDSM` is a Python interface for global diffuse sky models: all-sky maps in Healpix format of diffuse Galactic radio emission.

This package includes interfaces to:
 * **GSM2008:** A model of diffuse Galactic radio emission from 10 MHz to 100 GHz, [Oliveira-Costa et. al., (2008)](https://ui.adsabs.harvard.edu/abs/2008MNRAS.388..247D/abstract). 
 * **GSM2016:** An improved model of diffuse galactic radio emission from 10 MHz to 5 THz, [Zheng et. al., (2016)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.464.3486Z/abstract).
 * **LFSS:** The LWA1 Low Frequency Sky Survey (10-408 MHz) [Dowell et. al. (2017)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.469.4537D/abstract). 
 * **Haslam:** A frequency-scaled model using a spectral index, based on the [Haslam 408 MHz](https://lambda.gsfc.nasa.gov/product/foreground/fg_2014_haslam_408_info.cfm) all-sky map.

In general, these are *not* wrappers of the original code (GSM2008 was written in Fortan and GSM2016 in C); instead they provides a uniform API with some additional features and advantages, such as healpy integration for imaging, and sky rotation for observed skies. 


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

        git clone https://github.com/telegraphic/pygdsm
       
An run `pip install .`. On first run, the sky model data will be downloaded from Zenodo / LAMBDA. These are about 500 MB total, and will be downloaded into your astropy cache (`~/.astropy/`). The data are hosted on [Zenodo](https://zenodo.org/record/3479985#.XaASx79S-AY).

Examples
---------

To get a quick feel of what `PyGDSM` does, have a look at the 
[GSM2008 quickstart guide](http://nbviewer.ipython.org/github/telegraphic/PyGDSM/blob/master/docs/pygdsm_quickstart.ipynb), and the new
[GSM2016 quickstart guide](http://nbviewer.ipython.org/github/telegraphic/PyGDSM/blob/master/docs/pygdsm2016_quickstart.ipynb).

Q & A
-----

**Q. What's the difference between this and the `gsm.f` from the main GSM2008 website?**
     The `gsm.f` is a very basic Fortran code, which reads and writes values to and from
     ASCII files, and uses a command line interface for input. If you want to run this code
     on an ancient computer with nothing but Fortran installed, then `gsm.f` is the way to go. 
     In contrast, `PyGDSM` is a Python code that leverages a lot of other Packages so that you 
     can do more stuff more efficiently. For example: you can view a sky model in a healpy 
     image; you can write a sky model to a Healpix FITS file; and believe it or not, the 
     Python implementation is *much faster*. Have a look at the 
     [quickstart guide](http://nbviewer.ipython.org/github/telegraphic/PyGDSM/blob/master/docs/pygdsm_quickstart.ipynb)
     to get a feel for what `PyGDSM` does.

**Q. Are the outputs of `gsm.f` and `pygdsm` identical?**.  
     **NO**. The cubic  spline interpolation implementation differs, so values will differ by as 
     much as a few percent. The interpolation code used in `gsm.f` does not have an open-source
     license (it's from [Numerical Recipes](http://www.nr.com/licenses/) ), so we haven't 
     implemented it (one could probably come up with an equivalent that didn't infringe).
     Nevertheless, the underlying PCA data are identical, and I've run tests to check that
     the two outputs are indeed comparable. 

**Q. What's the difference between this and the [Zheng et. al. github repo](https://github.com/jeffzhen/gsm2016)?**
     `pygdsm` provides two classes: `GlobalSkyModel16()` and `GSMObserver16()`, which once instantiated
     provide methods for programatically generating sky models. The Zheng et. al. github repo is a 
     simple, low-dependency, command line tool. As of PyGDSM 1.4.0, we have implemented improved interpolation
     via cubic spline or PCHIP, which avoids discontinuities identified using the 2016 method. Have a look at the 
     [GSM2016 quickstart guide](http://nbviewer.ipython.org/github/telegraphic/PyGDSM/blob/master/docs/pygdsm2016_quickstart.ipynb)
     to get a feel for what `PyGDSM` does.

**Q. Why does this package download so much data when first run?**
     The package size is dominated by the PCA healpix maps, which have about 3 million points each.
     They're compressed using HDF5 LZF, so are actually about 3x smaller than the `*.dat`
     files that come in the original `gsm.tar.gz` file. The next biggest thing is test data,
     so that the output can be compared against precomputed output from `gsm.f`. The package now also includes
     the Zheng et. al. data, which is another ~300 MB.
   

References
----------

The sky model data contained here is from:
* GSM2008 http://space.mit.edu/~angelica/gsm/index.html (link no longer active)
* GSM2016 https://github.com/jeffzhen/gsm2016
* LFSS https://lda10g.alliance.unm.edu/LWA1LowFrequencySkySurvey/
* Haslam https://lambda.gsfc.nasa.gov/product/foreground/fg_2014_haslam_408_info.cfm


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
https://ui.adsabs.harvard.edu/abs/2017MNRAS.469.4537D/abstract

An improved source-subtracted and destriped 408-MHz all-sky map 
M. Remazeilles, C. Dickinson,A.J. Banday,  M. Bigot-Sazy, T. Ghosh
MNRAS 451, 4, 4311-4327 (2014)
https://ui.adsabs.harvard.edu/abs/2015MNRAS.451.4311R/abstract
 
```

PyGSDM has an [ascl.net entry](https://ascl.net/1603.013):

```
D. C. Price, 2016, 2.0.0, Astrophysics Source Code Library, 1603.013
```

License
-------

All *code* in PyGDSM is licensed under the MIT license (not the underlying *data*). 
The PCA data, by Zheng et. al. is licensed under MIT also (see https://github.com/jeffzhen/gsm2016).

