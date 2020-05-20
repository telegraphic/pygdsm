[![Build Status](https://travis-ci.org/telegraphic/pygdsm.svg?branch=master)](https://travis-ci.org/telegraphic/pygdsm)
[![codecov](https://codecov.io/gh/telegraphic/pygdsm/branch/master/graph/badge.svg)](https://codecov.io/gh/telegraphic/pygdsm)
[![astropy](http://img.shields.io/badge/powered%20by-AstroPy-orange.svg?style=flat)](http://www.astropy.org/) 
 
PyGDSM
=====

**Note: LFSS methods are in beta!!**

`PyGDSM` is a Python interface for global diffuse sky models: all-sky maps in Healpix format of diffuse Galactic radio emission.

This package includes interfaces to:
 * The Global Sky Model (GSM2008) of [Oliveira-Costa et. al., (2008)](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2966.2008.13376.x/abstract), 
 * [Zheng et. al., (2016)](http://arxiv.org/abs/1605.04920), GSM2016. 
 * [The LWA1 Low Frequency Sky Model](https://lda10g.alliance.unm.edu/LWA1LowFrequencySkySurvey/), LFSS.

from 10 MHz to 94 GHz (GSM2008), 10 MHz to 5 THz (GSM2016), and 10 MHz to 408 MHz (LFSS.)

This is *not* a wrapper of the original code, it is a python-based equivalent
that provides a useful API which has some additional features and advantages, such as healpy integration for imaging. 

The GSM2008, GSM2016 and LFSS classes provided here are designed to have the same API (i.e. function names & usage).
Instead of the original ASCII DAT files that contain the principal component analysis
(PCA), from the GSM2008, data are stored in HDF5, which can be quickly read into memory, and takes up less space to boot.
Similarly, the GSM2016 and LFSS data is converted into HDF5. 

Quickstart
----------

The first thing to do will be to make sure you've got the dependencies: 

* [numpy](http://www.numpy.org/)
* [scipy](http://www.scipy.org/install.html)
* [healpy](http://healpy.readthedocs.org/en/latest/)
* [h5py](http://www.h5py.org/)
* [astropy](http://www.astropy.org/)

Then clone the directory

        git clone https://github.com/telegraphic/PyGSM
        
You need to have git-lfs installed to download the data files. If you don't, you can get these via wget:

      wget -O gsm2016_components.h5 https://zenodo.org/record/3479985/files/gsm2016_components.h5?download=1
      wget -O gsm_components.h5 https://zenodo.org/record/3479985/files/gsm_components.h5?download=1

Which are hosted on [Zenodo](https://zenodo.org/record/3479985#.XaASx79S-AY).

You may then install this by running `python setup.py install`.

To get a quick feel of what `PyGSM` does, have a look at the 
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

The PCA data contained here is from http://space.mit.edu/~angelica/gsm/index.html and
https://github.com/jeffzhen/gsm2016.

The original GSM2008 paper is:

```
A. de Oliveira-Costa, M. Tegmark, B.M. Gaensler, J. Jonas, T.L. Landecker and P. Reich
A model of diffuse Galactic radio emission from 10 MHz to 100 GHz
Mon. Not. R. Astron. Soc. 388, 247-260 (2008)
doi:10.111/j.1365-2966.2008.13376.x
```

Which is published in [MNRAS](http://onlinelibrary.wiley.com/doi/10.1111/j.1365-2966.2008.13376.x/abstract)
and is also available on the [arXiv](http://arxiv.org/abs/0802.1525).

And the GSM2016 paper is:

```
H. Zheng (MIT), M. Tegmark, J. Dillon, A. Liu, A. Neben, J. Jonas, P. Reich, W.Reich
An Improved Model of Diffuse Galactic Radio Emission from 10 MHz to 5 THz
```

which is available on the [arXiv](http://arxiv.org/abs/1605.04920).

PyGSM has an [ascl.net entry](https://ascl.net/1603.013):

```
D. C. Price, 2016, 2.0.0, Astrophysics Source Code Library, 1603.013
```

License
-------

All *code* in PyGSM is licensed under the MIT license (not the underlying *data*). The PCA data, by Zheng et. al. is licensed under MIT also (see https://github.com/jeffzhen/gsm2016).
