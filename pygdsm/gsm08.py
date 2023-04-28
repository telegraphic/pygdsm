"""
gsm.py
======

Python interface for the Global Sky Model (GSM) or Oliveira-Costa et. al.
This is a python-based equivalent to the Fortran `gsm.f` that comes with the
original data. Instead of the original ASCII DAT files that contain the PCA
data, data are stored in HDF5, which is more efficient.

References
----------
A. de Oliveira-Costa, M. Tegmark, B.M. Gaensler, J. Jonas, T.L. Landecker and P. Reich
A model of diffuse Galactic radio emission from 10 MHz to 100 GHz
Mon. Not. R. Astron. Soc. 388, 247-260 (2008)
doi:10.111/j.1365-2966.2008.13376.x

PCA data from: space.mit.edu/home/angelica/gsm

"""

import numpy as np
from scipy.interpolate import interp1d, pchip
import h5py
from astropy import units
import healpy as hp


from .component_data import GSM_FILEPATH
from .plot_utils import show_plt
from .base_observer import BaseObserver
from .base_skymodel import BaseSkyModel

T_CMB = 2.725

class GlobalSkyModel(BaseSkyModel):
    """ Global sky model (GSM) class for generating sky models.
    """

    def __init__(self, freq_unit='MHz', basemap='haslam', interpolation='pchip', include_cmb=False):
        """ Global sky model (GSM) class for generating sky models.

        Upon initialization, the map PCA data are loaded into memory and interpolation
        functions are pre-computed.

        Parameters
        ----------
        freq_unit: 'Hz', 'MHz', or 'GHz'
            Unit of frequency. Defaults to 'MHz'.
        basemap: 'haslam', 'wmap' or '5deg'
            GSM version to generate. The 5deg map has 5.1 degree resolution.
            This is a synthesized map made of all the maps considered in
            the de Oliveira-Costa et. al. paper
            At frequencies below 1GHz, haslam is preferred; you get higher
            resolution (1 degree) by locking to the Haslam 408 MHz map.
            At CMB frequencies, it is best to lock to the WMAP 23 GHz
            map, which is presented denoised with 2 degree resolution.
        interpolation: 'cubic' or 'pchip'
            Choose whether to use cubic spline interpolation or
            piecewise cubic hermitian interpolating polynomial (PCHIP).
            PCHIP is designed to never locally overshoot data, whereas
            splines are designed to have smooth first and second derivatives.
        include_cmb: bool (default False)
            Choose whether to include the CMB. Defaults to False. A value of
            T_CMB = 2.725 K is used.

        Notes
        -----
        The scipy `interp1d` function does not allow one to explicitly
        set second derivatives to zero at the endpoints, as is done in
        the original GSM. As such, results will differ. Further, we default
        to use PCHIP interpolation.

        """

        try:
            assert basemap in {'5deg', 'wmap', 'haslam'}
        except AssertionError:
            raise RuntimeError("GSM basemap unknown: %s. Choose '5deg', 'haslam' or 'wmap'" % basemap)

        try:
            assert interpolation in {'cubic', 'pchip'}
        except AssertionError:
            raise RuntimeError("Interpolation must be set to either 'cubic' or 'pchip'")

        data_unit = 'K'
        super(GlobalSkyModel, self).__init__('GSM2008', GSM_FILEPATH, freq_unit, data_unit, basemap)

        self.interpolation_method = interpolation
        self.include_cmb = include_cmb


        self.pca_map_data = None
        self.interp_comps = None
        self.update_interpolants()

        self.generated_map_data = None
        self.generated_map_freqs = None
        self.nside = 512

    def update_interpolants(self):
       # Choose the PCA map to load from the HDF5 file
        pca_map_dict = {"5deg": "component_maps_5deg",
                        "haslam": "component_maps_408locked",
                        "wmap": "component_maps_23klocked"}
        pca_map_key = pca_map_dict[self.basemap]
        self.pca_map_data = self.h5[pca_map_key][:]

        # Now, load the PCA eigenvalues
        pca_table = self.h5["components"][:]
        pca_freqs_mhz = pca_table[:, 0]
        pca_scaling   = pca_table[:, 1]
        pca_comps     = pca_table[:, 2:].T

        # Interpolate to the desired frequency values
        ln_pca_freqs = np.log(pca_freqs_mhz)

        if self.interpolation_method == 'cubic':
            spl_scaling = interp1d(ln_pca_freqs, np.log(pca_scaling), kind='cubic')
            spl1 = interp1d(ln_pca_freqs,   pca_comps[0],   kind='cubic')
            spl2 = interp1d(ln_pca_freqs,   pca_comps[1],   kind='cubic')
            spl3 = interp1d(ln_pca_freqs,   pca_comps[2],   kind='cubic')

        else:
            spl_scaling = pchip(ln_pca_freqs, np.log(pca_scaling))
            spl1 = pchip(ln_pca_freqs,   pca_comps[0])
            spl2 = pchip(ln_pca_freqs,   pca_comps[1])
            spl3 = pchip(ln_pca_freqs,   pca_comps[2])
        self.interp_comps = (spl_scaling, spl1, spl2, spl3)

    def generate(self, freqs):
        """ Generate a global sky model at a given frequency or frequencies

        Parameters
        ----------
        freqs: float or np.array
            Frequency for which to return GSM model

        Returns
        -------
        gsm: np.array
            Global sky model in healpix format, with NSIDE=512. Output map
            is in galactic coordinates, and in antenna temperature units (K).

        """
        # convert frequency values into Hz
        freqs = np.array(freqs) * units.Unit(self.freq_unit)
        freqs_mhz = freqs.to('MHz').value

        if isinstance(freqs_mhz, float):
            freqs_mhz = np.array([freqs_mhz])

        try:
            assert np.min(freqs_mhz) >= 10
            assert np.max(freqs_mhz) <= 94000
        except AssertionError:
            raise RuntimeError("Frequency values lie outside 10 MHz < f < 94 GHz")

        # Load interpolators and do interpolation
        ln_freqs     = np.log(freqs_mhz)
        spl_scaling, spl1, spl2, spl3 = self.interp_comps
        comps = np.row_stack((spl1(ln_freqs), spl2(ln_freqs), spl3(ln_freqs)))
        scaling = np.exp(spl_scaling(ln_freqs))

        # Finally, compute the dot product via einsum (awesome function)
        # c=comp, f=freq, p=pixel. We want to dot product over c for each freq
        #print comps.shape, self.pca_map_data.shape, scaling.shape
        map_out = np.einsum('cf,pc,f->fp', comps, self.pca_map_data, scaling)
        
        if self.include_cmb:
            map_out += T_CMB

        if map_out.shape[0] == 1:
            map_out = map_out[0]
        self.generated_map_data = map_out
        self.generated_map_freqs = freqs



        return map_out

    def set_basemap(self, new_basemap):
        self.basemap = new_basemap
        self.update_interpolants()
        if self.generated_map_freqs is not None:
            self.generate(self.generated_map_freqs)

    def set_freq_unit(self, new_unit):
        self.freq_unit = new_unit
        self.update_interpolants()
        if self.generated_map_freqs is not None:
            self.generate(self.generated_map_freqs)

    def set_interpolation_method(self, new_method):
        self.interpolation_method = new_method
        self.update_interpolants()
        if self.generated_map_freqs is not None:
            self.generate(self.generated_map_freqs)


class GSMObserver(BaseObserver):
    def __init__(self):
        """ Initialize the Observer object.

        Calls ephem.Observer.__init__ function and adds on gsm
        """
        super(GSMObserver, self).__init__(gsm=GlobalSkyModel)