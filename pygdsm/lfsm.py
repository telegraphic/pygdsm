import numpy as np
from astropy import units
import healpy as hp
from scipy.interpolate import interp1d

from .component_data import LFSM_FILEPATH
from .base_skymodel import BaseSkyModel
from .base_observer import BaseObserver

T_CMB = 2.725

def rotate_equatorial_to_galactic(map):
    """
    Given a map in equatorial coordinates, convert it to Galactic coordinates.
    """
    rotCG = hp.rotator.Rotator(coord=('C', 'G'))
    nSides = hp.pixelfunc.npix2nside(map.size)
    theta, phi = hp.pixelfunc.pix2ang(nSides, range(map.size))
    theta_new, phi_new = rotCG(theta, phi, inv=True)
    map2 = hp.get_interp_val(map, theta_new, phi_new)
    return map2


class LowFrequencySkyModel(BaseSkyModel):
    """ LWA1 Low Frequency Sky Model
    """

    def __init__(self,  freq_unit='MHz', include_cmb=False):
        """ Global sky model (GSM) class for generating sky models.

        Parameters
        ----------
        freq_unit (str): Frequency unit to use, defaults to MHz
        include_cmb (bool):  Choose whether to include the CMB. Defaults to False. A value of
                             T_CMB = 2.725 K is used if True.
        """
        data_unit = 'K'
        basemap = 'LFSS'
        super(LowFrequencySkyModel, self).__init__('LFSM', LFSM_FILEPATH, freq_unit, data_unit, basemap)

        self.pca_map = self.h5['lfsm_component_maps_3.0deg.dat'][:]
        self.pca_components =  self.h5['lfsm_components.dat'][:]
        self.nside = 256

        self.include_cmb = include_cmb

        freqs  = self.pca_components[:, 0]
        sigmas = self.pca_components[:, 1]
        comps  = self.pca_components[:, 2:]

        self.scaleFunc = interp1d(np.log(freqs), np.log(sigmas), kind='slinear')

        self.compFuncs = []
        for i in range(comps.shape[1]):
            self.compFuncs.append(interp1d(np.log(freqs), comps[:, i], kind='cubic'))

    def generate(self, freqs):
        """ Generate a global sky model at a given frequency or frequencies

        Parameters
        ----------
        freqs: float or np.array
            Frequency for which to return GSM model

        Returns
        -------
        gsm: np.array
            Global sky model in healpix format, with NSIDE=256. Output map
            is in galactic coordinates, and in antenna temperature units (K).

        """
        # convert frequency values into Hz
        freqs = np.array(freqs) * units.Unit(self.freq_unit)
        freqs_mhz = freqs.to('MHz').value

        if isinstance(freqs_mhz, float):
            freqs_mhz = np.array([freqs_mhz])

        try:
            assert np.min(freqs_mhz) >= 10
            assert np.max(freqs_mhz) <= 408
        except AssertionError:
            raise RuntimeError("Frequency values lie outside 10 MHz < f < 408 MHz")

        map_out = 0.0
        if isinstance(freqs, np.ndarray):
            if freqs.ndim > 0:
                map_out = np.zeros(shape=(freqs.shape[0], hp.nside2npix(self.nside)))
            else:
                map_out = np.zeros(shape=(1, hp.nside2npix(self.nside))) 
        else:
            map_out = np.zeros(shape=(1, hp.nside2npix(self.nside)))

        for ff in range(map_out.shape[0]):
            for i, compFunc in enumerate(self.compFuncs):
                map_out[ff] += compFunc(np.log(freqs_mhz[ff])) * self.pca_map[:, i]
            map_out[ff] *= np.exp(self.scaleFunc(np.log(freqs_mhz[ff])))

            map_out[ff] =  rotate_equatorial_to_galactic(map_out[ff])

        map_out = map_out.squeeze()

        if self.include_cmb == False:
            map_out -= T_CMB

        self.generated_map_data = map_out
        self.generated_map_freqs = freqs
        return map_out


class LFSMObserver(BaseObserver):
    def __init__(self):
        """ Initialize the Observer object.

        Calls ephem.Observer.__init__ function and adds on gsm
        """
        super(LFSMObserver, self).__init__(gsm=LowFrequencySkyModel)