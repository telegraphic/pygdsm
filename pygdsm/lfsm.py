import healpy as hp
import numpy as np
from astropy import units
from scipy.interpolate import interp1d, pchip

from .base_observer import BaseObserver
from .base_skymodel import BaseSkyModel
from .component_data import LFSM_DATA_URL, download_from_url_list

T_CMB = 2.725


def equatorial_to_galactic_coords(nside):
    """
    Precompute the pixel coordinates (in equatorial frame) needed to resample
    a map into Galactic coordinates, for a given nside. This depends only on
    nside, so it can be computed once and reused for every frequency/map.
    """
    rotCG = hp.rotator.Rotator(coord=("C", "G"))
    npix = hp.pixelfunc.nside2npix(nside)
    theta, phi = hp.pixelfunc.pix2ang(nside, np.arange(npix))
    theta_new, phi_new = rotCG(theta, phi, inv=True)
    return theta_new, phi_new


class LowFrequencySkyModel(BaseSkyModel):
    """LWA1 Low Frequency Sky Model"""

    def __init__(self, freq_unit="MHz", include_cmb=False, interpolation="cubic"):
        """Global sky model (GSM) class for generating sky models.

        Parameters
        ----------
        freq_unit (str): Frequency unit to use, defaults to MHz
        include_cmb (bool):  Choose whether to include the CMB. Defaults to False. A value of
                             T_CMB = 2.725 K is used if True.
        interpolation (str): 'cubic' or 'pchip'. Choose whether to use cubic spline
                             interpolation or piecewise cubic hermitian interpolating
                             polynomial (PCHIP) for the PCA component coefficients.
                             PCHIP is designed to never locally overshoot the data,
                             whereas splines are designed to have smooth first and
                             second derivatives; near the sparse, low-frequency edge
                             of the model's tabulated frequencies, cubic splines can
                             overshoot and produce unphysical negative temperatures.
                             Defaults to 'cubic' for backwards compatibility.
        """
        data_unit = "K"
        basemap = "LFSS"

        if interpolation not in ("cubic", "pchip"):
            raise RuntimeError(
                "INTERPOLATION ERROR: %s not supported. Only cubic, pchip are allowed."
                % interpolation
            )

        # download component data as needed using astropy cache
        LFSM_FILEPATH = download_from_url_list(LFSM_DATA_URL)

        super(LowFrequencySkyModel, self).__init__(
            "LFSM", LFSM_FILEPATH, freq_unit, data_unit, basemap
        )

        self.pca_map = self.h5["lfsm_component_maps_3.0deg.dat"][:]
        self.pca_components = self.h5["lfsm_components.dat"][:]
        self.nside = 256
        self._eq2gal_theta, self._eq2gal_phi = equatorial_to_galactic_coords(self.nside)

        self.include_cmb = include_cmb
        self.interpolation_method = interpolation

        freqs = self.pca_components[:, 0]
        sigmas = self.pca_components[:, 1]
        comps = self.pca_components[:, 2:]

        self.scaleFunc = interp1d(np.log(freqs), np.log(sigmas), kind="slinear")

        self.compFuncs = []
        for i in range(comps.shape[1]):
            if self.interpolation_method == "pchip":
                self.compFuncs.append(pchip(np.log(freqs), comps[:, i]))
            else:
                self.compFuncs.append(interp1d(np.log(freqs), comps[:, i], kind="cubic"))

    def generate(self, freqs):
        """Generate a global sky model at a given frequency or frequencies

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
        freqs_mhz = freqs.to("MHz").value

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

            map_out[ff] = hp.get_interp_val(
                map_out[ff], self._eq2gal_theta, self._eq2gal_phi
            )

        map_out = map_out.squeeze()

        if self.include_cmb == False:
            map_out -= T_CMB

        self.generated_map_data = map_out
        self.generated_map_freqs = freqs
        return map_out


class LFSMObserver(BaseObserver):
    def __init__(self, gsm=None, **kwargs):
        """Initialize the Observer object.

        Parameters
        ----------
        gsm: LowFrequencySkyModel instance, optional
            A pre-instantiated sky model. If not provided, one is created
            using any supplied keyword arguments.
        **kwargs:
            Keyword arguments passed to LowFrequencySkyModel (e.g. freq_unit,
            include_cmb).
        """
        if gsm is None:
            gsm = LowFrequencySkyModel(**kwargs)
        super(LFSMObserver, self).__init__(gsm=gsm)
