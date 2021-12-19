import numpy as np
import h5py
from astropy import units
import healpy as hp
import ephem

from .component_data import GSM2016_FILEPATH
from .plot_utils import show_plt
from .base_observer import BaseObserver
from .base_skymodel import BaseSkyModel

kB = 1.38065e-23
C = 2.99792e8
h = 6.62607e-34
T = 2.725
hoverk = h / kB


def K_CMB2MJysr(K_CMB, nu):#in Kelvin and Hz
    B_nu = 2 * (h * nu)* (nu / C)**2 / (np.exp(hoverk * nu / T) - 1)
    conversion_factor = (B_nu * C / nu / T)**2 / 2 * np.exp(hoverk * nu / T) / kB
    return  K_CMB * conversion_factor * 1e20#1e-26 for Jy and 1e6 for MJy

def K_RJ2MJysr(K_RJ, nu):#in Kelvin and Hz
    conversion_factor = 2 * (nu / C)**2 * kB
    return  K_RJ * conversion_factor * 1e20#1e-26 for Jy and 1e6 for MJy


def rotate_map(hmap, rot_theta, rot_phi, nest=True):
    nside = hp.npix2nside(len(hmap))

    # Get theta, phi for non-rotated map
    t, p = hp.pix2ang(nside, np.arange(hp.nside2npix(nside)), nest= nest)  # theta, phi

    # Define a rotator
    r = hp.Rotator(deg=False, rot=[rot_phi, rot_theta])

    # Get theta, phi under rotated co-ordinates
    trot, prot = r(t, p)

    # Inerpolate map onto these co-ordinates
    rot_map = hp.get_interp_val(hmap, trot, prot, nest= nest)

    return rot_map


class GlobalSkyModel2016(BaseSkyModel):
    """ Global sky model (GSM) class for generating sky models.
    """

    def __init__(self, freq_unit='MHz', data_unit='TCMB', resolution='hi', theta_rot=0, phi_rot=0):
        """ Global sky model (GSM) class for generating sky models.

        Upon initialization, the map PCA data are loaded into memory and interpolation
        functions are pre-computed.

        Parameters
        ----------
        freq_unit: 'Hz', 'MHz', or 'GHz'
            Unit of frequency. Defaults to 'MHz'.
        data_unit: 'MJysr', 'TCMB', 'TRJ'
            Unit of output data. MJy/Steradian, T_CMB in Kelvin, or T_RJ.
        resolution: 'hi' or 'low'
            Resolution of output map. Either 300 arcmin (low) or 24 arcmin (hi).
            For frequencies under 10 GHz, output is 48 arcmin.

        Notes
        -----

        """

        if data_unit not in ['MJysr', 'TCMB', 'TRJ']:
            raise RuntimeError("UNIT ERROR: %s not supported. Only MJysr, TCMB, TRJ are allowed." % data_unit)

        if resolution.lower() in ('hi', 'high', 'h'):
            resolution = 'hi'
        elif resolution.lower() in ('low', 'lo', 'l'):
            resolution = 'low'
        else:
            raise RuntimeError("RESOLUTION ERROR: Must be either hi or low, not %s" % resolution)

        super(GlobalSkyModel2016, self).__init__('GSM2016', GSM2016_FILEPATH, freq_unit, data_unit, basemap='')

        self.resolution = resolution

        # Map data to load
        labels = ['Synchrotron', 'CMB', 'HI', 'Dust1', 'Dust2', 'Free-Free']
        self.n_comp = len(labels)

        if resolution=='hi':
            self.nside = 1024
            self.map_ni = np.array([self.h5['highres_%s_map'%lb][:] for lb in labels])
        else:
            self.nside = 64
            self.map_ni = np.array(self.h5['lowres_maps'])

        self.spec_nf = self.h5['spectra'][:]

        if theta_rot or phi_rot:
            for i,map in enumerate(self.map_ni):
                self.map_ni[i] = rotate_map(map, theta_rot, phi_rot, nest=True)

    def generate(self, freqs):
        """ Generate a global sky model at a given frequency or frequencies

        Parameters
        ----------
        freqs: float or np.array
            Frequency for which to return GSM model

        Returns
        -------
        gsm: np.array
            Global sky model in healpix format, with NSIDE=1024. Output map
            is in galactic coordinates, ring format.

        """

        # convert frequency values into Hz
        freqs = np.array(freqs) * units.Unit(self.freq_unit)
        freqs_ghz = freqs.to('GHz').value

        if isinstance(freqs_ghz, float):
            freqs_ghz = np.array([freqs_ghz])

        try:
            assert np.min(freqs_ghz) >= 0.01
            assert np.max(freqs_ghz) <= 5000
        except AssertionError:
            raise RuntimeError("Frequency values lie outside 10 MHz < f < 5 THz: %s")

        map_ni = self.map_ni
        # if self.resolution == 'hi':
        #     map_ni = self.map_ni_hr
        # else:
        #     map_ni = self.map_ni_lr

        spec_nf = self.spec_nf
        nfreq = spec_nf.shape[1]

        output = np.zeros((len(freqs_ghz), map_ni.shape[1]), dtype='float32')
        for ifreq, freq in enumerate(freqs_ghz):

            left_index = -1
            for i in range(nfreq - 1):
                if spec_nf[0, i] <= freq <= spec_nf[0, i + 1]:
                    left_index = i
                    break

            # Do the interpolation
            interp_spec_nf = np.copy(spec_nf)
            interp_spec_nf[0:2] = np.log10(interp_spec_nf[0:2])
            x1 = interp_spec_nf[0, left_index]
            x2 = interp_spec_nf[0, left_index + 1]
            y1 = interp_spec_nf[1:, left_index]
            y2 = interp_spec_nf[1:, left_index + 1]
            x = np.log10(freq)
            interpolated_vals = (x * (y2 - y1) + x2 * y1 - x1 * y2) / (x2 - x1)
            output[ifreq] = np.sum(10.**interpolated_vals[0] * (interpolated_vals[1:, None] * map_ni), axis=0)

            output[ifreq] = hp.pixelfunc.reorder(output[ifreq], n2r=True)

            if self.data_unit == 'TCMB':
                conversion = 1. / K_CMB2MJysr(1., 1e9 * freq)
            elif self.data_unit == 'TRJ':
                conversion = 1. / K_RJ2MJysr(1., 1e9 * freq)
            else:
                conversion = 1.
            output[ifreq] *= conversion

#            output.append(result)

        if len(output) == 1:
            output = output[0]
        #else:
        #    map_data = np.row_stack(output)

        self.generated_map_freqs = freqs
        self.generated_map_data = output

        return output


class GSMObserver2016(BaseObserver):
    def __init__(self):
        """ Initialize the Observer object.

        Calls ephem.Observer.__init__ function and adds on gsm
        """
        super(GSMObserver2016, self).__init__(gsm=GlobalSkyModel2016)

