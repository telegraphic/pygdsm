import numpy as np
from scipy.interpolate import interp1d, pchip
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


class GlobalSkyModel16(BaseSkyModel):
    """ Global sky model (GSM) class for generating sky models.
    """

    def __init__(self, freq_unit='MHz', data_unit='TCMB', resolution='hi', theta_rot=0, phi_rot=0, interpolation='pchip', include_cmb=False):
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

        if data_unit not in ['MJysr', 'TCMB', 'TRJ']:
            raise RuntimeError("UNIT ERROR: %s not supported. Only MJysr, TCMB, TRJ are allowed." % data_unit)

        if resolution.lower() in ('hi', 'high', 'h'):
            resolution = 'hi'
        elif resolution.lower() in ('low', 'lo', 'l'):
            resolution = 'low'
        else:
            raise RuntimeError("RESOLUTION ERROR: Must be either hi or low, not %s" % resolution)

        super(GlobalSkyModel16, self).__init__('GSM2016', GSM2016_FILEPATH, freq_unit, data_unit, basemap='')

        self.interpolation_method = interpolation
        self.resolution = resolution
        self.include_cmb = include_cmb

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
        
        # Now borrow code from the orignal GSM2008 model to do a sensible interpolation

        pca_freqs_ghz = spec_nf[0]
        pca_scaling   = spec_nf[1]
        pca_comps     = spec_nf[2:]
         # Interpolate to the desired frequency values
        ln_pca_freqs = np.log(pca_freqs_ghz)
        if self.interpolation_method == 'cubic':
            spl_scaling = interp1d(ln_pca_freqs, np.log(pca_scaling), kind='cubic')
            spl1 = interp1d(ln_pca_freqs,   pca_comps[0],   kind='cubic')
            spl2 = interp1d(ln_pca_freqs,   pca_comps[1],   kind='cubic')
            spl3 = interp1d(ln_pca_freqs,   pca_comps[2],   kind='cubic')
            spl4 = interp1d(ln_pca_freqs,   pca_comps[3],   kind='cubic')
            spl5 = interp1d(ln_pca_freqs,   pca_comps[4],   kind='cubic')
            spl6 = interp1d(ln_pca_freqs,   pca_comps[5],   kind='cubic')

        else:
            spl_scaling = pchip(ln_pca_freqs, np.log(pca_scaling))
            spl1 = pchip(ln_pca_freqs,   pca_comps[0])
            spl2 = pchip(ln_pca_freqs,   pca_comps[1])
            spl3 = pchip(ln_pca_freqs,   pca_comps[2])
            spl4 = pchip(ln_pca_freqs,   pca_comps[3])
            spl5 = pchip(ln_pca_freqs,   pca_comps[4])
            spl6 = pchip(ln_pca_freqs,   pca_comps[5])
            
        self.interp_comps = (spl_scaling, spl1, spl2, spl3, spl4, spl5, spl6)
     
        ln_freqs = np.log(freqs_ghz)
        comps = np.row_stack((spl1(ln_freqs), spl2(ln_freqs), spl3(ln_freqs), spl4(ln_freqs), spl5(ln_freqs), spl6(ln_freqs)))
        scaling = np.exp(spl_scaling(ln_freqs))
        
        # Finally, compute the dot product via einsum (awesome function)
        # c=comp, f=freq, p=pixel. We want to dot product over c for each freq
        #print comps.shape, self.pca_map_data.shape, scaling.shape
        
        output = np.single(np.einsum('cf,pc,f->fp', comps, map_ni.T, scaling))
        

        for ifreq, freq in enumerate(freqs_ghz):
            
            output[ifreq] = hp.pixelfunc.reorder(output[ifreq], n2r=True)

            # DCP 2024.03.29 - Add CMB if requested
            if self.include_cmb:
                output[ifreq] +=  K_CMB2MJysr(T, 1e9 * freq)
                
            if self.data_unit == 'TCMB':
                conversion = 1. / K_CMB2MJysr(1., 1e9 * freq)
            elif self.data_unit == 'TRJ':
                conversion = 1. / K_RJ2MJysr(1., 1e9 * freq)
            else:
                conversion = 1.
            output[ifreq] *= conversion


        if len(output) == 1:
            output = output[0]
        #else:
        #    map_data = np.row_stack(output)

        self.generated_map_freqs = freqs
        self.generated_map_data = output

        return output


class GSMObserver16(BaseObserver):
    def __init__(self):
        """ Initialize the Observer object.

        Calls ephem.Observer.__init__ function and adds on gsm
        """
        super(GSMObserver16, self).__init__(gsm=GlobalSkyModel16)
