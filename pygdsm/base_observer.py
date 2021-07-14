import ephem
import healpy as hp
import numpy as np

from pygdsm.plot_utils import show_plt


class BaseObserver(ephem.Observer):
    """ Observer of the Global Sky Model.

    Generates the Observed sky, for a given point on Earth.
    Applies the necessary rotations and coordinate transformations
    so that the observed 'sky' can be returned, instead of the
    full galaxy-centered GSM.

    This class is based on pyephem's Observer(). The GSM bit can be thought of
    as an 'add on' to ephem.Observer, adding the methods generate() and  view(),
    which allows all-sky images for a given point on earth to be produced.
    """

    def __init__(self, gsm):
        """ Initialize the Observer object.

        Calls ephem.Observer.__init__ function and adds on gsm
        """
        super(BaseObserver, self).__init__()
        self.observed_sky = None
        self.gsm = gsm()
        self._setup()

    def _setup(self):
        self._freq = 100
        # Generate mapping from pix <-> angles
        self.gsm.generate(self._freq)
        self._n_pix  = hp.get_map_size(self.gsm.generated_map_data)
        self._n_side = hp.npix2nside(self._n_pix)
        self._theta, self._phi = hp.pix2ang(self._n_side, np.arange(self._n_pix))

    def generate(self, freq=None):
        """ Generate the observed sky for the observer, based on the GSM.

        Parameters
        ----------
        freq: float
            Frequency of map to generate, in units of MHz (default).

        Returns
        -------
        observed_sky: np.array
            Numpy array representing the healpix image, centered on zenith,
            with below the horizon masked.
        """

        # Only regenerate if freq has changed
        if freq is not None:
            if np.isclose(freq, self._freq):
                pass
            else:
                self.gsm.generate(freq)
                self._freq = freq

        sky = self.gsm.generated_map_data

        # Get RA and DEC of zenith
        ra_rad, dec_rad = self.radec_of(0, np.pi/2)
        ra_deg  = ra_rad / np.pi * 180
        dec_deg = dec_rad / np.pi * 180

        # Apply rotation
        hrot = hp.Rotator(rot=[ra_deg, dec_deg], coord=['G', 'C'], inv=True)
        g0, g1 = hrot(self._theta, self._phi)
        pix0 = hp.ang2pix(self._n_side, g0, g1)
        sky_rotated = sky[pix0]

        # Generate a mask for below horizon
        mask1 = g1 + np.pi / 2 > 2 * np.pi
        mask2 = g1 < np.pi / 2
        mask = np.invert(np.logical_or(mask1, mask2))

        self.observed_sky = hp.ma(sky_rotated)
        self.observed_sky.mask = mask

        return self.observed_sky

    def view(self, logged=False, show=False, **kwargs):
        """ View the local sky, in orthographic projection.

        Parameters
        ----------
        logged: bool
            Default False, return the log2 image
        """
        sky = self.observed_sky
        if logged:
            sky = np.log2(sky)

        hp.orthview(sky, half_sky=True, **kwargs)

        if show:
            show_plt()
        return sky

    def view_observed_gsm(self, logged=False, show=False, **kwargs):
        """ View the GSM (Mollweide), with below-horizon area masked. """
        sky = self.observed_sky
        if logged:
            sky = np.log2(sky)

        # Get RA and DEC of zenith
        ra_rad, dec_rad = self.radec_of(0, np.pi / 2)
        ra_deg  = ra_rad / np.pi * 180
        dec_deg = dec_rad / np.pi * 180

        # Apply rotation
        derotate = hp.Rotator(rot=[ra_deg, dec_deg])
        g0, g1 = derotate(self._theta, self._phi)
        pix0 = hp.ang2pix(self._n_side, g0, g1)
        sky = sky[pix0]

        coordrotate = hp.Rotator(coord=['C', 'G'], inv=True)
        g0, g1 = coordrotate(self._theta, self._phi)
        pix0 = hp.ang2pix(self._n_side, g0, g1)
        sky = sky[pix0]

        hp.mollview(sky, coord='G', **kwargs)

        if show:
            show_plt()

        return sky