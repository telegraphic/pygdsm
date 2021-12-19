import ephem
import healpy as hp
import numpy as np
from astropy.time import Time
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
        self._time = Time(self.date.datetime())
        # Generate mapping from pix <-> angles
        self.gsm.generate(self._freq)
        self._n_pix  = hp.get_map_size(self.gsm.generated_map_data)
        self._n_side = hp.npix2nside(self._n_pix)
        self._theta, self._phi = hp.pix2ang(self._n_side, np.arange(self._n_pix))

        self._pix0 = None
        self._mask = None

    def generate(self, freq=None, obstime=None):
        """ Generate the observed sky for the observer, based on the GSM.

        Parameters
        ----------
        freq: float
            Frequency of map to generate, in units of MHz (default).
        obstime: astropy.time.Time
            Time of observation to generate

        Returns
        -------
        observed_sky: np.array
            Numpy array representing the healpix image, centered on zenith,
            with below the horizon masked.
        """

        # Check if freq or obstime has changed
        if freq is not None:
            if not np.isclose(freq, self._freq):
                self.gsm.generate(freq)
                self._freq = freq

        # Check if time has changed -- astropy allows None == Time() comparison
        if obstime == self._time or obstime is None:
            time_has_changed = False
        else:
            time_has_changed = True
            self._time = Time(obstime)  # This will catch datetimes, but Time() object should be passed
            self.date  = obstime.to_datetime()

        # Rotation is quite slow, only recompute if time has changed, or it has never been run
        if time_has_changed or self.observed_sky is None or self._pix0 is None or self._mask is None:
            # Get RA and DEC of zenith
            ra_rad, dec_rad = self.radec_of(0, np.pi/2)
            ra_deg  = ra_rad / np.pi * 180
            dec_deg = dec_rad / np.pi * 180

            # Apply rotation
            hrot = hp.Rotator(rot=[ra_deg, dec_deg], coord=['G', 'C'], inv=True)
            g0, g1 = hrot(self._theta, self._phi)
            pix0 = hp.ang2pix(self._n_side, g0, g1)

            # Generate a mask for below horizon
            mask1 = self._phi + np.pi / 2 > 2 * np.pi
            mask2 = self._phi < np.pi / 2
            mask = np.invert(np.logical_or(mask1, mask2))
            self._pix0 = pix0
            self._mask = mask

        # Apply rotation and mask
        sky = self.gsm.generated_map_data
        sky_rotated = sky[self._pix0]
        self.observed_sky = hp.ma(sky_rotated)
        self.observed_sky.mask = self._mask

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
