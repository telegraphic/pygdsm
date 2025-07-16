import ephem
import healpy as hp
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time

from pygdsm.plot_utils import show_plt
from pygdsm.utils import hpix2sky, sky2hpix


class BaseObserver(ephem.Observer):
    """Observer of the Global Sky Model.

    Generates the Observed sky, for a given point on Earth.
    Applies the necessary rotations and coordinate transformations
    so that the observed 'sky' can be returned, instead of the
    full galaxy-centered GSM.

    This class is based on pyephem's Observer(). The GSM bit can be thought of
    as an 'add on' to ephem.Observer, adding the methods generate() and  view(),
    which allows all-sky images for a given point on earth to be produced.
    """

    def __init__(self, gsm):
        """Initialize the Observer object.

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
        self._n_pix = hp.get_map_size(self.gsm.generated_map_data)
        self._n_side = hp.npix2nside(self._n_pix)
        self._theta, self._phi = hp.pix2ang(self._n_side, np.arange(self._n_pix))

        self._pix0 = None
        self._mask = None
        self._horizon_elevation = 0.0
        self._observed_ra = None
        self._observed_dec = None

    def generate(self, freq=None, obstime=None, horizon_elevation=None):
        """ Generate the observed sky for the observer, based on the GSM.

        Parameters
        ----------
        freq: float
            Frequency of map to generate, in units of MHz (default).
        obstime: astropy.time.Time
            Time of observation to generate
        horizon_elevation: float
            Elevation of the artificial horizon (default 0.0)

        Returns
        -------
        observed_sky: np.array
            Numpy array representing the healpix image, centered on zenith,
            with below the horizon masked.
        """
        # Check to see if frequency has changed.
        if freq is not None:
            if not np.isclose(freq, self._freq):
                self.gsm.generate(freq)
                self._freq = freq

        sky = self.gsm.generated_map_data

        # Check if time has changed -- astropy allows None == Time() comparison
        if obstime == self._time or obstime is None:
            time_has_changed = False
        else:
            time_has_changed = True
            self._time = Time(obstime)  # This will catch datetimes, but Time() object should be passed
            self.date = obstime.to_datetime()

        # Match pyephem convertion -- string is degrees, int/float is rad
        horizon_elevation = ephem.degrees(horizon_elevation or 0.0)
        if self._horizon_elevation == horizon_elevation:
            horizon_has_changed = False
        else:
            self._horizon_elevation = horizon_elevation
            horizon_has_changed = True

        # Checking this separately encase someone tries to be smart and sets both the attribute and kwarg
        if self._horizon_elevation < 0:
            raise ValueError(f"Horizon elevation must be greater or equal to 0 degrees (currently {np.rad2deg(horizon_elevation)}).")

        # Match pyephem convertion -- string is degrees, int/float is rad
        horizon_elevation = ephem.degrees(horizon_elevation or 0.0)
        if self._horizon_elevation == horizon_elevation:
            horizon_has_changed = False
        else:
            self._horizon_elevation = horizon_elevation
            horizon_has_changed = True

        # Checking this separately encase someone tries to be smart and sets both the attribute and kwarg
        if self._horizon_elevation < 0:
            raise ValueError(f"Horizon elevation must be greater or equal to 0 degrees (currently {np.rad2deg(horizon_elevation)}).")

        # Rotation is quite slow, only recompute if time or frequency has changed, or it has never been run
        if time_has_changed or self.observed_sky is None or horizon_has_changed:
            # Get RA and DEC of zenith
            ra_zen, dec_zen = self.radec_of(0, np.pi / 2)
            sc_zen = SkyCoord(ra_zen, dec_zen, unit=("rad", "rad"))
            pix_zen = sky2hpix(self._n_side, sc_zen)
            vec_zen = hp.pix2vec(self._n_side, pix_zen)

            # Convert to degrees
            ra_zen *= 180 / np.pi
            dec_zen *= 180 / np.pi

            # Generate below-horizon mask using query_disc
            mask = np.ones(shape=self._n_pix, dtype='bool')
            pix_visible = hp.query_disc(self._n_side, vec=vec_zen, radius=np.pi/2 - self._horizon_elevation)
            mask[pix_visible] = 0
            self._mask = mask

            # Transform from Galactic coordinates to Equatorial
            rot = hp.Rotator(coord=["G", "C"])
            eq_theta, eq_phi = rot(self._theta, self._phi)

            # Convert from Equatorial colatitude and longitude to normal RA and DEC
            dec = 90.0 - np.abs(eq_theta * (180 / np.pi))
            ra = ((eq_phi + 2 * np.pi) % (2 * np.pi)) * (180 / np.pi)

            # Apply rotation to convert from Galactic to Equatorial and center on zenith
            hrot = hp.Rotator(rot=[ra_zen, dec_zen], coord=["G", "C"], inv=True)
            g0, g1 = hrot(self._theta, self._phi)
            pix0 = hp.ang2pix(self._n_side, g0, g1)
            self._pix0 = pix0

            dec_rotated = dec[self._pix0]
            ra_rotated = ra[self._pix0]
            self._observed_ra = ra_rotated
            self._observed_dec = dec_rotated

        sky_rotated = sky[self._pix0]
        mask_rotated = self._mask[self._pix0]

        self.observed_sky = hp.ma(sky_rotated)
        self.observed_sky.mask = mask_rotated

        return self.observed_sky

    def view(self, logged=False, show=False, **kwargs):
        """View the local sky, in orthographic projection.

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

    @property
    def observed_gsm(self):
        """Return the GSM (Mollweide), with below-horizon area masked."""
        sky = self.observed_sky

        # Get RA and DEC of zenith
        ra_rad, dec_rad = self.radec_of(0, np.pi / 2)
        ra_deg = ra_rad / np.pi * 180
        dec_deg = dec_rad / np.pi * 180

        # Apply rotation
        derotate = hp.Rotator(rot=[ra_deg, dec_deg])
        g0, g1 = derotate(self._theta, self._phi)
        pix0 = hp.ang2pix(self._n_side, g0, g1)
        sky = sky[pix0]

        coordrotate = hp.Rotator(coord=["C", "G"], inv=True)
        g0, g1 = coordrotate(self._theta, self._phi)
        pix0 = hp.ang2pix(self._n_side, g0, g1)
        sky = sky[pix0]
        return sky

    def view_observed_gsm(self, logged=False, show=False, **kwargs):
        """View the GSM (Mollweide), with below-horizon area masked.

        Args:
            logged (bool): Apply log2 to data (default False)
            show (bool): Call plt.show() (default False)

        Returns:
            sky (np.array): Healpix map of observed GSM.
        """
        sky = self.observed_gsm

        if logged:
            sky = np.log2(sky)

        hp.mollview(sky, coord="G", **kwargs)
        if show:
            show_plt()
        return sky
