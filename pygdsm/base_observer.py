import healpy as hp
import numpy as np
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, SkyCoord
from astropy.time import Time

from pygdsm.plot_utils import show_plt
from pygdsm.utils import hpix2sky, sky2hpix


def _parse_angle_radians(value):
    """Parse an angle given as a degree-string or a radian-float (pyephem convention)."""
    return np.deg2rad(float(value)) if isinstance(value, str) else float(value)


class BaseObserver:
    """Observer of the Global Sky Model.

    Generates the Observed sky, for a given point on Earth.
    Applies the necessary rotations and coordinate transformations
    so that the observed 'sky' can be returned, instead of the
    full galaxy-centered GSM.

    Exposes the same .lat / .lon / .elev / .date / radec_of() interface that
    this class previously inherited from pyephem's Observer(), now backed by
    astropy.coordinates.EarthLocation instead.
    """

    def __init__(self, gsm):
        """Initialize the Observer object.

        Parameters
        ----------
        gsm: sky model class or instance
            A pre-instantiated sky model object, or a sky model class to be
            instantiated with default arguments.
        """
        self._lat = 0.0
        self._lon = 0.0
        self._elev = 0.0
        self._time = Time.now()
        self.observed_sky = None
        self.gsm = gsm() if isinstance(gsm, type) else gsm
        self._setup()

    @property
    def lat(self):
        return self._lat

    @lat.setter
    def lat(self, value):
        self._lat = float(value)

    @property
    def lon(self):
        return self._lon

    @lon.setter
    def lon(self, value):
        self._lon = float(value)

    @property
    def elev(self):
        return self._elev

    @elev.setter
    def elev(self, value):
        self._elev = float(value)

    @property
    def date(self):
        return self._time

    @date.setter
    def date(self, value):
        self._time = value if isinstance(value, Time) else Time(value)

    def radec_of(self, az, alt):
        """Compute the ICRS RA/Dec (radians) of a given topocentric az/alt."""
        location = EarthLocation.from_geodetic(lon=self.lon * u.deg, lat=self.lat * u.deg, height=self.elev * u.m)
        aa = SkyCoord(az=az * u.rad, alt=alt * u.rad, frame=AltAz(obstime=self.date, location=location))
        sc = aa.transform_to("icrs")
        return sc.ra.rad, sc.dec.rad

    def _setup(self):
        self._freq = 100
        # Generate mapping from pix <-> angles
        self.gsm.generate(self._freq)
        self._n_pix = hp.get_map_size(self.gsm.generated_map_data)
        self._n_side = hp.npix2nside(self._n_pix)
        self._theta, self._phi = hp.pix2ang(self._n_side, np.arange(self._n_pix))

        # Galactic -> Equatorial transform is fixed for the object's lifetime
        # (depends only on nside), so it's computed once here rather than on
        # every generate() call.
        rot = hp.Rotator(coord=["G", "C"])
        eq_theta, eq_phi = rot(self._theta, self._phi)
        self._dec_allsky = 90.0 - np.abs(eq_theta * (180 / np.pi))
        self._ra_allsky = ((eq_phi + 2 * np.pi) % (2 * np.pi)) * (180 / np.pi)

        self._pix0 = None
        self._mask = None
        self._horizon_elevation = 0.0
        self._observed_ra = None
        self._observed_dec = None
        self._date_cache = self.date.jd
        self._location_cache = (self.lat, self.lon)

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
            with below the horizon masked. See `.galactic_map` for the
            un-rotated galactic-frame map.
        """
        # Check to see if frequency has changed.
        if freq is not None:
            if not np.isclose(freq, self._freq):
                self.gsm.generate(freq)
                self._freq = freq

        sky = self.gsm.generated_map_data

        # Check if time has changed, either via obstime kwarg or direct assignment to self.date
        if obstime is not None:
            obstime_astropy = obstime if isinstance(obstime, Time) else Time(obstime)
            if obstime_astropy != self._time:
                time_has_changed = True
                self._time = obstime_astropy
                self.date = obstime_astropy.to_datetime()
            else:
                time_has_changed = False
        else:
            # Detect changes to self.date set directly (e.g. ov.date = datetime(...))
            time_has_changed = (self.date.jd != self._date_cache)

        if time_has_changed:
            self._date_cache = self.date.jd

        # Detect changes to self.lat / self.lon set directly
        location_cache = (self.lat, self.lon)
        location_has_changed = (location_cache != self._location_cache)
        if location_has_changed:
            self._location_cache = location_cache

        # Match pyephem convention -- string is degrees, int/float is rad
        horizon_elevation = _parse_angle_radians(horizon_elevation or 0.0)
        if self._horizon_elevation == horizon_elevation:
            horizon_has_changed = False
        else:
            self._horizon_elevation = horizon_elevation
            horizon_has_changed = True

        if self._horizon_elevation < 0:
            raise ValueError(f"Horizon elevation must be greater or equal to 0 degrees (currently {np.rad2deg(horizon_elevation)}).")

        # Rotation is quite slow, only recompute if time, location, or horizon has changed, or it has never been run
        if time_has_changed or location_has_changed or self.observed_sky is None or horizon_has_changed:
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

            # Apply rotation to convert from Galactic to Equatorial and center on zenith
            hrot = hp.Rotator(rot=[ra_zen, dec_zen], coord=["G", "C"], inv=True)
            g0, g1 = hrot(self._theta, self._phi)
            pix0 = hp.ang2pix(self._n_side, g0, g1)
            self._pix0 = pix0

            self._observed_ra = self._ra_allsky[self._pix0]
            self._observed_dec = self._dec_allsky[self._pix0]

        sky_rotated = sky[self._pix0]
        mask_rotated = self._mask[self._pix0]

        self.observed_sky = hp.ma(sky_rotated)
        self.observed_sky.mask = mask_rotated

        return self.observed_sky

    @property
    def galactic_map(self):
        """The un-rotated galactic-frame map from the last generate() call.

        Unlike `observed_sky`, this is not rotated to the observer's zenith
        or masked below the horizon -- it's the raw GSM output, equivalent to
        calling `self.gsm.generate(freq)` directly.
        """
        return self.gsm.generated_map_data

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
