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
        self._observed_ra = None
        self._observed_dec = None

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
            self.date  = obstime.to_datetime()

        # Rotation is quite slow, only recompute if time or frequency has changed, or it has never been run
        if time_has_changed or self.observed_sky is None:
            # Get RA and DEC of zenith
            ra_zen, dec_zen = self.radec_of(0, np.pi/2)
            ra_zen  *= (180 / np.pi)
            dec_zen *= (180 / np.pi)

            # Transform from Galactic coordinates to Equatorial
            rot = hp.Rotator(coord=['G', 'C'])
            eq_theta, eq_phi = rot(self._theta, self._phi)

            # Convert from Equatorial colatitude and longitude to normal RA and DEC
            dec = 90.0 - np.abs(eq_theta*(180/np.pi))
            ra = ( (eq_phi + 2*np.pi) % (2*np.pi) )*(180/np.pi)

            # Generate a mask for below horizon (distance from zenith > 90 deg)
            dist = 2 * np.arcsin( np.sqrt( np.sin((dec_zen-dec)*(np.pi/180)/2)**2 + \
                    np.cos(dec_zen*(np.pi/180))*np.cos(dec*(np.pi/180))*np.sin((ra_zen-ra)*(np.pi/180)/2)**2 ) )
            mask = dist > (np.pi/2)
            self._mask = mask

            # Apply rotation to convert from Galactic to Equatorial and center on zenith
            hrot = hp.Rotator(rot=[ra_zen, dec_zen], coord=['G', 'C'], inv=True)
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
