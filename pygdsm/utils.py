import healpy as hp
import numpy as np
from astropy.coordinates import SkyCoord


def hpix2sky(nside: int, pix_ids: np.ndarray) -> SkyCoord:
    """Convert a healpix pixel_id into a SkyCoord

    Args:
        nside (int): Healpix NSIDE parameter
        pix_ids (np.array): Array of pixel IDs

    Returns:
        sc (SkyCoord): Corresponding SkyCoordinates
    """
    gl, gb = hp.pix2ang(nside, pix_ids, lonlat=True)
    sc = SkyCoord(gl, gb, frame="galactic", unit=("deg", "deg"))
    return sc


def sky2hpix(nside: int, sc: SkyCoord) -> np.ndarray:
    """Convert a SkyCoord into a healpix pixel_id

    Args:
        nside (int): Healpix NSIDE parameter
        sc (SkyCoord): Astropy sky coordinates array

    Returns:
        pix (np.array): Array of healpix pixel IDs
    """
    gl, gb = sc.galactic.l.to("deg").value, sc.galactic.b.to("deg").value
    pix = hp.ang2pix(nside, gl, gb, lonlat=True)
    return pix
