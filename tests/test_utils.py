import healpy as hp
import numpy as np

from pygdsm.utils import hpix2sky, sky2hpix


def test_pix2sky():
    """Small test routine for converting healpix pixel_id to and from SkyCoords"""
    NSIDE = 32
    pix = np.arange(hp.nside2npix(NSIDE))
    sc = hpix2sky(NSIDE, pix)
    pix_roundtrip = sky2hpix(NSIDE, sc)
    assert np.allclose(pix, pix_roundtrip)


if __name__ == "__main__":
    test_pix2sky()
