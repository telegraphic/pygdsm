from datetime import datetime

import healpy as hp
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
import pytest

from pygdsm import GSMObserver, GSMObserver16, LFSMObserver

matplotlib.use("Agg")

BASELINE_DIR = "baseline"
MPL_KWARGS = dict(baseline_dir=BASELINE_DIR, remove_text=True, tolerance=10)

_OBS_HOURS = [0, 4, 8, 12, 16, 20]


def _make_gsm_observer(hour):
    ov = GSMObserver()
    ov.lon = "-118.2"
    ov.lat = "37.2"
    ov.elev = 1222
    ov.date = datetime(2000, 1, 1, hour, 0)
    ov.generate(50)
    return ov


def test_gsm_observer(show=False):
    """Test GSMObserver() is working"""
    (latitude, longitude, elevation) = ("37.2", "-118.2", 1222)
    ov = GSMObserver()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation
    ov.date = datetime(2000, 1, 1, 23, 0)
    ov.generate(50)
    ov.view(logged=True)
    ov.view_observed_gsm(logged=True)
    if show:
        plt.show()

    (latitude, longitude, elevation) = ("37.2", "-118.2", 1222)
    ov = GSMObserver16()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation
    ov.date = datetime(2000, 1, 1, 23, 0)
    ov.generate(50)
    ov.view(logged=True)
    ov.view_observed_gsm(logged=True)
    if show:
        plt.show()

    (latitude, longitude, elevation) = ("37.2", "-118.2", 1222)
    ov = LFSMObserver()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation
    ov.date = datetime(2000, 1, 1, 23, 0)
    ov.generate(50)
    ov.view(logged=True)
    ov.view_observed_gsm(logged=True)

    with pytest.raises(ValueError) as e:
        ov.generate(horizon_elevation=-1e-3)

    if show:
        plt.show()


@pytest.mark.parametrize("hour", _OBS_HOURS)
@pytest.mark.mpl_image_compare(**MPL_KWARGS)
def test_observed_mollview_galactic(hour):
    """Galactic-frame Mollweide view of GSM with horizon mask."""
    ov = _make_gsm_observer(hour)
    ov.view_observed_gsm(logged=True, min=9, max=20)
    return plt.gcf()


@pytest.mark.parametrize("hour", _OBS_HOURS)
@pytest.mark.mpl_image_compare(**MPL_KWARGS)
def test_observed_mollview_ecliptic(hour):
    """Ecliptic-frame Mollweide view of GSM with horizon mask."""
    ov = _make_gsm_observer(hour)
    sky = np.log2(ov.observed_gsm)
    hp.mollview(sky, coord=["G", "E"], min=9, max=20)
    return plt.gcf()


@pytest.mark.parametrize("hour", _OBS_HOURS)
@pytest.mark.mpl_image_compare(**MPL_KWARGS)
def test_observed_mollview_equatorial(hour):
    """Equatorial-frame Mollweide view of GSM with horizon mask."""
    ov = _make_gsm_observer(hour)
    sky = np.log2(ov.observed_gsm)
    hp.mollview(sky, coord=["G", "C"], min=9, max=20)
    return plt.gcf()


@pytest.mark.parametrize("hour", _OBS_HOURS)
@pytest.mark.mpl_image_compare(**MPL_KWARGS)
def test_observed_ortho(hour):
    """Orthographic view of observer sky."""
    ov = _make_gsm_observer(hour)
    ov.view(logged=True, min=9, max=20)
    return plt.gcf()


@pytest.mark.parametrize("hour", _OBS_HOURS)
@pytest.mark.mpl_image_compare(**MPL_KWARGS)
def test_observed_ortho_85deg_horizon(hour):
    """Orthographic view with 85-degree artificial horizon."""
    ov = _make_gsm_observer(hour)
    ov.generate(freq=50, horizon_elevation='85.0')
    ov.view(logged=True, min=9, max=20)
    return plt.gcf()



def test_generate_with_and_without_args():
    """Test generating without frequency argument"""
    (latitude, longitude, elevation) = ("37.2", "-118.2", 1222)
    ov = GSMObserver()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation
    ov.date = datetime(2000, 1, 1, 23, 0)

    ov.generate(50)

    ov.date = datetime(2000, 1, 1, 22, 0)
    ov.generate()

    ov.generate(51)
    ov.generate(52)
    ov.generate()
    now = Time(datetime.now())
    ov.generate(obstime=now)
    ov.generate(obstime=now)
    ov.generate(obstime=now, freq=53)
    ov.generate(obstime=now, freq=52)
    ov.generate(obstime=now, freq=53, horizon_elevation=0.0)
    ov.generate(obstime=now, freq=52, horizon_elevation='0.0')
    ov.generate(obstime=now, freq=53, horizon_elevation=np.deg2rad(85.0))
    ov.generate(obstime=now, freq=52, horizon_elevation='85.0')

if __name__ == "__main__":
    test_gsm_observer(show=True)
    test_observed_mollview()
    test_generate_with_and_without_args()
