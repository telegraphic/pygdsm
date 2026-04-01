from datetime import datetime

import healpy as hp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.time import Time
import pytest

from pygdsm import GSMObserver, GSMObserver16, LFSMObserver
from pygdsm import GlobalSkyModel, GlobalSkyModel16, LowFrequencySkyModel, HaslamSkyModel, HaslamObserver
from pygdsm.mckay26 import McKaySkyModel, McKayObserver

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


def _setup_observer(ov):
    ov.lon = "-118.2"
    ov.lat = "37.2"
    ov.elev = 1222
    ov.date = datetime(2000, 1, 1, 12, 0)
    return ov


def test_observer_kwargs_gsm16():
    """GSMObserver16 should forward kwargs to GlobalSkyModel16."""
    ov = GSMObserver16(include_cmb=True, resolution="low", freq_unit="MHz")
    assert ov.gsm.include_cmb is True
    assert ov.gsm.resolution == "low"


def test_observer_kwargs_gsm08():
    """GSMObserver should forward kwargs to GlobalSkyModel."""
    ov = GSMObserver(basemap="haslam")
    assert ov.gsm.basemap == "haslam"


def test_observer_kwargs_lfsm():
    """LFSMObserver should forward kwargs to LowFrequencySkyModel."""
    ov = LFSMObserver(include_cmb=True)
    assert ov.gsm.include_cmb is True


def test_observer_kwargs_haslam():
    """HaslamObserver should forward kwargs to HaslamSkyModel."""
    ov = HaslamObserver(spectral_index=-2.8, freq_unit="GHz")
    assert ov.gsm.spectral_index == -2.8
    assert ov.gsm.freq_unit == "GHz"


def test_observer_kwargs_mckay():
    """McKayObserver should forward kwargs to McKaySkyModel."""
    ov = McKayObserver(resolution="low")
    assert ov.gsm.resolution == "low"


def test_observer_prebuilt_instance_gsm16():
    """GSMObserver16 should accept a pre-built GlobalSkyModel16 instance."""
    gsm = GlobalSkyModel16(include_cmb=True)
    ov = GSMObserver16(gsm=gsm)
    assert ov.gsm is gsm
    assert ov.gsm.include_cmb is True


def test_observer_prebuilt_instance_gsm08():
    """GSMObserver should accept a pre-built GlobalSkyModel instance."""
    gsm = GlobalSkyModel(basemap="haslam")
    ov = GSMObserver(gsm=gsm)
    assert ov.gsm is gsm
    assert ov.gsm.basemap == "haslam"


def test_observer_gsm16_kwargs_generates_correctly():
    """GSMObserver16 instantiated with kwargs should generate a valid map."""
    ov = _setup_observer(GSMObserver16(include_cmb=True, resolution="low"))
    sky = ov.generate(100)
    assert sky is not None
    assert np.any(np.isfinite(sky))


def test_observer_default_still_works():
    """Observers with no arguments should still work as before (backward compatibility)."""
    for cls in (GSMObserver, GSMObserver16, LFSMObserver, HaslamObserver, McKayObserver):
        ov = _setup_observer(cls())
        sky = ov.generate(100)
        assert sky is not None



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
    for _hour in _OBS_HOURS:
        test_observed_mollview_galactic(_hour)
        test_observed_mollview_ecliptic(_hour)
        test_observed_mollview_equatorial(_hour)
        test_observed_ortho(_hour)
        test_observed_ortho_85deg_horizon(_hour)
    test_generate_with_and_without_args()
