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


# --- Tests for issue #38: observer location correctness ---

def _make_observer_at(lat, lon, obs_date):
    """Helper: create a GSMObserver16 at a given lat/lon string."""
    ov = GSMObserver16()
    ov.lat = lat
    ov.lon = lon
    ov.elev = 0
    ov.date = obs_date
    return ov


def test_different_locations_produce_different_skies():
    """Observers at widely separated lat/lon should see different skies (issue #38).

    Uses the correct pyephem attributes: .lat and .lon.
    """
    obs_date = datetime(2025, 12, 22, 15, 0)
    ov_south = _make_observer_at("-30", "-30", obs_date)
    ov_north = _make_observer_at("60", "60", obs_date)

    sky_south = ov_south.generate(1000)
    sky_north = ov_north.generate(1000)

    # Compare full healpix maps (always the same size regardless of masking)
    assert not np.allclose(sky_south.data, sky_north.data), (
        "Observers at different locations produced identical skies"
    )


def test_changing_location_invalidates_cache():
    """Updating .lat/.lon on an existing observer should produce a new sky (issue #38).

    The generate() cache is keyed on time and horizon, but NOT on lat/lon.
    This test will FAIL if location changes are not detected.
    """
    obs_date = datetime(2025, 12, 22, 15, 0)
    ov = GSMObserver16()
    ov.lat = "-30"
    ov.lon = "-30"
    ov.elev = 0
    ov.date = obs_date
    sky_first = ov.generate(1000).data.copy()

    # Move the observer to a very different location
    ov.lat = "60"
    ov.lon = "60"
    sky_second = ov.generate(1000).data.copy()

    assert not np.allclose(sky_first, sky_second), (
        "Changing .lat/.lon did not invalidate the generate() cache"
    )


def test_wrong_attribute_names_are_ignored_by_ephem():
    """Setting .latitude/.longitude (not .lat/.lon) is silently ignored by pyephem (issue #38).

    pyephem's Observer uses .lat and .lon.  Setting .latitude/.longitude creates
    plain Python instance attributes that ephem never reads, so the observer
    position stays at its default (0 deg) regardless of the values assigned.
    This test documents the footgun: two observers 'at' different latitudes via
    the wrong attrs will produce the same sky as an observer at lat=0.
    """
    obs_date = datetime(2025, 12, 22, 15, 0)

    # Observer using CORRECT attrs
    ov_correct = GSMObserver16()
    ov_correct.lat = "45"
    ov_correct.lon = "0"
    ov_correct.elev = 0
    ov_correct.date = obs_date
    sky_correct = ov_correct.generate(1000).data.copy()

    # Observer using WRONG attrs (replicates the issue reporter's code)
    ov_wrong = GSMObserver16()
    ov_wrong.latitude = "45"   # ignored by ephem -- stays at default lat
    ov_wrong.longitude = "0"   # ignored by ephem -- stays at default lon
    ov_wrong.elev = 0
    ov_wrong.date = obs_date
    sky_wrong = ov_wrong.generate(1000).data.copy()

    # The wrong-attr observer will NOT be at lat=45; its effective lat is the
    # ephem default (0).  So the skies should differ.
    assert not np.allclose(sky_correct, sky_wrong), (
        ".latitude/.longitude appear to be respected by ephem (unexpected) — "
        "if this assertion fails the footgun may have been fixed upstream"
    )

if __name__ == "__main__":
    test_gsm_observer(show=True)
    for _hour in _OBS_HOURS:
        test_observed_mollview_galactic(_hour)
        test_observed_mollview_ecliptic(_hour)
        test_observed_mollview_equatorial(_hour)
        test_observed_ortho(_hour)
        test_observed_ortho_85deg_horizon(_hour)
    test_generate_with_and_without_args()
