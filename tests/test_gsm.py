"""
test_gsm.py
===========

Tests for GSM module.
"""

import os
import time

import h5py
import healpy as hp
import numpy as np
import pytest
from astropy.coordinates import SkyCoord

from pygdsm import GlobalSkyModel
from pygdsm.component_data import GSM2008_TEST_DATA_URL, download_from_url_list


# Download test data as needed (stored in astropy cache)
GSM2008_TEST_DATA = download_from_url_list(GSM2008_TEST_DATA_URL)

def test_gsm_generate():
    """Test maps generate successfully, and that output shapes are correct."""

    freq = 40
    gsm = GlobalSkyModel()
    map = gsm.generate(freq)
    assert map.shape == (3145728,)

    freqs = [40, 80, 120, 160]
    gsm = GlobalSkyModel(interpolation="cubic")
    map = gsm.generate(freqs)
    assert map.shape == (4, 3145728)

    freqs = np.linspace(1, 90, 10)
    freq_unit = "GHz"
    for map_id in ["5deg", "haslam", "wmap"]:
        gsm = GlobalSkyModel(basemap=map_id, interpolation="pchip", freq_unit=freq_unit)
        map = gsm.generate(freqs)
        assert map.shape == (10, 3145728)


def test_compare_to_gsm():
    """Compare output of python version to fortran version.

    Compares against output of original GSM. Note that the interpolation
    functions used in this and in the original differ. So, the outputs
    differ too, by as much as 5%. This just goes to show that there's only
    so much accuracy that one can get via interpolation.

    If the frequency matches a map used to generate the GSM, then the
    fortran and python should match up pretty much exactly as no interpolation
    is required. In between maps is where the differences will become more apparent.

    Note: The gsm.f that is supplied only ever uses the haslam map. It also
          has no option to change the interpolation method.

    Note: Because each output
    """
    gsm = GlobalSkyModel(basemap="haslam", interpolation="pchip")
    gsm_fortran = h5py.File(GSM2008_TEST_DATA, "r")
    for freq in (10, 100, 408, 1000, 1420, 2326, 10000, 23000, 40000, 90000, 94000):
        print("\nComparing at %i MHz..." % freq)
        a = gsm_fortran["map_%imhz.out" % freq][:]
        b = gsm.generate(freq)
        frac_avg = np.average(np.abs(1 - a / b))
        if frac_avg > 0.01:
            print("FORTRAN: ", a)
            print("PYTHON:  ", b)
            print(f"FRAC AVG: {frac_avg:.6f}")
        assert frac_avg < 0.035


def test_set_methods():
    gsm = GlobalSkyModel()
    gsm.generate(40)
    # gsm.view()

    gsm.set_basemap("haslam")
    gsm.view()

    gsm.set_basemap("5deg")
    gsm.view()

    gsm.set_basemap("wmap")
    gsm.view()

    gsm.set_freq_unit("GHz")
    gsm.view()


def test_speed():
    gsm = GlobalSkyModel(basemap="haslam")

    t1 = time.time()
    for ff in np.linspace(10, 10000, 100):
        gsm.generate(ff)
    t2 = time.time()
    print("Time taken: %2.2fs" % (t2 - t1))


def test_write_fits():
    gsm = GlobalSkyModel()
    gsm.generate(1000)
    gsm.write_fits("test_write_fits.fits")

    d_fits = hp.read_map("test_write_fits.fits")
    d_gsm = gsm.generated_map_data

    assert d_fits.shape == d_gsm.shape
    assert np.allclose(d_fits, d_gsm)

    os.remove("test_write_fits.fits")


def test_cmb_removal():
    g = GlobalSkyModel(freq_unit="MHz", include_cmb=False)
    sky_no_cmb = g.generate(400)
    g = GlobalSkyModel(freq_unit="MHz", include_cmb=True)
    sky_with_cmb = g.generate(400)

    T_cmb = (sky_with_cmb - sky_no_cmb).mean()
    print(T_cmb)
    assert np.isclose(T_cmb, 2.725)


def test_get_sky_temperature():
    gc = SkyCoord(0, 0, unit="deg", frame="galactic")
    freqs = (50, 100, 150)
    g = GlobalSkyModel()
    T = g.get_sky_temperature(gc, freqs)

    T_gold = np.array([258368.7463149, 50466.45058671, 18968.12555978])
    assert np.allclose(T, T_gold)


def test_stupid_values():
    with pytest.raises(RuntimeError):
        g = GlobalSkyModel(basemap="haslamalan")
    with pytest.raises(RuntimeError):
        g = GlobalSkyModel(interpolation="linear")
    with pytest.raises(RuntimeError):
        g = GlobalSkyModel()
        g.generate(9999999999)


def test_set_interpolation_method():
    g = GlobalSkyModel(interpolation="pchip")
    assert g.interpolation_method == "pchip"
    g.set_interpolation_method("cubic")
    assert g.interpolation_method == "cubic"


if __name__ == "__main__":
    test_stupid_values()
    test_set_interpolation_method()
    test_gsm_generate()
    test_compare_to_gsm()
    test_speed()
    test_write_fits()
    test_set_methods()
    test_cmb_removal()
    test_get_sky_temperature()
