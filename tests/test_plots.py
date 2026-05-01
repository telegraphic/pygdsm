"""test_plots.py — pytest-mpl baseline image tests for sky model views."""

from datetime import datetime

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pytest

from pygdsm import (
    GlobalSkyModel,
    GlobalSkyModel16,
    GSMObserver,
    GSMObserver16,
    HaslamSkyModel,
    LFSMObserver,
    LowFrequencySkyModel,
)
from pygdsm.mckay25 import McKaySkyModel

BASELINE_DIR = "baseline"
MPL_KWARGS = dict(baseline_dir=BASELINE_DIR, remove_text=True, tolerance=10)


@pytest.mark.mpl_image_compare(**MPL_KWARGS)
def test_plot_gsm08_408mhz():
    gsm = GlobalSkyModel(freq_unit="MHz")
    gsm.generate(408)
    gsm.view()
    return plt.gcf()


@pytest.mark.mpl_image_compare(**MPL_KWARGS)
def test_plot_gsm16_408mhz():
    gsm = GlobalSkyModel16(freq_unit="MHz", data_unit="TCMB", resolution="low")
    gsm.generate(408)
    gsm.view()
    return plt.gcf()


@pytest.mark.mpl_image_compare(**MPL_KWARGS)
def test_plot_haslam_408mhz():
    gsm = HaslamSkyModel(freq_unit="MHz")
    gsm.generate(408)
    gsm.view()
    return plt.gcf()


@pytest.mark.mpl_image_compare(**MPL_KWARGS)
def test_plot_lfsm_200mhz():
    gsm = LowFrequencySkyModel(freq_unit="MHz")
    gsm.generate(200)
    gsm.view()
    return plt.gcf()


@pytest.mark.mpl_image_compare(**MPL_KWARGS)
def test_plot_mckay25_150mhz():
    gsm = McKaySkyModel(freq_unit="MHz", resolution="low")
    gsm.generate(150)
    gsm.view()
    return plt.gcf()


@pytest.mark.mpl_image_compare(**MPL_KWARGS)
def test_plot_gsm_observer_orthview():
    ov = GSMObserver()
    ov.lon = "148.263659"
    ov.lat = "-32.998370"
    ov.elev = 100
    ov.date = datetime(2000, 1, 1, 23, 0)
    ov.generate(200)
    ov.view(logged=True)
    return plt.gcf()


@pytest.mark.mpl_image_compare(**MPL_KWARGS)
def test_plot_gsm_observer_mollview():
    ov = GSMObserver()
    ov.lon = "148.263659"
    ov.lat = "-32.998370"
    ov.elev = 100
    ov.date = datetime(2000, 1, 1, 23, 0)
    ov.generate(200)
    ov.view_observed_gsm(logged=True)
    return plt.gcf()


@pytest.mark.mpl_image_compare(**MPL_KWARGS)
def test_plot_gsm16_observer_orthview():
    ov = GSMObserver16()
    ov.lon = "148.263659"
    ov.lat = "-32.998370"
    ov.elev = 100
    ov.date = datetime(2000, 1, 1, 23, 0)
    ov.generate(200)
    ov.view(logged=True)
    return plt.gcf()
