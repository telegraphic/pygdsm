import sys
sys.path.append("/workdata/pygdsm")

import pylab as plt
import healpy as hp
import numpy as np
from datetime import datetime

import pytest

from pygdsm import GlobalSkyModel16, GSMObserver16
from pygdsm import GlobalSkyModel, GSMObserver


def test_compare_gsm_to_old():
    g = GlobalSkyModel16(freq_unit='GHz', resolution='hi', data_unit='TCMB')
    d = g.generate(0.408)
    g.view()

    g_old = GlobalSkyModel(freq_unit='GHz', basemap='haslam')
    d_old = g_old.generate(0.408)
    g_old.view()


def test_observer_test():

    # Setup observatory location - in this case, Parkes Australia
    (latitude, longitude, elevation) = ('-32.998370', '148.263659', 100)
    ov = GSMObserver16()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation
    ov.date = datetime(2000, 1, 1, 23, 0)

    ov.generate(1400)
    d = ov.view(logged=True)

    ov = GSMObserver()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation
    ov.date = datetime(2000, 1, 1, 23, 0)

    ov.generate(1400)
    d = ov.view(logged=True)
    plt.show()

def test_interp():
    f = np.arange(40, 80, 5)
    for interp in ('pchip', 'cubic'):
        for SkyModel in (GlobalSkyModel, GlobalSkyModel16):
            name = str(SkyModel).strip("<>").split('.')[-1].strip("' ")
            gsm = SkyModel(freq_unit='MHz', interpolation=interp)
            d = gsm.generate(f)

            sky_spec = d.mean(axis=1)
            fit = np.poly1d(np.polyfit(f, sky_spec, 5))(f)

            plt.plot(f, sky_spec - fit, label=f'{name}: {interp}')

    plt.xlabel("Frequency [MHz]")
    plt.ylabel("Residual [K]")
    plt.legend()
    plt.show()



def test_gsm_opts():
    g = GlobalSkyModel16(freq_unit='GHz', resolution='hi', data_unit='TCMB')
    d = g.generate(0.408)
    g = GlobalSkyModel16(freq_unit='GHz', resolution='lo', data_unit='MJysr')
    d = g.generate(0.408)
    g = GlobalSkyModel16(freq_unit='MHz', resolution='lo', data_unit='TRJ')
    d = g.generate(408)

    with pytest.raises(RuntimeError):
        g.generate(5e12)

    with pytest.raises(RuntimeError):
        g = GlobalSkyModel16(resolution='oh_hai')

    with pytest.raises(RuntimeError):
        g = GlobalSkyModel16(data_unit='furlongs/fortnight')

def test_cmb_removal():
    g = GlobalSkyModel16(freq_unit='GHz', resolution='lo', data_unit='TCMB', include_cmb=False)
    sky_no_cmb = g.generate(400)
    g = GlobalSkyModel16(freq_unit='GHz', resolution='lo', data_unit='TCMB', include_cmb=True)
    sky_with_cmb = g.generate(400)    

    T_cmb = (sky_with_cmb - sky_no_cmb).mean()
    print(T_cmb)
    assert np.isclose(T_cmb, 2.725)


if __name__ == "__main__":
    test_compare_gsm_to_old()
    test_observer_test()
    test_gsm_opts()
    test_interp()