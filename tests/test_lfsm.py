import sys
import numpy as np

sys.path.append("/workdata/pygdsm")

import pylab as plt
import healpy as hp
from datetime import datetime

from pygdsm import GlobalSkyModel16, GlobalSkyModel, LowFrequencySkyModel
from pygdsm import GSMObserver16, GSMObserver, LFSMObserver


def test_compare_gsm_to_old():

    gl = LowFrequencySkyModel(freq_unit='MHz')
    dl = gl.generate(408)
    gl.view()

    import pylab as plt
    plt.show()


def test_observer_test():

    # Setup observatory location - in this case, Parkes Australia
    (latitude, longitude, elevation) = ('-32.998370', '148.263659', 100)
    ov = LFSMObserver()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation
    ov.date = datetime(2000, 1, 1, 23, 0)

    ov.generate(200)
    d = ov.view(logged=True)

    ov = GSMObserver()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation
    ov.date = datetime(2000, 1, 1, 23, 0)

    ov.generate(200)
    d = ov.view(logged=True)
    plt.show()

def test_cmb_removal():
    g = LowFrequencySkyModel(freq_unit='MHz', include_cmb=False)
    sky_no_cmb = g.generate(400)
    g = LowFrequencySkyModel(freq_unit='MHz',  include_cmb=True)
    sky_with_cmb = g.generate(400)    

    T_cmb = (sky_with_cmb - sky_no_cmb).mean()
    print(T_cmb)
    assert np.isclose(T_cmb, 2.725)

if __name__ == "__main__":
    test_compare_gsm_to_old()
    # test_observer_test()
