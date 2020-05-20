import sys
sys.path.append("/workdata/pygdsm")

import pylab as plt
import healpy as hp
from datetime import datetime

from pygdsm import GlobalSkyModel2016, GSMObserver
from pygdsm import GlobalSkyModel
from pygdsm import GSMObserver2016


def test_compare_gsm_to_old():
    g = GlobalSkyModel2016(freq_unit='GHz', resolution='hi', data_unit='TCMB')
    d = g.generate(0.408)
    g.view()

    g_old = GlobalSkyModel(freq_unit='GHz', basemap='haslam')
    d_old = g_old.generate(0.408)
    g_old.view()

def test_observer_test():

    # Setup observatory location - in this case, Parkes Australia
    (latitude, longitude, elevation) = ('-32.998370', '148.263659', 100)
    ov = GSMObserver2016()
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

if __name__ == "__main__":
    test_compare_gsm_to_old()
    test_observer_test()