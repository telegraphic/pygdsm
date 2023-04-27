from pygdsm import GSMObserver, GSMObserver16, LFSMObserver
import pylab as plt
import healpy as hp
from datetime import datetime
import numpy as np
import os
from astropy.time import Time

def test_gsm_observer(show=False):
    """ Test GSMObserver() is working
    """
    (latitude, longitude, elevation) = ('37.2', '-118.2', 1222)
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

    (latitude, longitude, elevation) = ('37.2', '-118.2', 1222)
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

    (latitude, longitude, elevation) = ('37.2', '-118.2', 1222)
    ov = LFSMObserver()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation
    ov.date = datetime(2000, 1, 1, 23, 0)
    ov.generate(50)
    ov.view(logged=True)
    ov.view_observed_gsm(logged=True)
    if show:
        plt.show()

def test_observed_mollview():
    """ Generate animated maps of observing coverage over 24 hours """

    (latitude, longitude, elevation) = ('37.2', '-118.2', 1222)
    ov = GSMObserver()
    ov.lon = longitude
    ov.lat = latitude
    ov.elev = elevation

    obs = []
    if not os.path.exists('generated_sky'):
        os.mkdir('generated_sky')
    for ii in range(0, 24, 4):
        ov.date = datetime(2000, 1, 1, ii, 0)
        ov.generate(50)
        sky = ov.view_observed_gsm(logged=True, show=False, min=9, max=20)
        plt.savefig('generated_sky/galactic-%02d.png' % ii)
        plt.close()

        hp.mollview(sky, coord=['G', 'E'], min=9, max=20)
        plt.savefig('generated_sky/ecliptic-%02d.png' % ii)
        plt.close()

        hp.mollview(sky, coord=['G', 'C'], min=9, max=20)
        plt.savefig('generated_sky/equatorial-%02d.png' % ii)
        plt.close()

        ov.view(logged=True, show=False, min=9, max=20)
        plt.savefig('generated_sky/ortho-%02d.png' % ii)
        plt.close()

        print(ii)

    os.system('convert -delay 20 generated_sky/ortho-*.png ortho.gif')
    os.system('convert -delay 20 generated_sky/galactic-*.png galactic.gif')
    os.system('convert -delay 20 generated_sky/ecliptic-*.png ecliptic.gif')
    os.system('convert -delay 20 generated_sky/equatorial-*.png equatorial.gif')


def test_generate_with_and_without_args():
    """ Test generating without frequency argument """
    (latitude, longitude, elevation) = ('37.2', '-118.2', 1222)
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

if __name__ == "__main__":
    test_gsm_observer(show=True)
    test_observed_mollview()
    test_generate_with_and_without_args()
