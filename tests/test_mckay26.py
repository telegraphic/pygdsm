import healpy as hp
import numpy as np
import pytest

from pygdsm.mckay26 import McKaySkyModel, generate_mckay26_correction

FREQS = np.linspace(0.06, 0.35, 9)

F_SCALE_EXPECTED = np.array([1.17254264, 1.24851271, 1.20509035, 1.19931498, 1.22788176,
       1.28063423, 1.35123555, 1.43627054, 1.53396276])

T_OFFSET_EXPECTED = np.array([183.90395204,  82.42881943,  13.11002324,  13.88901167,
        19.59145165,  18.25899612,  17.17209542,  18.97684782,
        15.0879146 ])


def test_generate_mckay26_correction_t_offset():
    T_offset, _ = generate_mckay26_correction(FREQS)
    np.testing.assert_allclose(T_offset, T_OFFSET_EXPECTED, rtol=1e-6)


def test_generate_mckay26_correction_f_scale():
    _, F_scale = generate_mckay26_correction(FREQS)
    np.testing.assert_allclose(F_scale, F_SCALE_EXPECTED, rtol=1e-6)


def test_mckay26_data_unit_is_tcmb():
    m = McKaySkyModel(freq_unit="MHz")
    assert m.data_unit == "TCMB"


def test_mckay26_single_freq_shape():
    m = McKaySkyModel(freq_unit="MHz", resolution="low")
    sky = m.generate(150)
    assert sky.ndim == 1
    assert sky.shape[0] == hp.nside2npix(64)


def test_mckay26_multi_freq_shape():
    m = McKaySkyModel(freq_unit="MHz", resolution="low")
    freqs = [100, 150, 200]
    sky = m.generate(freqs)
    assert sky.shape == (len(freqs), hp.nside2npix(64))


def test_mckay26_out_of_range():
    m = McKaySkyModel(freq_unit="MHz")
    with pytest.raises((RuntimeError, ValueError)):
        m.generate(1000)