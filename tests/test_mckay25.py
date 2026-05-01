import healpy as hp
import numpy as np
import pytest

from pygdsm.mckay25 import McKaySkyModel, generate_mckay25_correction

FREQS = np.linspace(0.06, 0.35, 9)

F_SCALE_EXPECTED = np.array([1.17254264, 1.24851271, 1.20509035, 1.19931498, 1.22788176,
       1.28063423, 1.35123555, 1.43627054, 1.53396276])

T_OFFSET_EXPECTED = np.array([183.90395204,  82.42881943,  13.11002324,  13.88901167,
        19.59145165,  18.25899612,  17.17209542,  18.97684782,
        15.0879146 ])


def test_generate_mckay25_correction_t_offset():
    T_offset, _ = generate_mckay25_correction(FREQS)
    np.testing.assert_allclose(T_offset, T_OFFSET_EXPECTED, rtol=1e-6)


def test_generate_mckay25_correction_f_scale():
    _, F_scale = generate_mckay25_correction(FREQS)
    np.testing.assert_allclose(F_scale, F_SCALE_EXPECTED, rtol=1e-6)


def test_mckay25_data_unit_is_tcmb():
    m = McKaySkyModel(freq_unit="MHz")
    assert m.data_unit == "TCMB"


def test_mckay25_single_freq_shape():
    m = McKaySkyModel(freq_unit="MHz", resolution="low")
    sky = m.generate(150)
    assert sky.ndim == 1
    assert sky.shape[0] == hp.nside2npix(64)


def test_mckay25_multi_freq_shape():
    m = McKaySkyModel(freq_unit="MHz", resolution="low")
    freqs = [100, 150, 200]
    sky = m.generate(freqs)
    assert sky.shape == (len(freqs), hp.nside2npix(64))


def test_mckay25_out_of_range():
    m = McKaySkyModel(freq_unit="MHz")
    with pytest.raises((RuntimeError, ValueError)):
        m.generate(1000)


# --- CMB handling ---

def test_include_cmb_false_is_true_minus_2725_single_freq():
    """With include_cmb=False, output must equal include_cmb=True output minus 2.725 K."""
    m_with = McKaySkyModel(freq_unit="MHz", resolution="low", include_cmb=True)
    m_without = McKaySkyModel(freq_unit="MHz", resolution="low", include_cmb=False)
    sky_with = m_with.generate(150)
    sky_without = m_without.generate(150)
    np.testing.assert_allclose(sky_without, sky_with - 2.725, rtol=1e-6)


def test_include_cmb_false_is_true_minus_2725_multi_freq():
    """Same check over multiple frequencies; subtraction is broadcast across pixels."""
    freqs = [100, 150, 200]
    m_with = McKaySkyModel(freq_unit="MHz", resolution="low", include_cmb=True)
    m_without = McKaySkyModel(freq_unit="MHz", resolution="low", include_cmb=False)
    sky_with = m_with.generate(freqs)
    sky_without = m_without.generate(freqs)
    np.testing.assert_allclose(sky_without, sky_with - 2.725, rtol=1e-6)


def test_cmb_subtracted_after_scale_not_before():
    """Verify CMB is subtracted AFTER the scale correction, not before.

    If the 2.725 K were subtracted before (i.e. from the raw GSM16 map), the
    difference between include_cmb=True and include_cmb=False outputs would be
    2.725 * F_scale, not 2.725.  We confirm it is exactly 2.725.
    """
    freq_mhz = 150.0
    freqs_ghz = np.array([freq_mhz / 1e3])
    _, F_scale = generate_mckay25_correction(freqs_ghz)
    # F_scale at 150 MHz is noticeably different from 1.0
    assert not np.isclose(F_scale[0], 1.0, atol=0.01), "F_scale unexpectedly 1 — test is trivial"

    m_with = McKaySkyModel(freq_unit="MHz", resolution="low", include_cmb=True)
    m_without = McKaySkyModel(freq_unit="MHz", resolution="low", include_cmb=False)
    diff = m_with.generate(freq_mhz) - m_without.generate(freq_mhz)
    np.testing.assert_allclose(diff, 2.725, rtol=1e-6)


def test_include_cmb_attribute_restored_after_generate():
    """generate() must restore self.include_cmb to its original value."""
    m_false = McKaySkyModel(freq_unit="MHz", resolution="low", include_cmb=False)
    m_false.generate(150)
    assert m_false.include_cmb is False

    m_true = McKaySkyModel(freq_unit="MHz", resolution="low", include_cmb=True)
    m_true.generate(150)
    assert m_true.include_cmb is True