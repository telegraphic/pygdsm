"""
test_init.py
===========

Tests for GSM init commands
"""

from pygdsm import (
    GlobalSkyModel,
    GlobalSkyModel08,
    GlobalSkyModel16,
    GSMObserver,
    GSMObserver16,
    HaslamObserver,
    HaslamSkyModel,
    LFSMObserver,
    LowFrequencySkyModel,
    init_gsm,
    init_observer,
)


def test_init():
    assert isinstance(init_gsm("gsm"), GlobalSkyModel)
    assert isinstance(init_gsm("gsm08"), GlobalSkyModel)
    assert isinstance(init_gsm("gsm16"), GlobalSkyModel16)
    assert isinstance(init_gsm("lfsm"), LowFrequencySkyModel)
    assert isinstance(init_gsm("haslam"), HaslamSkyModel)

    assert isinstance(init_observer("gsm"), GSMObserver)
    assert isinstance(init_observer("gsm08"), GSMObserver)
    assert isinstance(init_observer("gsm16"), GSMObserver16)
    assert isinstance(init_observer("lfsm"), LFSMObserver)
    assert isinstance(init_observer("haslam"), HaslamObserver)

    gsm16 = init_gsm(
        "gsm16",
        freq_unit="MHz",
        include_cmb=True,
        interpolation="PCHIP",
    )
    assert isinstance(gsm16, GlobalSkyModel16)

if __name__ == "__main__":
    test_init()
