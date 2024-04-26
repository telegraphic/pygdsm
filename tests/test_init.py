"""
test_init.py
===========

Tests for GSM init commands
"""

from pygdsm import GlobalSkyModel, GlobalSkyModel08, GlobalSkyModel16, LowFrequencySkyModel, HaslamSkyModel
from pygdsm import GSMObserver, GSMObserver16, LFSMObserver, HaslamObserver
from pygdsm import init_gsm, init_observer


def test_init():
    assert isinstance(init_gsm('gsm'), GlobalSkyModel)
    assert isinstance(init_gsm('gsm08'), GlobalSkyModel)
    assert isinstance(init_gsm('gsm16'), GlobalSkyModel16)
    assert isinstance(init_gsm('lfsm'), LowFrequencySkyModel)
    assert isinstance(init_gsm('haslam'), HaslamSkyModel)

    assert isinstance(init_observer('gsm'), GSMObserver)
    assert isinstance(init_observer('gsm08'), GSMObserver)
    assert isinstance(init_observer('gsm16'), GSMObserver16)
    assert isinstance(init_observer('lfsm'), LFSMObserver)
    assert isinstance(init_observer('haslam'), HaslamObserver)


if __name__ == "__main__":
    test_init()