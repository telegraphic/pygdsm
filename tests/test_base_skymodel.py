"""
test_base_skymodel.py
=====================

Tests for GSM base skymodel
"""

import pytest
from pygdsm.base_skymodel import BaseSkyModel
from pygdsm.component_data import GSM_FILEPATH

from pathlib import Path


def test_base_skymodel_init():

    gsm = BaseSkyModel('TEST_GSM', GSM_FILEPATH,
                       freq_unit='MHz', basemap='haslam', data_unit='K')

    with pytest.raises(RuntimeError):
        gsm = BaseSkyModel('TEST_GSM', 'not_a_file.madeup',
                       freq_unit='MHz', basemap='haslam', data_unit='K')

    with pytest.raises(RuntimeError):
        gsm = BaseSkyModel('TEST_GSM', Path(__file__).absolute(),
                       freq_unit='MHz', basemap='haslam', data_unit='K')

if __name__ == "__main__":
    test_base_skymodel_init()