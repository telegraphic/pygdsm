from __future__ import annotations
from astropy.utils.data import conf, download_file
from astropy.config import get_cache_dir

# Zenodo can be slow to respond, so set timeout to 30 seconds
conf.remote_timeout = 30

FILE_HOST = "datacentral.org.au"
GSM_DATA_URL          = [
    "https://apps.datacentral.org.au/pygdsm/data/gsm_components.h5",
    "https://zenodo.org/record/3835582/files/gsm_components.h5?download=1"
    ]

GSM2016_DATA_URL      = [
    "https://apps.datacentral.org.au/pygdsm/data/gsm2016_components.h5",
    "https://zenodo.org/record/3835582/files/gsm2016_components.h5?download=1"
]

GSM2008_TEST_DATA_URL = [
    "https://apps.datacentral.org.au/pygdsm/data/gsm_fortran_test_data.h5",
    "https://zenodo.org/record/3835582/files/gsm_fortran_test_data.h5?download=1"
]

LFSM_DATA_URL         = [
    "https://apps.datacentral.org.au/pygdsm/data/lfsm.h5",
    "https://zenodo.org/record/3835582/files/lfsm.h5?download=1"
]

HASLAM_DATA_URL       = [
    "https://lambda.gsfc.nasa.gov/data/foregrounds/haslam_2014/haslam408_dsds_Remazeilles2014.fits"
]

DATA_URLS = [
            GSM_DATA_URL,
            GSM2016_DATA_URL,
            GSM2008_TEST_DATA_URL,
            LFSM_DATA_URL,
            HASLAM_DATA_URL
            ]

def download_from_url_list(url_list: list[str] | str):
    """Download from a list of URLs, trying each in turn.

    Args:
        url_list (list[str] | str): List of URLs to try, or a single URL.

    Returns:
        str: Path to the downloaded file.
    """
    if isinstance(url_list, str):
        url_list = [url_list]

    for url in url_list:
        try:
            return download_file(url, cache=True, show_progress=True)
        except Exception as e:
            print(f"Failed to download {url}: {e}")
    raise RuntimeError("All URLs failed to download.")


def download_map_data():
    """Download component data."""
    print("Checking for component data files...")
    for URL in DATA_URLS:
        download_from_url_list(URL)

    CACHE_PATH = get_cache_dir()
    print(f"Data saved in {CACHE_PATH}")