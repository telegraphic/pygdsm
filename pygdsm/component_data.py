from astropy.utils.data import download_file
from astropy.config import get_cache_dir_path

FILE_HOST = "datacentral.org.au"
GSM_DATA_URL = "https://apps.datacentral.org.au/pygdsm/data/gsm_components.h5"
GSM2016_DATA_URL = "https://apps.datacentral.org.au/pygdsm/data/gsm2016_components.h5"
GSM2008_TEST_DATA_URL = "https://apps.datacentral.org.au/pygdsm/data/gsm_fortran_test_data.h5"
LFSM_DATA_URL = "https://apps.datacentral.org.au/pygdsm/data/lfsm.h5"
HASLAM_DATA_URL = "https://lambda.gsfc.nasa.gov/data/foregrounds/haslam_2014/haslam408_dsds_Remazeilles2014.fits"

DATA_URLS = [
            GSM_DATA_URL,
            GSM2016_DATA_URL,
            GSM2008_TEST_DATA_URL,
            LFSM_DATA_URL,
            HASLAM_DATA_URL
            ]

def download_map_data():
    """Download component data."""
    print("Checking for component data files...")
    for URL in DATA_URLS:
        download_file(URL, cache=True, show_progress=True)
    CACHE_PATH = get_cache_dir_path()
    print(f"Data saved in {CACHE_PATH}")