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

def download_map_data():
    """Download component data."""
    print("Checking for component data files...")
    for URL in DATA_URLS:
        try:
            download_file(URL[0], cache=True, show_progress=True)
        except Exception as e:
            print(f"Failed to download {URL[0]}: {e}. Trying alternative URL...")
            try:
                # Try the alternative URL
                download_file(URL[1], cache=True, show_progress=True)
            except Exception as e:
                print(f"Failed to download {URL[1]}: {e}")
                raise RuntimeError(f"Failed to download data from both URLs: {URL[0]} and {URL[1]}")

    CACHE_PATH = get_cache_dir()
    print(f"Data saved in {CACHE_PATH}")