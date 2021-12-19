from astropy.utils.data import conf, download_file

# Zenodo can be slow to respond, so set timeout to 30 seconds
conf.remote_timeout = 30

GSM_FILEPATH      = download_file('https://zenodo.org/record/3835582/files/gsm_components.h5?download=1',
                                  cache=True, show_progress=True)

GSM2016_FILEPATH  = download_file('https://zenodo.org/record/3835582/files/gsm2016_components.h5?download=1',
                                  cache=True, show_progress=True)

GSM2008_TEST_DATA = download_file('https://zenodo.org/record/3835582/files/gsm_fortran_test_data.h5?download=1',
                                  cache=True, show_progress=True)

LFSM_FILEPATH = download_file('https://zenodo.org/record/3835582/files/lfsm.h5?download=1',
                              cache=True, show_progress=True)

HASLAM_FILEPATH = download_file('https://lambda.gsfc.nasa.gov/data/foregrounds/haslam_2014/haslam408_dsds_Remazeilles2014.fits',
                                cache=True, show_progress=True)
