from .gsm08 import GlobalSkyModel, GSMObserver
from .gsm16 import GlobalSkyModel16, GSMObserver16
from .haslam import HaslamObserver, HaslamSkyModel
from .lfsm import LFSMObserver, LowFrequencySkyModel
from .component_data import download_map_data as download_map_data

# Add aliases
GlobalSkyModel08 = GlobalSkyModel
GSMObserver08 = GSMObserver

def init_gsm(gsm_name: str = "gsm08", *args, **kwargs):
    """Initialize a GDSM object by ID/name

    Returns a diffuse sky model (subclass of BaseSkyModel), based on one of:
      * **GSM08:** A model of diffuse Galactic radio emission from 10 MHz to 100 GHz,
                   [Oliveira-Costa et. al., (2008)](https://ui.adsabs.harvard.edu/abs/2008MNRAS.388..247D/abstract).
      * **GSM16:** An improved model of diffuse galactic radio emission from 10 MHz to 5 THz,
                   [Zheng et. al., (2016)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.464.3486Z/abstract).
      * **LFSM:** The LWA1 Low Frequency Sky Survey (10-408 MHz)
                  [Dowell et. al. (2017)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.469.4537D/abstract).
      * **Haslam:** A frequency-scaled model using a spectral index, based on the Haslam 408 MHz map
                   [Haslam 408 MHz](https://lambda.gsfc.nasa.gov/product/foreground/fg_2014_haslam_408_info.cfm).

    Args:
        gsm_name (str): One of 'gsm08', 'gsm16', 'lfsm' or 'haslam'

    Returns:
        sky_model (various): Corresponding sky model
    """
    gsm_name = gsm_name.lower().strip()

    if gsm_name == 'gsm': # Shorthand for GSM08
        return GlobalSkyModel(*args, **kwargs)
    elif gsm_name == 'gsm08':
        return GlobalSkyModel(*args, **kwargs)
    elif gsm_name == 'gsm16':
        return GlobalSkyModel16(*args, **kwargs)
    elif gsm_name == 'lfsm':
        return LowFrequencySkyModel(*args, **kwargs)
    elif gsm_name == 'haslam':
        return HaslamSkyModel(*args, **kwargs)
    else:
        raise ValueError(f'Invalid model specification "{gsm_name}"')



def init_observer(gsm_name: str = "gsm08", *args, **kwargs):
    """Initialize a GDSM Observer object by ID/name

    Returns an observer (subclass of BaseObserver), where the diffuse sky is created from one of:
      * **GSM08:** A model of diffuse Galactic radio emission from 10 MHz to 100 GHz,
                   [Oliveira-Costa et. al., (2008)](https://ui.adsabs.harvard.edu/abs/2008MNRAS.388..247D/abstract).
      * **GSM16:** An improved model of diffuse galactic radio emission from 10 MHz to 5 THz,
                   [Zheng et. al., (2016)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.464.3486Z/abstract).
      * **LFSM:** The LWA1 Low Frequency Sky Survey (10-408 MHz)
                  [Dowell et. al. (2017)](https://ui.adsabs.harvard.edu/abs/2017MNRAS.469.4537D/abstract).
      * **Haslam:** A frequency-scaled model using a spectral index, based on the Haslam 408 MHz map
                   [Haslam 408 MHz](https://lambda.gsfc.nasa.gov/product/foreground/fg_2014_haslam_408_info.cfm).

    Args:
        gsm_name (str): One of 'gsm08', 'gsm16', 'lfsm' or 'haslam'

    Returns:
        observer (various): Corresponding sky model observer
    """
    gsm_name = gsm_name.lower().strip()

    if gsm_name == 'gsm': # Shorthand for GSM08
        return GSMObserver(*args, **kwargs)
    elif gsm_name == 'gsm08':
        return GSMObserver(*args, **kwargs)
    elif gsm_name == 'gsm16':
        return GSMObserver16(*args, **kwargs)
    elif gsm_name == 'lfsm':
        return LFSMObserver(*args, **kwargs)
    elif gsm_name == 'haslam':
        return HaslamObserver(*args, **kwargs)
    else:
        raise ValueError(f'Invalid model specification "{gsm_name}"')
