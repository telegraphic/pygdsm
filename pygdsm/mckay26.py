"""mckay26.py - Scaled version of GSM 2016 across 60-350 MHz

Reference:
McKay, L. et al., Precise Measurement of the Absolute Sky Brightness at 60–350 MHz, arXiv:2509.11846v3 (2026)
"""
import numpy as np

from .base_observer import BaseObserver
from .gsm16 import GlobalSkyModel16


offset_coeffs = [
    10.0086158,
    -9.03925211,
    1523.27697,
    -7270.93946,
    -7492.34793,
    68896.0788,
    8619.02374,
    -198752.637
]

scale_coeffs = [
    1.19851,
    -0.108045,
    1.65211,
    2.86844,
    -4.04712,
    7.41508
]

def generate_mckay26_correction(freqs_ghz: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Generate the T_offset and F_scale from McKay et al (2026).

    """
    f_145 = np.log10(freqs_ghz / .145)

    T_offset = np.sum([offset_coeffs[n] * f_145**n for n in range(len(offset_coeffs))], axis=0)
    F_scale  = np.sum([scale_coeffs[n] * f_145**n for n in range(len(scale_coeffs))], axis=0)

    return T_offset, F_scale


class McKaySkyModel(GlobalSkyModel16):
    """Sky model applying the McKay et al. (2026) correction to GSM 2016.

    GSM2016 has an offset exceeding 100 K below 100 MHz and must be scaled up by a factor
    increasing from 1.2 below 200 MHz to 1.5 at 350 MHz.
    https://arxiv.org/abs/2509.11846

    This subclasses GlobalSkyModel16 and applies a per-frequency T_offset and
    F_scale correction after generation:

        T_sky = (T_sky_gsm16 - T_offset) * F_scale

    Valid over 60–350 MHz.
    """

    def __init__(self, freq_unit="MHz", resolution="hi", theta_rot=0, phi_rot=0,
                 interpolation="pchip", include_cmb=False):
        super().__init__(
            freq_unit=freq_unit,
            data_unit="TCMB",
            resolution=resolution,
            theta_rot=theta_rot,
            phi_rot=phi_rot,
            interpolation=interpolation,
            include_cmb=include_cmb,
        )
        self.name = "McKay26"

    def generate(self, freqs):
        """Generate the McKay26-corrected sky model.

        Parameters
        ----------
        freqs: float or np.array
            Frequency for which to return the sky model (in units of self.freq_unit).

        Returns
        -------
        output: np.array
            Corrected healpix sky map, or array of maps for multiple frequencies.
        """
        output = super().generate(freqs)

        freqs_arr = np.array(freqs) * units.Unit(self.freq_unit)
        freqs_ghz = freqs_arr.to("GHz").value
        if isinstance(freqs_ghz, float):
            freqs_ghz = np.array([freqs_ghz])

        if np.min(freqs_ghz) < 0.06 or np.max(freqs_ghz) > 0.35:
            raise ValueError(
                "Frequency values lie outside 60 MHz < f < 350 MHz: "
                f"{freqs_ghz}"
            )

        T_offset, F_scale = generate_mckay26_correction(freqs_ghz)

        if output.ndim == 1:
            # single frequency: output is (npix,), corrections are length-1 arrays
            output = (output - T_offset[0]) * F_scale[0]
        else:
            # multiple frequencies: output is (nfreq, npix)
            output = (output - T_offset[:, np.newaxis]) * F_scale[:, np.newaxis]

        self.generated_map_data = output
        return output


class McKayObserver(BaseObserver):
    def __init__(self):
        """Initialize the Observer object.

        Calls ephem.Observer.__init__ function and adds on gsm
        """
        super(McKayObserver, self).__init__(gsm=McKaySkyModel)