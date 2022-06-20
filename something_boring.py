import numpy as np

from scipy.special import gamma, gammaincinv



def find_eff_intensity(t_lum, n, r_e, e=0):
    """
    Finds the intensity (surface mass density) at the effective radius (Ie) for a Sersic
    profile with known total luminosity (stellar mass).

    Parameters
    ----------
    t_lum : float
        The total luminosity/particles/mass of the simulation.
    n : float
        The Sersic index of the profile.
    r_e : float
        The effective radius.

    Returns
    -------
    float
        The intensity at the effective radius.

    Notes
    -----
        To compute the effective intensity, equation 2 from Graham & Driver
        2005 is implemented.
        https://ui.adsabs.harvard.edu/abs/2005PASA...22..118G/abstract

    """
    bn = bn_sersic(n)
    r_ep = r_e * (1 - e)
    I_e = t_lum * (bn ** (2 * n)) / (
        gamma(2 * n) * np.exp(bn) * 2 * np.pi * n * r_e * r_ep)
    return I_e


def bn_sersic(n):
    return gammaincinv(2 * n, .5)
