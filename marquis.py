# Reference: https://github.com/pybamm-team/PyBaMM/blob/develop/pybamm/input/parameters/lithium_ion/Marquis2019.py

import pybamm
import numpy as np


def graphite_mcmb2528_diffusivity_Dualfoil1998(sto, T):
    """
    Graphite MCMB 2528 diffusivity as a function of stochiometry, in this case the
    diffusivity is taken to be a constant. The value is taken from Dualfoil [1].

    References
    ----------
    .. [1] http://www.cchem.berkeley.edu/jsngrp/fortran.html

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stochiometry
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """

    D_ref = 3.9 * 10 ** (-14)
    E_D_s = 42770
    arrhenius = np.exp(E_D_s / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return D_ref * arrhenius


def graphite_mcmb2528_ocp_Dualfoil1998(sto):
    """
    Graphite MCMB 2528 Open-circuit Potential (OCP) as a function of the
    stochiometry. The fit is taken from Dualfoil [1]. Dualfoil states that the data
    was measured by Chris Bogatu at Telcordia and PolyStor materials, 2000. However,
    we could not find any other records of this measurment.

    References
    ----------
    .. [1] http://www.cchem.berkeley.edu/jsngrp/fortran.html
    """

    u_eq = (
        0.194
        + 1.5 * pybamm.exp(-120.0 * sto)
        + 0.0351 * pybamm.tanh((sto - 0.286) / 0.083)
        - 0.0045 * pybamm.tanh((sto - 0.849) / 0.119)
        - 0.035 * pybamm.tanh((sto - 0.9233) / 0.05)
        - 0.0147 * pybamm.tanh((sto - 0.5) / 0.034)
        - 0.102 * pybamm.tanh((sto - 0.194) / 0.142)
        - 0.022 * pybamm.tanh((sto - 0.9) / 0.0164)
        - 0.011 * pybamm.tanh((sto - 0.124) / 0.0226)
        + 0.0155 * pybamm.tanh((sto - 0.105) / 0.029)
    )

    return u_eq


def graphite_electrolyte_exchange_current_density_Dualfoil1998(
    c_e, c_s_surf, c_s_max
):
    """
    Exchange-current density for Butler-Volmer reactions between graphite and LiPF6 in
    EC:DMC.

    References
    ----------
    .. [2] http://www.cchem.berkeley.edu/jsngrp/fortran.html

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    c_s_max : :class:`pybamm.Symbol`
        Maximum particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """
    m_ref = 2 * 10 ** (-5)  # (A/m2)(m3/mol)**1.5 - includes ref concentrations
    # E_r = 37480
    # arrhenius = np.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return (
        m_ref * c_e**0.5 * c_s_surf**0.5 * (c_s_max - c_s_surf) ** 0.5
    )


def graphite_entropic_change_Moura2016(sto, c_s_max):
    """
    Graphite entropic change in open-circuit potential (OCP) at a temperature of
    298.15K as a function of the stochiometry taken from Scott Moura's FastDFN code
    [1].

    References
    ----------
    .. [1] https://github.com/scott-moura/fastDFN

    Parameters
    ----------
    sto : :class:`pybamm.Symbol`
        Stochiometry of material (li-fraction)

    """
    du_dT = (
        -1.5 * (120.0 / c_s_max) * np.exp(-120 * sto)
        + (0.0351 / (0.083 * c_s_max)) * ((np.cosh((sto - 0.286) / 0.083)) ** (-2))
        - (0.0045 / (0.119 * c_s_max)) * ((np.cosh((sto - 0.849) / 0.119)) ** (-2))
        - (0.035 / (0.05 * c_s_max)) * ((np.cosh((sto - 0.9233) / 0.05)) ** (-2))
        - (0.0147 / (0.034 * c_s_max)) * ((np.cosh((sto - 0.5) / 0.034)) ** (-2))
        - (0.102 / (0.142 * c_s_max)) * ((np.cosh((sto - 0.194) / 0.142)) ** (-2))
        - (0.022 / (0.0164 * c_s_max)) * ((np.cosh((sto - 0.9) / 0.0164)) ** (-2))
        - (0.011 / (0.0226 * c_s_max)) * ((np.cosh((sto - 0.124) / 0.0226)) ** (-2))
        + (0.0155 / (0.029 * c_s_max)) * ((np.cosh((sto - 0.105) / 0.029)) ** (-2))
    )

    return du_dT


def lico2_diffusivity_Dualfoil1998(sto, T):
    """
    LiCo2 diffusivity as a function of stochiometry, in this case the
    diffusivity is taken to be a constant. The value is taken from Dualfoil [1].

    References
    ----------
    .. [1] http://www.cchem.berkeley.edu/jsngrp/fortran.html

    Parameters
    ----------
    sto: :class:`pybamm.Symbol`
        Electrode stochiometry
    T: :class:`pybamm.Symbol`
        Dimensional temperature

    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """
    D_ref = 1 * 10 ** (-13)
    E_D_s = 18550
    arrhenius = np.exp(E_D_s / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return D_ref * arrhenius


def lico2_ocp_Dualfoil1998(sto):
    """
    Lithium Cobalt Oxide (LiCO2) Open-circuit Potential (OCP) as a a function of the
    stochiometry. The fit is taken from Dualfoil [1]. Dualfoil states that the data
    was measured by Oscar Garcia 2001 using Quallion electrodes for 0.5 < sto < 0.99
    and by Marc Doyle for sto<0.4 (for unstated electrodes). We could not find any
    other records of the Garcia measurements. Doyles fits can be found in his
    thesis [2] but we could not find any other record of his measurments.

    References
    ----------
    .. [1] http://www.cchem.berkeley.edu/jsngrp/fortran.html
    .. [2] CM Doyle. Design and simulation of lithium rechargeable batteries,
           1995.

    Parameters
    ----------
    sto : :class:`pybamm.Symbol`
       Stochiometry of material (li-fraction)

    """

    stretch = 1.062
    sto = stretch * sto

    u_eq = (
        2.16216
        + 0.07645 * pybamm.tanh(30.834 - 54.4806 * sto)
        + 2.1581 * pybamm.tanh(52.294 - 50.294 * sto)
        - 0.14169 * pybamm.tanh(11.0923 - 19.8543 * sto)
        + 0.2051 * pybamm.tanh(1.4684 - 5.4888 * sto)
        + 0.2531 * pybamm.tanh((-sto + 0.56478) / 0.1316)
        - 0.02167 * pybamm.tanh((sto - 0.525) / 0.006)
    )

    return u_eq


def lico2_electrolyte_exchange_current_density_Dualfoil1998(c_e, c_s_surf, c_s_max):
    """
    Exchange-current density for Butler-Volmer reactions between lico2 and LiPF6 in
    EC:DMC.

    References
    ----------
    .. [2] http://www.cchem.berkeley.edu/jsngrp/fortran.html

    Parameters
    ----------
    c_e : :class:`pybamm.Symbol`
        Electrolyte concentration [mol.m-3]
    c_s_surf : :class:`pybamm.Symbol`
        Particle concentration [mol.m-3]
    c_s_max : :class:`pybamm.Symbol`
        Maximum particle concentration [mol.m-3]
    T : :class:`pybamm.Symbol`
        Temperature [K]

    Returns
    -------
    :class:`pybamm.Symbol`
        Exchange-current density [A.m-2]
    """
    m_ref = 6 * 10 ** (-7)  # (A/m2)(m3/mol)**1.5 - includes ref concentrations
    # E_r = 39570
    # arrhenius = np.exp(E_r / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return (
        m_ref * c_e**0.5 * c_s_surf**0.5 * (c_s_max - c_s_surf) ** 0.5
    )


def lico2_entropic_change_Moura2016(sto, c_s_max):
    """
    Lithium Cobalt Oxide (LiCO2) entropic change in open-circuit potential (OCP) at
    a temperature of 298.15K as a function of the stochiometry. The fit is taken
    from Scott Moura's FastDFN code [1].

    References
    ----------
    .. [1] https://github.com/scott-moura/fastDFN

    Parameters
    ----------
    sto : :class:`pybamm.Symbol`
        Stochiometry of material (li-fraction)
    """
    # Since the equation for LiCo2 from this ref. has the stretch factor,
    # should this too? If not, the "bumps" in the OCV don't line up.
    stretch = 1.062
    sto = stretch * sto

    du_dT = (
        0.07645 * (-54.4806 / c_s_max) * ((1.0 / np.cosh(30.834 - 54.4806 * sto)) ** 2)
        + 2.1581 * (-50.294 / c_s_max) * ((np.cosh(52.294 - 50.294 * sto)) ** (-2))
        + 0.14169 * (19.854 / c_s_max) * ((np.cosh(11.0923 - 19.8543 * sto)) ** (-2))
        - 0.2051 * (5.4888 / c_s_max) * ((np.cosh(1.4684 - 5.4888 * sto)) ** (-2))
        - (0.2531 / 0.1316 / c_s_max) * ((np.cosh((-sto + 0.56478) / 0.1316)) ** (-2))
        - (0.02167 / 0.006 / c_s_max) * ((np.cosh((sto - 0.525) / 0.006)) ** (-2))
    )

    return du_dT


def electrolyte_diffusivity_Capiglia1999(c_e, T):
    """
    Diffusivity of LiPF6 in EC:DMC as a function of ion concentration. The original data
    is from [1]. The fit from Dualfoil [2].

    References
    ----------
    .. [1] C Capiglia et al. 7Li and 19F diffusion coefficients and thermal
    properties of non-aqueous electrolyte solutions for rechargeable lithium batteries.
    Journal of power sources 81 (1999): 859-862.
    .. [2] http://www.cchem.berkeley.edu/jsngrp/fortran.html

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature


    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """

    D_c_e = 5.34e-10 * np.exp(-0.65 * c_e / 1000)
    E_D_e = 37040
    arrhenius = np.exp(E_D_e / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return D_c_e * arrhenius


def electrolyte_conductivity_Capiglia1999(c_e, T):
    """
    Conductivity of LiPF6 in EC:DMC as a function of ion concentration. The original
    data is from [1]. The fit is from Dualfoil [2].

    References
    ----------
    .. [1] C Capiglia et al. 7Li and 19F diffusion coefficients and thermal
    properties of non-aqueous electrolyte solutions for rechargeable lithium batteries.
    Journal of power sources 81 (1999): 859-862.
    .. [2] http://www.cchem.berkeley.edu/jsngrp/fortran.html

    Parameters
    ----------
    c_e: :class:`pybamm.Symbol`
        Dimensional electrolyte concentration
    T: :class:`pybamm.Symbol`
        Dimensional temperature


    Returns
    -------
    :class:`pybamm.Symbol`
        Solid diffusivity
    """

    sigma_e = (
        0.0911
        + 1.9101 * (c_e / 1000)
        - 1.052 * (c_e / 1000) ** 2
        + 0.1554 * (c_e / 1000) ** 3
    )

    E_k_e = 34700
    arrhenius = np.exp(E_k_e / pybamm.constants.R * (1 / 298.15 - 1 / T))

    return sigma_e * arrhenius
