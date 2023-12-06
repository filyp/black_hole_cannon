# %%
import bisect
import csv
import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pygame
import scipy

# %%
# load sun statistics
# source: http://www.sns.ias.edu/~jnb/SNdata/Export/BP2004/bp2004stdmodel.dat

# Standard Solar Model (BP2004)
# astro-ph/0402114
# Columns in the Standard Model table (below) represent:
# 1)  Mass fraction in units of the solar mass
# 2)  Radius of the zone in units of one solar radius
# 3)  Temperature in units of deg (K)
# 4)  Density in units of g/cm^3
# 5)  Pressure in units of dyn/cm^2
# 6)  Luminosity fraction in units of the solar luminosity
# 7)  X(^1H): the hydrogen mass fraction
# 8)  X(^4He): the helium 4 mass fraction
# 9)  X(^3He): the helium 3 mass fraction
# 10) X(^12C): the carbon 12 mass fraction
# 11) X(^14N): the nitrogen 14 mass fraction
# 12) X(^16O): the oxygen 16 mass fraction

sun_stats = pd.read_csv(Path.cwd() / "data" / "sun_stats.csv", delimiter=",")
M_sun = 1.99e30  # kg
R_sun = 6.96e8  # m

r_frac_steps = sun_stats["R/Rsun"].to_numpy()


# %%
def r_frac_to_density(r_frac):
    """
    Given a radius in units of one solar radius,
    return the density in units of kg/m^3
    """
    if r_frac >= r_frac_steps[-1]:
        # we're outside the star
        return 0
    index = bisect.bisect_left(r_frac_steps, r_frac)

    density_in_g_per_cubic_cm = sun_stats.loc[index, "Rho"]
    density = density_in_g_per_cubic_cm * 1e3  # kg/m^3
    return density


def r_frac_to_m_frac(r_frac):
    """
    Given a radius in units of one solar radius,
    return the mass inside that radius in units of the solar mass
    """
    if r_frac >= r_frac_steps[-1]:
        # we're outside the star
        return 1
    index = bisect.bisect_left(r_frac_steps, r_frac)
    if index == 0:
        # if we're innermost, just to be safe return 0 to avoid numerical errors
        return 0

    m_frac = sun_stats.loc[index, "M/Msun"]
    return m_frac


def gravitational_acceleration(position):
    """
    Given a position in units of one solar radius,
    return the gravitational vector acceleration in units of m/s^2
    """
    r_frac = np.linalg.norm(position)
    m_frac = r_frac_to_m_frac(r_frac)
    mass_core = m_frac * M_sun

    _epsilon = 1  # m
    acc_strength = scipy.constants.G * mass_core / (r_frac * R_sun + _epsilon) ** 2
    acc_vector = -position / np.linalg.norm(position) * acc_strength
    return acc_vector
