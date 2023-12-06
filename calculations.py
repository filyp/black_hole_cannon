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

from utils import (
    M_sun,
    R_sun,
    gravitational_acceleration,
    r_frac_steps,
    r_frac_to_density,
    r_frac_to_m_frac,
)

seconds_in_year = 365.25 * 24 * 60 * 60

# %%

# for each radius, calculate the mass eaten per s
orbital_consumptions = []
for r_frac in r_frac_steps:
    density = r_frac_to_density(r_frac)
    mass_core = r_frac_to_m_frac(r_frac) * M_sun
    speed = np.sqrt(scipy.constants.G * mass_core / (r_frac * R_sun))
    consumption_per_s_per_m2 = density * speed
    orbital_consumptions.append(consumption_per_s_per_m2)

# plot the consumption per s per m^2
import matplotlib.pyplot as plt

plt.plot(r_frac_steps, orbital_consumptions)
plt.xlabel("r/R_sun")
plt.ylabel("mass consumption in kg per s per m^2")

# max(orbital_consumptions)
max_ind = np.argmax(orbital_consumptions)
max_r_frac = r_frac_steps[max_ind]
max_r_frac
print(f"max r_frac: {max_r_frac:.8f}")

# plot vertical line at max consumption
plt.axvline(max_r_frac, color="red")
# %%
print(f"max consumption: {orbital_consumptions[max_ind]:.2e} kg/s/m^2")
# %%
# mass_bh = 1.9e27  # jupiter mass
mass_bh = 1.9e27 * 13  # heaviest known planet

schwarzschild_radius = 2 * scipy.constants.G * mass_bh / scipy.constants.c**2
print(f"Schwarzschild radius: {schwarzschild_radius:.2e} m")

# %%
k = 2
eating_disc_surface = np.pi * (k * schwarzschild_radius) ** 2
print(f"eating disc surface: {eating_disc_surface:.2e} m^2")
# %%
# kg eaten per sec
eaten_mass_per_second = orbital_consumptions[max_ind] * eating_disc_surface
print(f"eaten mass per second: {eaten_mass_per_second:.2e} kg/s")

# %%
