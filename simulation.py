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
    r_frac_to_density,
    r_frac_to_m_frac,
)

seconds_in_year = 365.25 * 24 * 60 * 60

# %%
r_orbit = 0.09031  # in units of solar radius  # best orbit for consumption
mass_core = r_frac_to_m_frac(r_orbit) * M_sun
orbital_velocity = np.sqrt(scipy.constants.G * mass_core / (r_orbit * R_sun))
print(f"orbital velocity = {orbital_velocity:.2e} m/s")

# %%
# tweak these parameters:

# mass_bh = 6e24  # kg  # earth mass
# mass_bh = 1.9e27  # jupiter mass
mass_bh = 1.9e27 * 13  # heaviest known planet
# mass_bh = 1.9e27 * 90  # heaviest known brown dwarf
# source: https://www.livescience.com/space/astronomy/whats-the-largest-planet-in-the-universe

# velocity = np.array([orbital_velocity, 0])  # m/s, don't exceed 3e8 m/s
# position = np.array([0.0, r_orbit])  # in unints of solar radius

# velocity = np.array([0.0, 0.0])  # m/s, don't exceed 3e8 m/s
# position = np.array([0.0, 0.20])  # in unints of solar radius

velocity = np.array([orbital_velocity / 1, 0.0])  # m/s, don't exceed 3e8 m/s
position = np.array([-0.5, 0.5])  # in unints of solar radius

k = 2  # how much to multiply the Schwarzschild radius by, for calculating eating disc
time_step = 5  # s


# --------------------------------------------------------------------------------------
r_sun_px = 200
width = 4 * r_sun_px
height = 4 * r_sun_px
screen = pygame.display.set_mode((width, height))

pygame.font.init()
font = pygame.font.SysFont("monospace", 20)

t = 0  # s
_eaten_mass_total = 0  # needed to avoid numerical errors
_initial_mass_bh = mass_bh

run_sim = True
iteration = 0

# time.sleep(5)
start_time = time.time()
while run_sim:
    iteration += 1
    # check for quit
    for event in pygame.event.get():
        if event.type == pygame.QUIT:
            run_sim = False

    # ----------------------------------------------------------------------------------
    # the interesting stuff happens here

    # update the time
    t += time_step
    # update the velocity
    acceleration = gravitational_acceleration(position)
    velocity += acceleration * time_step
    # update the position
    position += velocity * time_step / R_sun  # position in units of solar radius

    # eat matter from the sun
    schwarzschild_radius = 2 * scipy.constants.G * mass_bh / scipy.constants.c**2
    eating_disc_surface = np.pi * (k * schwarzschild_radius) ** 2
    eaten_volume = eating_disc_surface * time_step * np.linalg.norm(velocity)
    eaten_mass = r_frac_to_density(np.linalg.norm(position)) * eaten_volume
    old_mass_bh = mass_bh
    _eaten_mass_total += eaten_mass
    mass_bh = _initial_mass_bh + _eaten_mass_total
    # momentum must be conserved
    # TODO this over long sim will introduce numerical errors, think of a better way
    velocity = velocity * old_mass_bh / mass_bh

    # ----------------------------------------------------------------------------------

    # draw the bodies
    if iteration % 10 == 0:
        # pace the animation
        _target_draw_time = start_time + (iteration / 10) * 0.02
        _to_sleep = _target_draw_time - time.time()
        # print(_to_sleep, iter)
        if _to_sleep > 0:
            time.sleep(_to_sleep)
        else:
            print("WARNING: drawing too slow", _to_sleep, iteration)

        screen.fill((0, 0, 0))
        pygame.draw.circle(screen, "yellow", (width // 2, height // 2), r_sun_px)
        x, y = position
        bh_px = r_sun_px / 20
        pygame.draw.circle(
            screen, "red", ((2 + x) * r_sun_px, (2 - y) * r_sun_px), bh_px
        )
        # print stats on screen
        texts = [
            f"mass = {mass_bh:.20e} kg",
            f"eaten mass total = {_eaten_mass_total:.2e} kg",
            f"Schw. radius = {schwarzschild_radius:.2e} m",
            f"avg eaten mass = {_eaten_mass_total/t:.2e} kg/s",
            f"    eaten mass = {eaten_mass/time_step:.2e} kg/s",
            f"t = {t / seconds_in_year:6.4f} years",
            f"r = {np.linalg.norm(position):4.2f} R_sun",
            f"v = {np.linalg.norm(velocity):5.2e} m/s",
            f"a = {np.linalg.norm(acceleration):5.2e} m/s^2",
        ]
        for i, text in enumerate(texts):
            text_surface = font.render(text, True, "white")
            text_rect = text_surface.get_rect(topleft=(25, 25 + 25 * i))
            screen.blit(text_surface, text_rect)
        pygame.display.update()


pygame.quit()
