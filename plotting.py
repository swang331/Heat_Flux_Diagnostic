"""
Plotting utilities for the Underworld3 convection rewrite.
"""

import os

import matplotlib.pyplot as plt
import numpy as np


def _ensure_parent(save_path):
    if save_path is not None:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)


def plot_flux_timeseries(rows, save_path=None):
    """Plot boundary heat flux through time."""
    time = np.array([row[1] for row in rows])
    inner_flux = np.array([row[3] for row in rows])
    outer_flux = np.array([row[4] for row in rows])

    plt.figure(figsize=(7, 4))
    plt.plot(time, inner_flux, label="Inner boundary")
    plt.plot(time, outer_flux, label="Outer boundary")
    plt.xlabel("Time")
    plt.ylabel("Average radial heat flux")
    plt.title("Boundary heat flux evolution")
    plt.legend()
    plt.tight_layout()

    _ensure_parent(save_path)
    if save_path is not None:
        plt.savefig(save_path, dpi=300)
    plt.close()


def plot_nusselt_timeseries(nu_rows, save_path=None):
    """Plot Nusselt-style diagnostics through time."""
    time = np.array([row[1] for row in nu_rows])
    nu_inner = np.array([row[2] for row in nu_rows])
    nu_outer = np.array([row[3] for row in nu_rows])

    plt.figure(figsize=(7, 4))
    plt.plot(time, nu_inner, label="Inner boundary Nu")
    plt.plot(time, nu_outer, label="Outer boundary Nu")
    plt.xlabel("Time")
    plt.ylabel("Nusselt number")
    plt.title("Nusselt number evolution")
    plt.legend()
    plt.tight_layout()

    _ensure_parent(save_path)
    if save_path is not None:
        plt.savefig(save_path, dpi=300)
    plt.close()


def plot_rms_velocity_timeseries(rows, save_path=None):
    """Plot RMS velocity through time."""
    time = np.array([row[1] for row in rows])
    vrms = np.array([row[5] for row in rows])

    plt.figure(figsize=(7, 4))
    plt.plot(time, vrms)
    plt.xlabel("Time")
    plt.ylabel("RMS velocity")
    plt.title("RMS velocity evolution")
    plt.tight_layout()

    _ensure_parent(save_path)
    if save_path is not None:
        plt.savefig(save_path, dpi=300)
    plt.close()


def plot_boundary_flux_vs_angle(theta, flux, title, save_path=None):
    """Plot boundary heat flux as a function of azimuth angle."""
    plt.figure(figsize=(7, 4))
    plt.plot(theta, flux)
    plt.xlabel("Azimuth angle (rad)")
    plt.ylabel("Radial heat flux")
    plt.title(title)
    plt.tight_layout()

    _ensure_parent(save_path)
    if save_path is not None:
        plt.savefig(save_path, dpi=300)
    plt.close()


def plot_temperature_field(temperature, save_path=None):
    """Plot temperature field on the annulus."""
    coords = temperature.coords
    values = temperature.data[:, 0]

    plt.figure(figsize=(6, 6))
    sc = plt.scatter(coords[:, 0], coords[:, 1], c=values, s=6)
    plt.gca().set_aspect("equal")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Annular temperature field")
    plt.colorbar(sc, label="Temperature")
    plt.tight_layout()

    _ensure_parent(save_path)
    if save_path is not None:
        plt.savefig(save_path, dpi=300)
    plt.close()


def plot_temperature_streamlines(temperature, velocity, save_path=None):
    """
    Plot temperature field with velocity vectors.

    This avoids matplotlib triangulation/interpolation, which can fail on
    annulus point clouds.
    """
    tcoords = temperature.coords
    tvals = temperature.data[:, 0]

    vcoords = velocity.coords
    u = velocity.data[:, 0]
    v = velocity.data[:, 1]

    plt.figure(figsize=(7, 7))
    sc = plt.scatter(tcoords[:, 0], tcoords[:, 1], c=tvals, s=6)

    stride = max(1, len(vcoords) // 300)
    plt.quiver(
        vcoords[::stride, 0],
        vcoords[::stride, 1],
        u[::stride],
        v[::stride],
    )

    plt.gca().set_aspect("equal")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Annular convection: temperature and velocity")
    plt.colorbar(sc, label="Temperature")
    plt.tight_layout()

    _ensure_parent(save_path)
    if save_path is not None:
        plt.savefig(save_path, dpi=300)

    plt.close()
