"""
Diagnostics for radial heat flux in annulus convection.
"""

import csv
import os

import numpy as np
import sympy
import underworld3 as uw

from config import ModelConfig


def compute_rayleigh_number(config: ModelConfig) -> float:
    """
    Compute the Earth-interpretation Rayleigh number:
        Ra = g * rho * alpha * DeltaT * H^3 / (eta * kappa)
    """
    earth = config.earth
    h = earth.mantle_thickness_m
    delta_t = earth.delta_temperature_K

    return (
        earth.gravity_m_s2
        * earth.density_kg_m3
        * earth.thermal_expansivity_1_K
        * delta_t
        * h**3
        / (earth.reference_viscosity_Pa_s * earth.thermal_diffusivity_m2_s)
    )


def compute_radial_heat_flux(mesh, temperature, velocity):
    """
    Compute conductive, advective, and total radial heat flux:
        q_cond = -∇T · e_r
        q_adv  = T (u · e_r)
        q_tot  = q_cond + q_adv
    """
    x_coord, y_coord = mesh.X
    radius = sympy.sqrt(x_coord**2 + y_coord**2)
    radial_unit = mesh.X / radius

    temp_scalar = temperature.sym[0]
    grad_t = sympy.Matrix([
        sympy.diff(temp_scalar, x_coord),
        sympy.diff(temp_scalar, y_coord),
    ])

    conductive_flux = -(grad_t.dot(radial_unit))
    advective_flux = temp_scalar * (velocity.sym.dot(radial_unit))
    total_flux = conductive_flux + advective_flux

    return conductive_flux, advective_flux, total_flux


def evaluate_boundary_average(mesh, flux_expression, boundary_radius, tolerance=0.03):
    """Average a flux expression near a specified boundary radius."""
    coords = mesh.data
    radius = np.sqrt(coords[:, 0] ** 2 + coords[:, 1] ** 2)
    mask = np.abs(radius - boundary_radius) < tolerance

    if not np.any(mask):
        return np.nan

    values = uw.function.evaluate(flux_expression, coords[mask]).reshape(-1)
    return float(np.mean(values))


def extract_boundary_flux_vs_angle(mesh, flux_expression, boundary_radius, tolerance=0.03):
    """Extract radial heat flux around a boundary as a function of azimuth."""
    coords = mesh.data
    radius = np.sqrt(coords[:, 0] ** 2 + coords[:, 1] ** 2)
    mask = np.abs(radius - boundary_radius) < tolerance

    if not np.any(mask):
        return np.array([]), np.array([])

    boundary_coords = coords[mask]
    theta = np.mod(np.arctan2(boundary_coords[:, 1], boundary_coords[:, 0]), 2.0 * np.pi)
    values = uw.function.evaluate(flux_expression, boundary_coords).reshape(-1)

    order = np.argsort(theta)
    return theta[order], values[order]


def compute_rms_velocity(velocity):
    """Compute RMS velocity from the current velocity field."""
    speed_sq = velocity.data[:, 0] ** 2 + velocity.data[:, 1] ** 2
    return float(np.sqrt(np.mean(speed_sq)))


def compute_nusselt_numbers(rows, config: ModelConfig):
    """Convert boundary fluxes into simple Nusselt-style diagnostics."""
    delta_t = config.temperature_inner - config.temperature_outer
    shell_thickness = config.outer_radius - config.inner_radius

    nu_rows = []
    for step, time, dt, inner_flux, outer_flux, vrms in rows:
        nu_inner = inner_flux * shell_thickness / delta_t
        nu_outer = outer_flux * shell_thickness / delta_t
        nu_rows.append([step, time, nu_inner, nu_outer])
    return nu_rows


def save_diagnostics_csv(output_path, rows):
    """Save time-dependent diagnostics to CSV."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["step", "time", "dt", "inner_flux", "outer_flux", "vrms"])
        writer.writerows(rows)


def save_nusselt_csv(output_path, nu_rows):
    """Save Nusselt-style diagnostics to CSV."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["step", "time", "nu_inner", "nu_outer"])
        writer.writerows(nu_rows)


def save_summary_csv(output_path, config: ModelConfig, rayleigh_number: float):
    """Save Earth reference values and estimated Rayleigh number."""
    earth = config.earth
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow([
            "planet",
            "planet_radius_m",
            "core_radius_m",
            "mantle_thickness_m",
            "gravity_m_s2",
            "surface_temperature_K",
            "cmb_temperature_K",
            "density_kg_m3",
            "thermal_expansivity_1_K",
            "thermal_diffusivity_m2_s",
            "thermal_conductivity_W_mK",
            "heat_capacity_J_kgK",
            "reference_viscosity_Pa_s",
            "estimated_rayleigh_number",
        ])
        writer.writerow([
            earth.planet_name,
            earth.planet_radius_m,
            earth.core_radius_m,
            earth.mantle_thickness_m,
            earth.gravity_m_s2,
            earth.surface_temperature_K,
            earth.cmb_temperature_K,
            earth.density_kg_m3,
            earth.thermal_expansivity_1_K,
            earth.thermal_diffusivity_m2_s,
            earth.thermal_conductivity_W_mK,
            earth.heat_capacity_J_kgK,
            earth.reference_viscosity_Pa_s,
            rayleigh_number,
        ])


def save_boundary_flux_csv(output_path, theta, flux):
    """Save boundary heat flux as a function of azimuth angle."""
    os.makedirs(os.path.dirname(output_path), exist_ok=True)

    with open(output_path, "w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle)
        writer.writerow(["theta_rad", "flux"])
        for angle, value in zip(theta, flux):
            writer.writerow([angle, value])
