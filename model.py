"""
Model setup utilities for Earth-like annulus thermal convection in Underworld3.
"""

import numpy as np
import sympy
import underworld3 as uw

from config import ModelConfig


def create_mesh(config: ModelConfig):
    """Create the annulus mesh using notebook-style resolution."""
    mantle_thickness = config.outer_radius - config.inner_radius
    cell_size = mantle_thickness / config.mesh_resolution

    return uw.meshing.Annulus(
        radiusOuter=config.outer_radius,
        radiusInner=config.inner_radius,
        cellSize=cell_size,
        qdegree=config.qdegree,
    )


def create_fields(mesh):
    """Create mesh variables for the model."""
    velocity = uw.discretisation.MeshVariable("V0", mesh, 2, degree=2)
    pressure = uw.discretisation.MeshVariable("p", mesh, 1, degree=1, continuous=False)
    temperature = uw.discretisation.MeshVariable("T", mesh, 1, degree=3, continuous=True)
    temperature_initial = uw.discretisation.MeshVariable("T0", mesh, 1, degree=3, continuous=True)
    return velocity, pressure, temperature, temperature_initial


def initialise_temperature(mesh, temperature, temperature_initial, config: ModelConfig):
    """
    Initialise temperature with a conductive profile plus perturbation:
        T = T_outer + ΔT (1 - r') + A ΔT sin(n θ) cos(pi r')
    """
    r, th = mesh.CoordinateSystem.R

    mantle_thickness = config.outer_radius - config.inner_radius
    r_prime = (r - config.inner_radius) / mantle_thickness
    delta_t = config.temperature_inner - config.temperature_outer

    init_t = (
        config.temperature_outer
        + delta_t * (1.0 - r_prime)
        + config.perturbation_amplitude
        * delta_t
        * sympy.sin(config.perturbation_wavenumber * th)
        * sympy.cos(np.pi * r_prime)
    )

    with mesh.access(temperature, temperature_initial):
        values = uw.function.evaluate(init_t, temperature.coords).reshape((-1, 1))
        temperature.data[...] = values
        temperature_initial.data[...] = values


def configure_stokes_solver(mesh, velocity, pressure, temperature, config: ModelConfig):
    """
    Configure the Stokes solver:
    - constant viscosity
    - radial thermal buoyancy
    - natural penalty BC on Upper/Lower boundaries
    """
    unit_rvec = mesh.CoordinateSystem.unit_e_0
    gamma_n = unit_rvec

    stokes = uw.systems.Stokes(
        mesh,
        velocityField=velocity,
        pressureField=pressure,
    )

    stokes.bodyforce = temperature.sym[0] * unit_rvec

    stokes.constitutive_model = uw.constitutive_models.ViscousFlowModel
    stokes.constitutive_model.Parameters.viscosity = config.viscosity
    stokes.tolerance = 1.0e-3
    stokes.petsc_options["fieldsplit_velocity_mg_coarse_pc_type"] = "svd"

    penalty = 1_000_000.0 * config.viscosity
    stokes.add_natural_bc(penalty * gamma_n.dot(velocity.sym) * gamma_n, "Upper")
    stokes.add_natural_bc(penalty * gamma_n.dot(velocity.sym) * gamma_n, "Lower")

    return stokes


def configure_temperature_solver(mesh, velocity, temperature, config: ModelConfig):
    """Configure the advection-diffusion solver."""
    adv_diff = uw.systems.AdvDiffusion(
        mesh,
        u_Field=temperature,
        V_fn=velocity.sym,
        order=2,
    )

    adv_diff.constitutive_model = uw.constitutive_models.DiffusionModel
    adv_diff.constitutive_model.Parameters.diffusivity = config.diffusivity

    adv_diff.add_dirichlet_bc(config.temperature_inner, "Lower")
    adv_diff.add_dirichlet_bc(config.temperature_outer, "Upper")

    adv_diff.petsc_options.setValue("snes_rtol", 0.001)
    adv_diff.petsc_options.setValue("ksp_rtol", 0.0001)

    return adv_diff
