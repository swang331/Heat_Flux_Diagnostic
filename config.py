"""
Configuration settings for the Underworld3 Earth-like annulus convection model.
"""

from dataclasses import dataclass, field
from pathlib import Path


PROJECT_ROOT = Path(__file__).resolve().parent


@dataclass
class EarthConfig:
    """Physical reference values for an Earth-like mantle analogue."""

    planet_name: str = "Earth"
    planet_radius_m: float = 6.3710084e6
    core_radius_m: float = 3.48e6
    gravity_m_s2: float = 9.8
    surface_temperature_K: float = 273.0
    cmb_temperature_K: float = 3600.0
    density_kg_m3: float = 3300.0
    thermal_expansivity_1_K: float = 2.0e-5
    thermal_diffusivity_m2_s: float = 1.1394e-6
    thermal_conductivity_W_mK: float = 4.7
    heat_capacity_J_kgK: float = 1250.0
    reference_viscosity_Pa_s: float = 1.02e21

    @property
    def mantle_thickness_m(self) -> float:
        """Return mantle thickness in meters."""
        return self.planet_radius_m - self.core_radius_m

    @property
    def delta_temperature_K(self) -> float:
        """Return the thermal contrast across the shell."""
        return self.cmb_temperature_K - self.surface_temperature_K


@dataclass
class ModelConfig:
    """Container for simulation and output parameters."""

    earth: EarthConfig = field(default_factory=EarthConfig)

    # Nondimensional annulus benchmark geometry
    inner_radius: float = 0.55
    outer_radius: float = 1.0
    mesh_resolution: int = 8
    qdegree: int = 3

    # Nondimensional thermal BCs
    temperature_inner: float = 1.0
    temperature_outer: float = 0.0

    # Constant-viscosity baseline, matching the notebook spirit
    viscosity: float = 1.0
    diffusivity: float = 1.0

    # Small thermal perturbation to seed convection
    perturbation_amplitude: float = 0.05
    perturbation_wavenumber: int = 3

    # Runtime
    max_steps: int = 200
    output_interval: int = 10
    timestep_safety: float = 0.25
    max_dt: float = 1.0e-3
    min_dt: float = 1.0e-6

    # Outputs in the project folder
    results_data_dir: str = str(PROJECT_ROOT / "results" / "data")
    results_figure_dir: str = str(PROJECT_ROOT / "results" / "figures")
