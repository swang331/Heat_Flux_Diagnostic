"""
Run the refactored Underworld3 annulus convection workflow.
"""

import os

from config import ModelConfig
from diagnostics import (
    compute_nusselt_numbers,
    compute_radial_heat_flux,
    compute_rayleigh_number,
    compute_rms_velocity,
    evaluate_boundary_average,
    extract_boundary_flux_vs_angle,
    save_boundary_flux_csv,
    save_diagnostics_csv,
    save_nusselt_csv,
    save_summary_csv,
)
from model import (
    configure_stokes_solver,
    configure_temperature_solver,
    create_fields,
    create_mesh,
    initialise_temperature,
)
from plotting import (
    plot_boundary_flux_vs_angle,
    plot_flux_timeseries,
    plot_nusselt_timeseries,
    plot_rms_velocity_timeseries,
    plot_temperature_field,
    plot_temperature_streamlines,
)
from tests import run_basic_checks


def main() -> None:
    """Run the annulus convection workflow."""
    config = ModelConfig()
    run_basic_checks(config)

    os.makedirs(config.results_data_dir, exist_ok=True)
    os.makedirs(config.results_figure_dir, exist_ok=True)

    rayleigh_number = compute_rayleigh_number(config)
    print(f"Planet analogue: {config.earth.planet_name}")
    print(f"Estimated Rayleigh number: {rayleigh_number:.3e}")

    mesh = create_mesh(config)
    velocity, pressure, temperature, temperature_initial = create_fields(mesh)

    initialise_temperature(mesh, temperature, temperature_initial, config)

    stokes = configure_stokes_solver(
        mesh,
        velocity,
        pressure,
        temperature,
        config,
    )
    adv_diff = configure_temperature_solver(
        mesh,
        velocity,
        temperature,
        config,
    )

    conductive_flux, advective_flux, total_flux = compute_radial_heat_flux(
        mesh,
        temperature,
        velocity,
    )

    elapsed_time = 0.0
    rows = []

    for step in range(config.max_steps):
        stokes.solve(zero_init_guess=True)

        dt_est = adv_diff.estimate_dt()
        dt = min(config.max_dt, max(config.min_dt, config.timestep_safety * dt_est))

        adv_diff.solve(timestep=dt, zero_init_guess=False)
        elapsed_time += dt

        inner_flux = evaluate_boundary_average(mesh, total_flux, config.inner_radius)
        outer_flux = evaluate_boundary_average(mesh, total_flux, config.outer_radius)
        vrms = compute_rms_velocity(velocity)

        rows.append([step, elapsed_time, dt, inner_flux, outer_flux, vrms])

        if step % config.output_interval == 0:
            print(
                f"Step {step:04d} | time = {elapsed_time:.4e} | "
                f"dt = {dt:.4e} | vrms = {vrms:.4e}"
            )

    save_diagnostics_csv(
        os.path.join(config.results_data_dir, "heat_flux_history.csv"),
        rows,
    )

    nusselt_rows = compute_nusselt_numbers(rows, config)
    save_nusselt_csv(
        os.path.join(config.results_data_dir, "nusselt_history.csv"),
        nusselt_rows,
    )

    save_summary_csv(
        os.path.join(config.results_data_dir, "summary_metrics.csv"),
        config,
        rayleigh_number,
    )

    theta_inner, flux_inner = extract_boundary_flux_vs_angle(
        mesh, total_flux, config.inner_radius
    )
    theta_outer, flux_outer = extract_boundary_flux_vs_angle(
        mesh, total_flux, config.outer_radius
    )

    save_boundary_flux_csv(
        os.path.join(config.results_data_dir, "inner_flux_vs_angle.csv"),
        theta_inner,
        flux_inner,
    )
    save_boundary_flux_csv(
        os.path.join(config.results_data_dir, "outer_flux_vs_angle.csv"),
        theta_outer,
        flux_outer,
    )

    plot_flux_timeseries(
        rows,
        save_path=os.path.join(config.results_figure_dir, "heat_flux_timeseries.png"),
    )
    plot_nusselt_timeseries(
        nusselt_rows,
        save_path=os.path.join(config.results_figure_dir, "nusselt_timeseries.png"),
    )
    plot_rms_velocity_timeseries(
        rows,
        save_path=os.path.join(config.results_figure_dir, "rms_velocity_timeseries.png"),
    )
    plot_boundary_flux_vs_angle(
        theta_inner,
        flux_inner,
        "Inner boundary heat flux vs angle",
        save_path=os.path.join(config.results_figure_dir, "inner_flux_vs_angle.png"),
    )
    plot_boundary_flux_vs_angle(
        theta_outer,
        flux_outer,
        "Outer boundary heat flux vs angle",
        save_path=os.path.join(config.results_figure_dir, "outer_flux_vs_angle.png"),
    )
    plot_temperature_field(
        temperature,
        save_path=os.path.join(config.results_figure_dir, "temperature_field.png"),
    )
    plot_temperature_streamlines(
        temperature,
        velocity,
        save_path=os.path.join(config.results_figure_dir, "temperature_streamlines.png"),
    )


if __name__ == "__main__":
    main()
