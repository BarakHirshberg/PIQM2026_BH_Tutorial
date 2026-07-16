"""Plotting helpers for the tutorial figures.

These live outside the notebook so the recipe reads as physics + analysis: each
function draws one figure from already-computed results (and the exact reference,
which it recomputes from :mod:`analysis`). They return the Matplotlib axes.
"""

import matplotlib.pyplot as plt
import numpy as np

import analysis


def plot_boson_energy_curve(bhw_points, e_points, e_err, bhw_range=(0.5, 5.5)):
    """Boson E(T): PIMD points with error bars on the exact bosonic curve.

    Energies are in units of hbar*omega0; the distinguishable curve and the
    ground-state line (4.5 hbar*omega0) are shown for reference.
    """
    omega0 = analysis.omega0()
    bhw = np.linspace(*bhw_range, 80)
    temps = [analysis.temperature_for(b) for b in bhw]

    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    ax.plot(
        bhw,
        [analysis.analytical_energy(T, "bosonic") / omega0 for T in temps],
        "-",
        color="C0",
        label="bosons (exact)",
    )
    ax.plot(
        bhw,
        [analysis.analytical_energy(T, "dist") / omega0 for T in temps],
        "--",
        color="gray",
        alpha=0.7,
        label="distinguishable (exact)",
    )
    ax.axhline(4.5, ls=":", color="k", alpha=0.5)
    ax.text(
        0.55,
        4.62,
        r"ground state $4.5\,\hbar\omega_0$",
        fontsize=8,
        color="k",
        alpha=0.7,
    )
    ax.errorbar(
        bhw_points,
        e_points,
        yerr=e_err,
        fmt="o",
        color="C0",
        ms=8,
        capsize=4,
        label="bosons (PIMD)",
    )
    ax.set_xlabel(r"$\beta\hbar\omega_0$")
    ax.set_ylabel(r"$\langle E\rangle\,/\,\hbar\omega_0$")
    ax.set_title("Three bosons in a harmonic trap")
    ax.legend()
    return ax


def plot_statistics_bars(labels, sim, err, ref, title):
    """Grouped bar chart of PIMD vs exact total energies (mHa) at one temperature."""
    x = np.arange(len(labels))
    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    ax.bar(
        x - 0.2,
        np.array(sim) * 1e3,
        0.4,
        yerr=np.array(err) * 1e3,
        capsize=4,
        label="PIMD",
    )
    ax.bar(x + 0.2, np.array(ref) * 1e3, 0.4, label="exact")
    ax.set_xticks(x, labels)
    ax.set_ylabel("total energy / mHa")
    ax.set_title(title)
    ax.legend()
    return ax


def plot_fermion_ensemble(traj_energies, mean, err, exact):
    """One-run-vs-many fermion picture: scattered per-trajectory energies (grey)
    and the sign-weighted mean with its error bar (blue), against the exact line."""
    traj = np.asarray(traj_energies)
    n = len(traj)
    jitter = 0.15 * (np.arange(n) - (n - 1) / 2) / max((n - 1) / 2, 1)

    fig, ax = plt.subplots(figsize=(6, 4), constrained_layout=True)
    ax.axhline(exact * 1e3, color="k", ls="--", label=f"exact ({exact * 1e3:.3f} mHa)")
    ax.plot(
        jitter,
        traj * 1e3,
        "o",
        color="gray",
        alpha=0.6,
        label="individual trajectories",
    )
    ax.errorbar(
        [0],
        [mean * 1e3],
        yerr=[err * 1e3],
        fmt="s",
        color="C0",
        ms=9,
        capsize=5,
        label="sign-weighted mean",
    )
    ax.set_xlim(-1, 1)
    ax.set_xticks([])
    ax.set_ylabel("fermionic total energy / mHa")
    ax.set_title("Three fermions at 30 K: one run vs. eight")
    ax.legend()
    return ax
