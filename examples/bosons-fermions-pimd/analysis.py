"""Analysis helpers for the bosonic / fermionic PIMD cookbook recipe.

This module is a trimmed and modernised version of the ``analysis.py`` /
``fermion.py`` pair shipped with the original PIQM 2023 tutorial
(https://github.com/i-pi/piqm2023-tutorial).  Two things changed for i-PI 3.x:

* the analytical reference energies (distinguishable / bosonic / fermionic) are
  kept essentially verbatim -- they are exact results for non-interacting
  particles in a harmonic trap and do not depend on the i-PI version;
* the fermionic reweighting no longer needs the custom ``ExchangePotential``
  class or ``MDAnalysis``.  i-PI 3.x writes the ``fermionic_sign`` directly to
  the output file, so recovering fermionic averages is a two-line numpy
  operation (see :func:`reweighted_fermionic_energy`).
"""

import re

import numpy as np


# --------------------------------------------------------------------------- #
#  Reading i-PI output                                                         #
# --------------------------------------------------------------------------- #
def read_ipi_output(filename):
    """Read an i-PI ``*.out`` file into a dict keyed by property name.

    The header lines (``# column N --> name``) are parsed to map each named
    property onto its column, so callers can index by name rather than by a
    hard-coded column number.
    """
    regex = re.compile(r".*column *([0-9]*) *--> ([^ {]*)")

    fields = []
    cols = []
    with open(filename, "r") as f:
        for line in f:
            if line.startswith("#"):
                match = regex.match(line)
                if match is None:
                    continue
                fields.append(match.group(2))
                cols.append(int(match.group(1)) - 1)
            else:
                break  # done with the header

    raw = np.loadtxt(filename)
    if raw.ndim == 1:  # a single output row
        raw = raw[np.newaxis, :]

    return {name: raw[:, col] for name, col in zip(fields, cols)}


def mean_energies(filename, skip_steps=0):
    """Return time-averaged kinetic and potential energies from an i-PI run.

    ``virial_fq`` is minus the quantum kinetic energy in i-PI's convention, so
    we flip its sign here.  The centroid-virial estimator (``kinetic_cv``) is
    only returned when it is present in the file (it is *not* valid under
    bosonic exchange, so we do not request it for exchange runs).
    """
    o = read_ipi_output(filename)
    out = {
        "virial": np.mean(-o["virial_fq"][skip_steps:]),
        "potential": np.mean(o["potential"][skip_steps:]),
    }
    if "kinetic_td" in o:
        out["kinetic_td"] = np.mean(o["kinetic_td"][skip_steps:])
    if "kinetic_cv" in o:
        out["kinetic_cv"] = np.mean(o["kinetic_cv"][skip_steps:])
    # Total energy from the (always-valid for bound systems) virial estimator
    out["total"] = out["virial"] + out["potential"]
    return out


def reweighted_fermionic_energy(filename, skip_steps=0):
    """Recover the fermionic total energy from a bosonic run by reweighting.

    Fermionic expectation values are obtained from a bosonic simulation via

        <A>_F = <A * s> / <s>,

    where ``s`` is the fermionic sign recorded by i-PI 3.x in the
    ``fermionic_sign`` column.  This replaces the entire hand-written
    reweighting machinery (``fermion.py`` + ``MDAnalysis``) of the 2023
    tutorial.

    Returns ``(mean_sign, fermionic_total_energy)``.
    """
    o = read_ipi_output(filename)
    sign = o["fermionic_sign"][skip_steps:]
    kinetic = -o["virial_fq"][skip_steps:]
    potential = o["potential"][skip_steps:]
    total = kinetic + potential

    mean_sign = np.mean(sign)
    fermionic_total = np.mean(total * sign) / mean_sign
    return mean_sign, fermionic_total


# --------------------------------------------------------------------------- #
#  Analytical reference energies (exact, non-interacting harmonic trap)        #
# --------------------------------------------------------------------------- #
def getZk(k, bhw, dim):
    return np.power((np.exp(0.5 * k * bhw) / (np.exp(k * bhw) - 1)), dim)


def getdZk(k, bhw, dim):
    return (
        -0.5 * k * dim * getZk(k, bhw, dim)
        * (1 + np.exp(-k * bhw)) / (1 - np.exp(-k * bhw))
    )


def get_harmonic_energy(n, bhw, dim, ptcl_type="dist"):
    """Mean energy of ``n`` non-interacting particles in a harmonic trap.

    ``ptcl_type`` is one of ``dist``, ``bosonic`` or ``fermionic`` (the last
    only implemented for ``n == 3``, which is all this tutorial needs).  The
    bosonic branch uses the recursion of Alg. 4.7 in Krauth's
    *Statistical Mechanics: Algorithms and Computations*.
    """
    if ptcl_type == "dist":
        return -n * getdZk(1, bhw, dim) / getZk(1, bhw, dim)

    if ptcl_type == "bosonic":
        z_arr = np.zeros(n + 1)
        dz_arr = np.zeros(n + 1)
        z_arr[0] = 1.0
        for m in range(1, n + 1):
            sig_z = 0.0
            sig_dz = 0.0
            for j in range(m, 0, -1):
                sig_z += getZk(j, bhw, dim) * z_arr[m - j]
                sig_dz += (
                    getdZk(j, bhw, dim) * z_arr[m - j]
                    + getZk(j, bhw, dim) * dz_arr[m - j]
                )
            z_arr[m] = sig_z / m
            dz_arr[m] = sig_dz / m
        return -dz_arr[n] / z_arr[n]

    if ptcl_type == "fermionic" and n == 3:
        num = (
            -3 - np.exp(bhw) + 8 * np.exp(2 * bhw)
            + 17 * np.exp(3 * bhw) + 15 * np.exp(4 * bhw)
        )
        denom = 2 * (-1 - np.exp(bhw) + np.exp(3 * bhw) + np.exp(4 * bhw))
        return num / denom

    return 0.0


# Physical constants / model parameters (atomic units unless noted)
SPRING_CONSTANT = 1.21647924e-8  # k for hbar*omega0 = 3 meV, in Ha/Bohr^2
MASS = 1.0
DIM = 3
KELVIN_TO_HARTREE = 3.1668152e-06  # kB in Ha/K


def omega0():
    """Trap frequency omega0 = sqrt(k / m) in atomic units."""
    return np.sqrt(SPRING_CONSTANT / MASS)


def analytical_energy(temp=17.4, sys_type="dist"):
    """Analytical mean total energy (Ha) at temperature ``temp`` (K).

    ``sys_type`` is one of ``dist``, ``bosonic``, ``mixed`` (2 bosons + 1
    distinguishable) or ``fermionic``.
    """
    omega = omega0()
    beta = 1.0 / (temp * KELVIN_TO_HARTREE)
    bhw = beta * omega  # hbar = kB = 1 in these reduced units

    if sys_type == "bosonic":
        e = get_harmonic_energy(3, bhw, DIM, "bosonic")
    elif sys_type == "mixed":
        e = get_harmonic_energy(2, bhw, DIM, "bosonic") + get_harmonic_energy(
            1, bhw, DIM, "dist"
        )
    elif sys_type == "fermionic":
        e = get_harmonic_energy(3, bhw, DIM, "fermionic")
    else:
        e = get_harmonic_energy(3, bhw, DIM, "dist")

    return e * omega  # convert from units of hbar*omega back to Hartree
