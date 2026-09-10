#!/usr/bin/env python3
"""Deterministic finite-model and elementary-integral checks.

Run with Python 3.12; NumPy is the only dependency.
This is not a discretization or a certified computation of continuum eigenvalues.
"""

from __future__ import annotations

import json
import math
import platform
from pathlib import Path

import numpy as np


def symmetric_norm(matrix):
    return float(np.max(np.abs(np.linalg.eigvalsh(matrix))))


def exponential_of_negative(matrix, time):
    values, vectors = np.linalg.eigh(matrix)
    return (vectors * np.exp(-time * values)) @ vectors.T


def laplacian(weights):
    return np.diag(weights.sum(axis=1)) - weights


def graph_checks():
    n = 8
    points = np.linspace(0.08, 0.92, n)
    exterior = np.array([-1.4, -0.6, 1.5, 2.2])
    exterior_weights = np.array([0.02, 0.03, 0.025, 0.02])
    s = 0.65
    c = 0.7
    kernel = np.abs(points[:, None] - exterior[None, :]) ** (-1 - 2 * s)
    m = kernel.sum(axis=0)
    edge_weights = c * kernel * exterior_weights
    local_edges = np.zeros((n, n))
    for i in range(n - 1):
        local_edges[i, i + 1] = local_edges[i + 1, i] = 1.3 + 0.1 * i
    local = laplacian(local_edges)
    interior_edges = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            interior_edges[i, j] = interior_edges[j, i] = (
                0.015 / abs(points[i] - points[j]) ** (1 + 2 * s)
            )
    interior = local + laplacian(interior_edges)
    full_ii = interior + np.diag(edge_weights.sum(axis=1))
    full_ie = -edge_weights
    full_ee = np.diag(edge_weights.sum(axis=0))
    schur = full_ii - full_ie @ np.linalg.solve(full_ee, full_ie.T)
    mean_reduction = np.zeros((n, n))
    for j in range(len(exterior)):
        mean_reduction += c * exterior_weights[j] * np.outer(kernel[:, j], kernel[:, j]) / m[j]
    a_n = full_ii - mean_reduction
    selected = np.array([1, 2])
    b = np.zeros((n, n))
    for j in selected:
        b += c * exterior_weights[j] * np.outer(kernel[:, j], kernel[:, j]) / m[j]
    eta0 = float(edge_weights[:, selected].sum(axis=1).max())
    values0, vectors0 = np.linalg.eigh(a_n)
    gap = float(values0[1])
    e = np.ones(n) / math.sqrt(n)
    projector = np.eye(n) - np.outer(e, e)
    inv_on_q = (vectors0[:, 1:] / values0[1:]) @ vectors0[:, 1:].T
    inverse_form_root = (vectors0 / np.sqrt(values0 + 1)) @ vectors0.T
    rng = np.random.default_rng(1729)
    u = rng.normal(size=n)
    free_values = kernel.T @ u / m
    full_energy = float(u @ interior @ u + np.sum(edge_weights * (u[:, None] - free_values) ** 2))
    schur_error = symmetric_norm(schur - a_n)
    energy_error = abs(full_energy - float(u @ a_n @ u))
    constant_error = float(np.linalg.norm(a_n @ e))
    alpha, time = 1.0, 0.8
    resolvent0 = np.linalg.inv(a_n + alpha * np.eye(n))
    semigroup0 = exponential_of_negative(a_n, time)
    rows = []
    assert eta0 < gap / 2
    for tau in [1.0, 0.5, 0.25, 0.125, 0.0625, 0.03125]:
        perturbation = tau * b
        eta = tau * eta0
        matrix = a_n + perturbation
        spectrum = np.linalg.eigvalsh(matrix)
        first = float(spectrum[0])
        a = float(e @ perturbation @ e)
        r = projector @ perturbation @ e
        r_squared = float(r @ r)
        quadratic = float(r @ inv_on_q @ r)
        remainder = first - a + quadratic
        remainder_bound = eta * r_squared / (gap * (gap - eta))
        first_bound = r_squared / (gap - eta)
        form_defect = symmetric_norm(inverse_form_root @ perturbation @ inverse_form_root)
        form_ground_bound = (gap + 1) * a * form_defect / (gap - a)
        resolution = symmetric_norm(np.linalg.inv(matrix + alpha * np.eye(n)) - resolvent0)
        resolution_bound = eta / (alpha * (alpha + eta))
        semigroup = symmetric_norm(exponential_of_negative(matrix, time) - semigroup0)
        shifts = spectrum - values0
        selected_free = edge_weights.copy()
        selected_free[:, selected] *= 1 - tau
        # The interior diagonal retains all original exterior edge mass.
        # Columns with zero free mass are omitted before the Schur complement.
        keep = selected_free.sum(axis=0) > 0
        free = selected_free[:, keep]
        split_schur = full_ii - (free / free.sum(axis=0)) @ free.T
        split_error = symmetric_norm(split_schur - matrix)
        tolerance = 3e-11
        assert symmetric_norm(perturbation) <= eta + tolerance
        assert np.min(shifts) >= -tolerance
        assert np.max(shifts) <= eta + tolerance
        assert 0 <= a - first + tolerance <= first_bound + 2 * tolerance
        assert a - first <= form_ground_bound + tolerance
        assert abs(remainder) <= remainder_bound + tolerance
        assert resolution <= resolution_bound + tolerance
        assert semigroup <= time * eta + tolerance
        assert split_error < tolerance
        rows.append({
            "tau": tau,
            "eta": eta,
            "lambda1": first,
            "first_order_a": a,
            "quadratic_correction_magnitude": quadratic,
            "first_order_error": a - first,
            "first_order_error_bound": first_bound,
            "relative_form_defect": form_defect,
            "form_ground_error_bound": form_ground_bound,
            "second_order_remainder": remainder,
            "second_order_remainder_bound": remainder_bound,
            "remainder_divided_by_tau_cubed": remainder / tau ** 3,
            "resolvent_error": resolution,
            "resolvent_error_bound": resolution_bound,
            "semigroup_error": semigroup,
            "semigroup_error_bound": time * eta,
            "max_eigenvalue_shift": float(np.max(shifts)),
            "split_schur_error": split_error,
        })
    assert schur_error < 3e-11 and energy_error < 3e-11 and constant_error < 3e-11
    return {
        "description": "Eight-interior-vertex weighted graph; no continuum approximation claimed.",
        "parameters": {"s": s, "c": c, "interior_points": points.tolist(), "exterior_points": exterior.tolist(), "exterior_weights": exterior_weights.tolist(), "selected_exterior_columns": selected.tolist(), "resolvent_alpha": alpha, "semigroup_time": time},
        "neumann_gap": gap,
        "local_neumann_gap": float(np.linalg.eigvalsh(local)[1]),
        "schur_identity_error": schur_error,
        "minimized_energy_identity_error": energy_error,
        "constant_nullvector_error": constant_error,
        "rows": rows,
    }


def collar_mass(s, epsilon, cutoff=0.0):
    if cutoff == 0 and s >= 0.5:
        return math.inf
    if s == 0.5:
        return 2 * math.log(epsilon * (1 + cutoff) / (cutoff * (1 + epsilon)))
    p = 1 - 2 * s
    return (epsilon ** p - cutoff ** p - (1 + epsilon) ** p + (1 + cutoff) ** p) / (s * p)


def power_integral(low, high, exponent):
    if abs(exponent + 1) < 1e-14:
        return math.log(high / low)
    q = exponent + 1
    return (high ** q - low ** q) / q


def boundary_zero_energy_exact(s, epsilon):
    first = 1 / (3 - 2 * s) - 2 / (4 - 2 * s) + 1 / (5 - 2 * s)
    second = 0.0
    for p, coefficient in [(2, 1.0), (3, -2.0), (4, 1.0)]:
        for k in range(p + 1):
            second += coefficient * math.comb(p, k) * (-epsilon) ** (p - k) * power_integral(epsilon, 1 + epsilon, k - 2 * s)
    return (first - second) / s


def boundary_zero_energy_quadrature(s, epsilon, order):
    nodes, weights = np.polynomial.legendre.leggauss(order)
    t = (nodes + 1) / 2
    x = t ** 4
    integrand = x ** 2 * (1 - x) ** 2 * (x ** (-2 * s) - (x + epsilon) ** (-2 * s)) * 4 * t ** 3 / s
    return float(np.dot(weights / 2, integrand))


def collar_checks():
    finite = []
    nodes, weights = np.polynomial.legendre.leggauss(160)
    t = (nodes + 1) / 2
    for s in [0.2, 0.35, 0.45]:
        for epsilon in [0.1, 0.01, 0.001]:
            exact = collar_mass(s, epsilon)
            leading = epsilon ** (1 - 2 * s) / (s * (1 - 2 * s))
            # Independently integrate 2*m(-r) after r=epsilon*t**q.
            # The singular term becomes t**3, avoiding endpoint quadrature error.
            alpha = 2 * s
            q = 4 / (1 - alpha)
            transformed = (epsilon ** (1 - alpha) * q * t ** 3
                           - epsilon * q * t ** (q - 1)
                           * (1 + epsilon * t ** q) ** (-alpha)) / s
            quadrature = float(np.dot(weights / 2, transformed))
            assert abs(exact - quadrature) < 2e-11
            assert 0 < exact < leading
            finite.append({"s": s, "epsilon": epsilon, "mass": exact, "independent_quadrature": quadrature, "quadrature_error": abs(exact-quadrature), "leading_asymptotic": leading, "ratio": exact / leading})
    divergent = []
    for s in [0.5, 0.65, 0.85]:
        previous = 0.0
        for cutoff in [1e-2, 1e-4, 1e-6, 1e-8]:
            mass = collar_mass(s, 0.1, cutoff)
            assert mass > previous
            previous = mass
            divergent.append({"s": s, "epsilon": 0.1, "cutoff": cutoff, "truncated_mass": mass})
    zero_profile = []
    for s in [0.25, 0.5, 0.75, 0.9]:
        for epsilon in [0.1, 0.01]:
            exact = boundary_zero_energy_exact(s, epsilon)
            quadrature = boundary_zero_energy_quadrature(s, epsilon, 160)
            refined = boundary_zero_energy_quadrature(s, epsilon, 320)
            beta_bound = math.gamma(3 - 2 * s) * math.gamma(3) / math.gamma(6 - 2 * s) / s
            assert 0 < exact < beta_bound
            assert abs(exact - refined) < 2e-11
            zero_profile.append({"s": s, "epsilon": epsilon, "energy_exact_elementary": exact, "quadrature_160": quadrature, "quadrature_320": refined, "absolute_check_error": abs(exact - refined), "finite_beta_upper_bound": beta_bound})
    return {"normalization": "c=1; Ω=(0,1); both exterior collars; u=x(1-x) for weighted energy", "subcritical_mass": finite, "truncated_divergence": divergent, "boundary_zero_profile": zero_profile}


def normalization_checks():
    rows = []
    for s in [0.2, 0.5, 0.85]:
        kappa1 = math.gamma(s + 0.5) / math.gamma(s + 0.5)
        kappa3 = math.pi * math.gamma(s + 0.5) / math.gamma(s + 1.5)
        direct3 = 2 * math.pi / (1 + 2 * s)
        assert kappa1 == 1.0
        assert abs(kappa3 - direct3) < 1e-13
        rows.append({"s": s, "kappa_n1": kappa1, "kappa_n3": kappa3,
                     "direct_radial_n3": direct3, "n3_error": abs(kappa3-direct3)})
    return rows


def write_markdown(result, output):
    graph = result["graph"]
    lines = [
        "# Reproducible verification results",
        "",
        "Finite weighted-graph sanity checks and elementary one-dimensional integrals. No continuum eigenvalue computation or numerical certification is claimed.",
        "",
        f"Python {result['runtime']['python']}; NumPy {result['runtime']['numpy']}. All assertions passed.",
        "",
        f"Graph Neumann gap: {graph['neumann_gap']:.10g}; local gap: {graph['local_neumann_gap']:.10g}.",
        f"Schur identity error: {graph['schur_identity_error']:.3e}; minimized-energy identity error: {graph['minimized_energy_identity_error']:.3e}.",
        "",
        "| τ | η | λ₁ | a−λ₁ | quadratic magnitude | second-order remainder | remainder bound |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in graph["rows"]:
        lines.append("| " + " | ".join(f"{row[key]:.7g}" for key in ["tau", "eta", "lambda1", "first_order_error", "quadratic_correction_magnitude", "second_order_remainder", "second_order_remainder_bound"]) + " |")
    lines += ["", "All eigenvalue shifts, shifted resolvent errors, semigroup errors, and exact split-vertex Schur identities satisfy the stated bounds. Full values are in `results.json`.", "", "| s | cutoff τ | truncated mass for ε=0.1 |", "|---:|---:|---:|"]
    for row in result["collar"]["truncated_divergence"]:
        lines.append(f"| {row['s']:.2f} | {row['cutoff']:.0e} | {row['truncated_mass']:.9g} |")
    lines += ["", "The exact integral proves divergence as τ↓0 for s≥1/2. The finite table illustrates that proof; it is not the justification for divergence.", "", "| s | ε | finite energy of x(1−x) | exact-versus-quadrature error |", "|---:|---:|---:|---:|"]
    for row in result["collar"]["boundary_zero_profile"]:
        lines.append(f"| {row['s']:.2f} | {row['epsilon']:.2g} | {row['energy_exact_elementary']:.10g} | {row['absolute_check_error']:.2e} |")
    output.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    output = Path(__file__).resolve().parent
    result = {
        "runtime": {"python": platform.python_version(), "numpy": np.__version__, "scipy": "not required"},
        "random_seed": 1729,
        "scope": "Finite-model sanity checks, not continuum proof or certified numerics.",
        "graph": graph_checks(),
        "tangential_normalization": normalization_checks(),
        "collar": collar_checks(),
        "status": "all assertions passed",
    }
    (output / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    write_markdown(result, output / "results.md")
    print(json.dumps({"status": result["status"], "graph_rows": len(result["graph"]["rows"]), "collar_rows": sum(len(v) for v in result["collar"].values() if isinstance(v, list)), "outputs": [str(output / "results.json"), str(output / "results.md")]}, indent=2))


if __name__ == "__main__":
    main()
