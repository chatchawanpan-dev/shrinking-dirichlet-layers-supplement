#!/usr/bin/env python3
"""Interval Galerkin illustration of the geometric and spectral correction scales.

Python + NumPy only. Shifted orthonormal Legendre polynomials form an H1 trial
space. Local stiffness is exact; nonlocal matrices use positive quadrature.
Degree and quadrature refinement measure numerical stability, not a certified
continuum error. The first eigenvalue is evaluated through a scalar Schur
complement to avoid subtracting it from a large, ill-conditioned stiffness matrix.
"""
from pathlib import Path
import json
import math
import platform
import numpy as np
from numpy.polynomial.legendre import leggauss, legvander

ROOT = Path(__file__).resolve().parent


def gauss(order):
    x, w = leggauss(order)
    return (x + 1) / 2, w / 2


def basis(x, degree):
    return legvander(2 * np.asarray(x) - 1, degree) * np.sqrt(2 * np.arange(degree + 1) + 1)


def local_stiffness(degree):
    i, j = np.meshgrid(np.arange(degree + 1), np.arange(degree + 1), indexing="ij")
    minimum = np.minimum(i, j)
    return 2 * np.sqrt((2 * i + 1) * (2 * j + 1)) * minimum * (minimum + 1) * ((i + j) % 2 == 0)


def single_mass(s, epsilon):
    p = 1 - 2 * s
    return (epsilon ** p - math.expm1(p * math.log1p(epsilon))) / (2 * s * p)


def neumann_matrix(s, degree, order):
    t, wt = gauss(order)
    x, wx = gauss(2 * order)
    phi = basis(x, degree)
    endpoint = (-1.0) ** np.arange(degree + 1) * np.sqrt(2 * np.arange(degree + 1) + 1)
    difference = phi - endpoint
    parity = (-1.0) ** np.arange(degree + 1)

    # Interior triangle: x=y+q, y=(1-q)t. The diagonal singularity cancels
    # against the two polynomial differences. q=a^(1/(2-2s)) removes its weight.
    q = t ** (1 / (2 - 2 * s))
    y = (1 - q[:, None]) * t[None, :]
    quotient = (basis(y + q[:, None], degree) - basis(y, degree)) / q[:, None, None]
    weight = (wt[:, None] * wt[None, :] * (1 - q[:, None]) / (2 - 2 * s)).ravel()
    values = quotient.reshape(-1, degree + 1)
    interior = values.T @ (weight[:, None] * values)

    def integrated_variance(kernel, normalized_mass, outer_weight):
        weighted = kernel * wx[None, :]
        moment = weighted @ difference
        first = difference.T @ ((outer_weight @ weighted)[:, None] * difference)
        second = moment.T @ ((outer_weight / normalized_mass)[:, None] * moment)
        return first - second

    # r in (0,1). Subtracting the endpoint value keeps the constant nullvector
    # exact and removes the singular constant part before numerical integration.
    r = t
    kernel = (r[:, None] + x[None, :]) ** (-1 - 2 * s)
    mass = r ** (-2 * s) * (-np.expm1(-2 * s * np.log1p(1 / r))) / (2 * s)
    near = integrated_variance(kernel, mass, wt)

    # r>1: q=1/r, followed by q=a^(1/(2s)). Factor q^(1+2s) out
    # of the whole variance so the transformed integration has a smooth weight.
    reciprocal = t ** (1 / (2 * s))
    kernel = (1 + reciprocal[:, None] * x[None, :]) ** (-1 - 2 * s)
    mass = -np.expm1(-2 * s * np.log1p(reciprocal)) / (2 * s * reciprocal)
    far = integrated_variance(kernel, mass, wt / (2 * s))
    left = near + far
    exterior = left + parity[:, None] * left * parity[None, :]
    local = local_stiffness(degree)
    matrix = local + interior + exterior
    matrix = (matrix + matrix.T) / 2
    assert np.max(np.abs(matrix[0])) < 1e-12
    assert np.linalg.eigvalsh(matrix[1:, 1:])[0] >= 9.0
    return matrix, local, x, wx, difference, endpoint, parity


def collar_matrix(s, epsilon, ingredients, order):
    _, _, x, wx, difference, endpoint, parity = ingredients
    t, wt = gauss(order)
    r = epsilon * t
    kernel = (r[:, None] + x[None, :]) ** (-1 - 2 * s)
    mass = r ** (-2 * s) * (-np.expm1(-2 * s * np.log1p(1 / r))) / (2 * s)
    moment = (kernel * wx[None, :]) @ difference
    integrated = (epsilon * wt) @ moment
    left = single_mass(s, epsilon) * np.outer(endpoint, endpoint)
    left += np.outer(endpoint, integrated) + np.outer(integrated, endpoint)
    left += moment.T @ ((epsilon * wt / mass)[:, None] * moment)
    matrix = left + parity[:, None] * left * parity[None, :]
    return (matrix + matrix.T) / 2


def principal_schur(neumann, collar):
    trial = float(collar[0, 0])
    vector = collar[1:, 0]
    complement = neumann[1:, 1:] + collar[1:, 1:]
    assert np.linalg.eigvalsh(complement)[0] > trial
    eigenvalue = trial
    for _ in range(60):
        relaxation = float(vector @ np.linalg.solve(complement - eigenvalue * np.eye(len(vector)), vector))
        updated = trial - relaxation
        if abs(updated - eigenvalue) < 3e-15 * max(1.0, trial):
            eigenvalue = updated
            break
        eigenvalue = updated
    h = -np.linalg.solve(complement - eigenvalue * np.eye(len(vector)), vector)
    relaxation = float(-vector @ h)
    residual = abs(eigenvalue - trial + relaxation)
    assert residual < 1e-11 and relaxation >= 0
    return eigenvalue, relaxation, h, residual


def run_case(s, degree, order, epsilons):
    items = neumann_matrix(s, degree, order)
    neumann, local, _, _, _, endpoint, parity = items
    b_s = 1 / (2 * s * (1 - 2 * s))
    boundary = b_s * (np.outer(endpoint, endpoint) + np.outer(parity * endpoint, parity * endpoint))
    rhs = boundary[1:, 0]
    corrector = np.linalg.solve(neumann[1:, 1:], rhs)
    spectral = float(rhs @ corrector)
    assert spectral > 0
    rows = []
    for epsilon in epsilons:
        rho = epsilon ** (1 - 2 * s)
        collar = collar_matrix(s, epsilon, items, order)
        value, relaxation, h, residual = principal_schur(neumann, collar)
        error = h / rho + corrector
        h1_error = math.sqrt(float(error @ (local[1:, 1:] + np.eye(degree)) @ error))
        # Compute the small correction directly; do not cancel two leading values.
        geometric = -math.expm1((1 - 2 * s) * math.log1p(epsilon)) / (s * (1 - 2 * s))
        first_correction = geometric - relaxation
        scale = epsilon if s <= 0.25 else rho ** 2
        target = 1 / s if s < 0.25 else (1 / s + spectral if s == 0.25 else spectral)
        rows.append({"epsilon": epsilon, "rho": rho, "lambda1_galerkin": value,
                     "constant_trial": float(collar[0, 0]), "trial_minus_lambda": relaxation,
                     "relaxation_over_rho_squared": relaxation / rho ** 2,
                     "geometric_correction": geometric,
                     "minus_first_correction_over_scale": -first_correction / scale,
                     "predicted_scaled_limit_galerkin": target,
                     "corrector_H1_error": h1_error, "schur_residual": residual})
    return {"s": s, "degree": degree, "quadrature_order": order,
            "interior_x_quadrature_order": 2 * order, "S_galerkin": spectral,
            "G_exact": 1 / s, "rows": rows}


def main():
    epsilons = [10.0 ** (-k) for k in range(1, 13)]
    cases = []
    refinement = []
    for s in [0.15, 0.25, 0.35]:
        print(f"Computing s={s}", flush=True)
        for degree, order in [(8, 128), (16, 128), (24, 128), (24, 256)]:
            result = run_case(s, degree, order, epsilons if (degree, order) == (24, 256) else [1e-12])
            refinement.append(result)
            if (degree, order) == (24, 256):
                cases.append(result)
    result = {"status": "all algebraic assertions passed", "runtime": {"python": platform.python_version(), "numpy": np.__version__},
              "domain": "(0,1)", "c": 1, "basis": "orthonormal shifted Legendre, including the constant",
              "scope": "Continuum variational Galerkin experiment with numerical quadrature. No certified continuum eigenvalue or rigorous discretization error.",
              "cases": cases, "refinement": refinement}
    (ROOT / "continuum_results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    lines = ["# Continuum Galerkin illustration", "", result["scope"], "",
             "| s | degree | quadrature | S in trial space | relaxation / rho² at epsilon=10^-12 |", "|---:|---:|---:|---:|---:|"]
    for item in refinement:
        lines.append(f"| {item['s']} | {item['degree']} | {item['quadrature_order']} | {item['S_galerkin']:.9g} | {item['rows'][-1]['relaxation_over_rho_squared']:.9g} |")
    (ROOT / "continuum_results.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
