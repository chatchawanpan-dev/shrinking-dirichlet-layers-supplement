#!/usr/bin/env python3
"""Render the manuscript table from recorded verification results."""
import json
from pathlib import Path

root = Path(__file__).resolve().parents[1]
data = json.loads((root / "verification/results.json").read_text())
def tex(value):
    mantissa, exponent = f"{value:.3e}".split("e")
    return rf"${mantissa}\times10^{{{int(exponent)}}}$"
lines = [r"\begin{table}[htbp]", r"\centering", r"\small",
         r"\begin{tabular}{rrrr}", r"\toprule",
         r"$\tau$ & $a_\tau-\lambda_1(\tau)$ & $|\mathcal R_\tau|$ & Bound for $|\mathcal R_\tau|$ \\",
         r"\midrule"]
for row in data["graph"]["rows"]:
    lines.append(f"{row['tau']:.5g} & {tex(row['first_order_error'])} & {tex(abs(row['second_order_remainder']))} & {tex(row['second_order_remainder_bound'])}" + r" \\")
lines += [r"\bottomrule", r"\end{tabular}",
          r"\caption{Errors for the finite weighted-graph comparison, with perturbation $\tau B$. The second-order remainder and its bound are defined in~\eqref{eq:graph-remainder}.  All bounds use the graph's own Neumann gap and defect parameter.}",
          r"\label{tab:graph-checks}", r"\end{table}",
          f"The computations used Python {data['runtime']['python']} and NumPy {data['runtime']['numpy']}. The single random vector used for the energy identity has fixed seed {data['random_seed']}; all graph weights and perturbations are deterministic."]
(root / "supplementary/graph_table.tex").write_text("\n".join(lines) + "\n")
