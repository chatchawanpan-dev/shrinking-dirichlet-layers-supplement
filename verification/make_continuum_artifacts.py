#!/usr/bin/env python3
"""Generate the manuscript's continuum table and PGFPlots figure from raw JSON."""
from pathlib import Path
import json

root = Path(__file__).resolve().parents[1]
data = json.loads((root / "verification/continuum_results.json").read_text())
manuscript = root / "manuscript"
figure_dir = manuscript / "figures"
figure_dir.mkdir(exist_ok=True)
table = [r"\begin{table}[htbp]", r"\centering\small",
         r"\begin{tabular}{rrrrr}", r"\toprule",
         r"$s$ & $\widehat S_{16,128}$ & $\widehat S_{24,128}$ & $\widehat S_{24,256}$ & $(t_\varepsilon-\widehat\lambda_\varepsilon)/\rho_\varepsilon^2$ \\", r"\midrule"]
for case in data["cases"]:
    s = case["s"]
    lookup = {(row["degree"], row["quadrature_order"]): row for row in data["refinement"] if row["s"] == s}
    values = [lookup[key]["S_galerkin"] for key in [(16,128), (24,128), (24,256)]]
    values.append(case["rows"][-1]["relaxation_over_rho_squared"])
    table.append(f"{s:.2f} & " + " & ".join(f"{v:.4f}" for v in values) + r" \\")
    lines = ["epsilon geometric spectral combined limit"]
    for row in case["rows"]:
        scale = row["epsilon"] if s <= 0.25 else row["rho"] ** 2
        lines.append(" ".join(f"{v:.12g}" for v in [row["epsilon"], -row["geometric_correction"] / scale,
                       row["trial_minus_lambda"] / scale, row["minus_first_correction_over_scale"], row["predicted_scaled_limit_galerkin"]]))
    (figure_dir / f"scales_{str(s).replace('.', '_')}.dat").write_text("\n".join(lines) + "\n")
table += [r"\bottomrule\end{tabular}",
          r"\caption{Galerkin estimates on $(0,1)$ with $c=1$. Here $N$ is the largest Legendre degree and $Q$ the outer quadrature order; the inner exterior moment rule uses $2Q$ nodes. The last column uses $N=24$, $Q=256$, and $\varepsilon=10^{-12}$.}",
          r"\label{tab:continuum}", r"\end{table}"]
(manuscript / "continuum_table.tex").write_text("\n".join(table) + "\n")
plot = [r"\begin{figure}[htbp]", r"\centering", r"\begin{tikzpicture}",
        r"\begin{groupplot}[group style={group size=3 by 1,horizontal sep=1.0cm},width=0.31\textwidth,height=4.9cm,xmode=log,xmin=1e-12,xmax=1e-1,xtick={1e-12,1e-6,1e-1},xlabel={$\varepsilon$},tick label style={font=\scriptsize},label style={font=\small},title style={font=\small},grid=major,grid style={gray!18},legend style={font=\scriptsize,draw=none},legend columns=4]" ]
for case in data["cases"]:
    s = case["s"]
    filename = f"figures/scales_{str(s).replace('.', '_')}.dat"
    legend_option = ",legend to name=scaleslegend" if s == data["cases"][0]["s"] else ""
    plot.append(r"\nextgroupplot[title={$s=" + f"{s:.2f}" + r"$},ymin=0" + legend_option + "]")
    for color, style, column in [("blue!75!black","thick","combined"),("green!40!black","thick,dashed","geometric"),("orange!90!black","thick,dotted","spectral"),("gray","thin,dashdotted","limit")]:
        plot.append(r"\addplot["+color+","+style+r",no marks] table[x=epsilon,y="+column+"] {"+filename+"};")
    if s == data["cases"][0]["s"]:
        plot.append(r"\legend{Combined,Geometric,Spectral,Predicted limit}")
plot += [r"\end{groupplot}", r"\end{tikzpicture}", r"\par\smallskip\ref{scaleslegend}",
         r"\caption{The Galerkin correction $L_s\rho_\varepsilon-\widehat\lambda_\varepsilon$ and its geometric and spectral components, divided by $\varepsilon$ when $s\leq1/4$ and by $\rho_\varepsilon^2$ when $s>1/4$. The dashed--dotted limits use $G_s=1/s$ and, where it enters, the computed $\widehat S_{24,256}$. All curves come from the same finite Galerkin approximation. OpenAI Codex assisted in preparing the PGFPlots source; the curves are generated from the supplied numerical data.}",
         r"\label{fig:scales}", r"\end{figure}"]
(manuscript / "continuum_plot.tex").write_text("\n".join(plot) + "\n")
print("Generated continuum_table.tex, continuum_plot.tex, and three data files.")
