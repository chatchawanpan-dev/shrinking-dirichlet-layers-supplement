# File descriptions

| File | Contents and role |
|---|---|
| `README.md` | Methods, software, commands, limitations, author, support and provenance. |
| `DATA_DICTIONARY.md` | JSON fields and plot-column/scaling definitions. |
| `build.py` | Executes both numerical programs and both artifact generators. |
| `verification/requirements.txt` | Pinned numerical dependency, NumPy 2.3.5. |
| `verification/verify.py` | Finite weighted-graph, kernel-normalization and collar-integral checks. |
| `verification/results.json` | Recorded graph/kernel output, parameters, bounds and software versions. |
| `verification/continuum_galerkin.py` | Interval Galerkin calculations with degree/quadrature refinement. |
| `verification/continuum_results.json` | Recorded eigenvalues, correction scales, residuals and refinement data. |
| `verification/make_table.py` | Generates the graph table from graph JSON output. |
| `verification/make_continuum_artifacts.py` | Generates the interval table, PGFPlots source and three plot-data files. |
| `manuscript/continuum_table.tex` | Editable table of interval refinement results. |
| `manuscript/continuum_plot.tex` | Data-driven PGFPlots fragment, caption and AI-assistance disclosure. |
| `manuscript/figures/scales_0_15.dat` | Twelve scaled correction rows for s = 0.15. |
| `manuscript/figures/scales_0_25.dat` | Twelve scaled correction rows for s = 0.25. |
| `manuscript/figures/scales_0_35.dat` | Twelve scaled correction rows for s = 0.35. |
| `supplementary/graph_table.tex` | Editable graph perturbation-error and bound table. |
| `SHA256SUMS.txt` | File hashes, excluding this checksum file itself. |

The `manuscript` and `supplementary` directories contain only computational fragments and data. Running `build.py` additionally produces `verification/results.md` and `verification/continuum_results.md` as human-readable summaries. The full article and submission correspondence are not included.
