# Separation of Azeotropic Mixtures using Ionic Liquids

Bachelor's thesis project (FCHPT) by **Aliaksandr Baran**.

The work models the separation of azeotropic mixtures (e.g. MTBE–methanol) using
ionic liquids as entrainers. It covers vapour–liquid (VLE) and liquid–liquid (LLE)
equilibria described with the **NRTL** activity-coefficient model, together with the
material/energy balances and economics of the proposed separation scheme
(extractor, evaporator, valve, IL regeneration).

## Migrating the simulation Excel → MATLAB → Python

The calculations were first built in **Excel**, then ported step by step. All
three implementations are now complete:

```
Excel  ──►  MATLAB  ──►  Python
(reference)  (complete)   (complete)
```

1. **Excel** — the original spreadsheets with all VLE/LLE/balance/economics
   calculations. This is the reference implementation.
2. **MATLAB** — modular re-implementation: NRTL/LLE/VLE solvers plus the full
   process model (extractor, valve, evaporators, heat exchangers, pumps,
   economics) orchestrated by `matlab/main_simulation.m`.
3. **Python** — clean, scriptable port in the `python/separation/` package, with
   an improved variant (`python/separation_improved/`) that replaces the damped
   iteration with a simultaneous Newton solver. See [python/README.md](python/README.md).

### Solvent scenarios (clean vs regenerated IL)

The process is run for two extractor-solvent cases — identical flowsheet, only the
solvent composition `x_S = [MTBE, MeOH, IL]` differs:

| Scenario | `x_S` | Where |
|----------|-------|-------|
| **Clean (fresh) IL** | `[0, 0, 1]` | `python/main_simulation_improved.py clean`, `matlab/main_simulation_clean.m` |
| **Regenerated IL** | `[0.001559, 0.019959, 0.978483]` | `python/main_simulation.py`, `python/main_simulation_improved.py regen`, `matlab/main_simulation.m` |

## Repository structure

```
.
├── excel/                 Original Excel implementation (reference)
│   ├── hlavna-schema/      Main process scheme (columns, evaporator, valve, VLE)
│   ├── nova-schema/        Revised scheme + material balance & economics
│   ├── regenerovana-IL/    Variant with regenerated ionic liquid
│   ├── LLE/                Liquid–liquid equilibrium sheets (incl. LLE macro)
│   └── skusim/             Working drafts / NRTL ternary VLE sheet
├── matlab/                MATLAB re-implementation
│   ├── main_simulation.m        Full process simulation (regenerated IL)
│   ├── main_simulation_clean.m  Full process simulation (clean/fresh IL)
│   ├── separation_process.m     Process driver
│   ├── bacalar_project.m        LLE equilibrium script (NRTL, ternary)
│   ├── LLE_solver.m             LLE_solver() + NRTL gamma calculation
│   └── functions/               Modular unit-operation functions
├── python/                Python port (separation/ package + tests)
│   ├── separation/              Core modules (NRTL, LLE, VLE, units, economics)
│   ├── separation_improved/     Newton-solver / Brent-VLE improved variant
│   ├── main_simulation.py       Orchestrator (regenerated IL, Excel-verified)
│   └── main_simulation_improved.py  Improved runner (both solvent scenarios)
└── docs/                  Author's own documents
    ├── thesis/             Thesis text, tables, abbreviations, template
    ├── schema/             Process scheme diagrams
    ├── presentation/       Defence presentation
    └── images/             Figures / photos
```

## Model notes

- **NRTL** is used for activity coefficients; interaction parameters are given as
  `delta_g` [kJ/mol] with non-randomness factor `alpha`.
- The **LLE** equilibrium is solved by successive substitution on the distribution
  coefficients `K = gamma_R / gamma_E`, with damping for stability
  (see [matlab/LLE_solver.m](matlab/LLE_solver.m)). The improved Python variant
  solves the multistage extractor simultaneously with a Newton/TRF solver.
- **VLE** uses a Rachford–Rice flash with Antoine `Psat`; the ionic liquid is
  treated as non-volatile. The improved variant uses Brent's method for the
  low-pressure bubble point (≈2.45 kPa).

## Note on sources

This repository contains only the author's own work (simulations, code, thesis
documents). External literature — published papers, textbooks and third-party
theses used as references during the project — is intentionally **not** included
here for copyright reasons. They are cited in the thesis bibliography
([docs/thesis/Literatura.docx](docs/thesis/Literatura.docx)).
