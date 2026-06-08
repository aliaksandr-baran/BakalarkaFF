# Separation of Azeotropic Mixtures using Ionic Liquids

Bachelor's thesis project (FCHPT) by **Aliaksandr Baran**.

The work models the separation of azeotropic mixtures (e.g. MTBE–methanol) using
ionic liquids as entrainers. It covers vapour–liquid (VLE) and liquid–liquid (LLE)
equilibria described with the **NRTL** activity-coefficient model, together with the
material/energy balances and economics of the proposed separation scheme
(extractor, evaporator, valve, IL regeneration).

## Goal: migrating the simulation Excel → MATLAB → Python

The calculations were first built in **Excel** (current, complete implementation).
The plan is to port the simulation step by step:

```
Excel  ──►  MATLAB  ──►  Python
(current)   (in progress)  (planned)
```

1. **Excel** — the original spreadsheets with all VLE/LLE/balance/economics
   calculations. This is the reference implementation.
2. **MATLAB** — re-implementation of the equilibrium solvers (NRTL, LLE
   successive-substitution) as scripts/functions. Started.
3. **Python** — final goal: a clean, scriptable port of the model. Not started yet.

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
│   ├── bacalar_project.m   Main LLE equilibrium script (NRTL, ternary)
│   └── LLE_solver.m        LLE_solver() function + NRTL gamma calculation
├── python/                Planned Python port (placeholder)
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
  (see [matlab/LLE_solver.m](matlab/LLE_solver.m)).

## Note on sources

This repository contains only the author's own work (simulations, code, thesis
documents). External literature — published papers, textbooks and third-party
theses used as references during the project — is intentionally **not** included
here for copyright reasons. They are cited in the thesis bibliography
([docs/thesis/Literatura.docx](docs/thesis/Literatura.docx)).
