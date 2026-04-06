# NTL — Non-uniform Transmission Line Design Toolkit

A C++20 library and suite of applications for designing, optimising, and simulating non-uniform transmission lines (NTLs) and microwave power-splitting networks built from them: T-junctions (TJ) and Wilkinson power dividers (WPD).

---

## Overview

Transmission line impedance tapering is a core technique in microwave engineering. This toolkit models a transmission line whose characteristic impedance varies continuously along its length as a Fourier series of cosine and sine terms:

```
Z(z) = Z₀ · exp( Σ Cₙ cos(2πnz/d) + Σ Cₙ sin(2πmz/d) )
```

The Fourier coefficients (`Cn`) are the free parameters. The core solver propagates signals through the line by chaining 2×2 ABCD (transmission) matrices evaluated at discrete segments, then converts to S-parameters or Y-parameters as needed. Optimisation finds the `Cn` values that minimise input reflection across a set of target frequencies subject to impedance bounds.

---

## Repository layout

```
NTL_core/               Static library — all physics, models, and optimisers
  src/
    models/
      ntl.h / ntl.cpp         NTL model and matrix calculations
      tj.h  / tj.cpp          T-Junction (3-port) model
      wpd.h / wpd.cpp         Wilkinson Power Divider (3-port) model
      flt.h / flt.cpp         Filter model (stub)
    optimisation/
      optimiser.h / .cpp      Abstract NLopt-based optimiser base class
      NTL_opt.h / .cpp        Single NTL optimiser
      TJ_opt.h  / .cpp        T-Junction optimiser
      WPD_opt.h / .cpp        Wilkinson Power Divider optimiser
      flt_opt.h / .cpp        Filter optimiser (stub)
    common/
      file_handler.h / .cpp   JSON serialisation (nlohmann::json)
      helpers.h               Random start, ostream helpers
      enums.h                 phase, mag enums
    validation.h              Structural sanity checks on results
    verification.h            S-parameter performance checks on results

NTL_sim/                Static library — Qt-based simulation and plotting
  src/simulation/
    NTL_sim.h / .cpp    Frequency-sweep simulation for NTLs
    TJ_sim.h  / .cpp    Frequency-sweep simulation for T-Junctions
    WPD_sim.h / .cpp    Frequency-sweep simulation for WPDs

NTL_optimiser/          Executable — optimise a single NTL from a JSON setup file
NTL_simulator/          Executable — simulate a saved NTL
TJ_optimiser/           Executable — optimise a T-Junction from a JSON setup file
TJ_simulator/           Executable — simulate a saved T-Junction
WPD_optimiser/          Executable — optimise a Wilkinson divider from a JSON setup file
WPD_simulator/          Executable — simulate a saved Wilkinson divider
NTL_AutoCAD/            Executable — export NTL geometry as AutoCAD .scr scripts
NTL_data_extract/       Executable — export impedance profiles and S-params to CSV
```

---

## Dependencies

| Dependency | Role |
|---|---|
| [Eigen 3](https://eigen.tuxfamily.org) | Dense complex matrix arithmetic |
| [NLopt](https://nlopt.readthedocs.io) | Global (GN_ISRES) and local (LN_COBYLA / LD_MMA) optimisation |
| [nlohmann/json](https://github.com/nlohmann/json) | JSON serialisation for setup files and results |
| [Qt 6](https://www.qt.io) (core, gui, widgets, charts) | Plotting and windowed simulation output |
| OpenMP | Parallel frequency-loop evaluation in the objective function |

All projects target **x64 Windows**, toolset **v145** (MSVC), C++20 standard. Property sheets in `my_props/` centralise include paths and library links.

---

## Core concepts

### NTL model

An `NTL::NTL` object holds:

- `Z0` — reference characteristic impedance (Ω)
- `er` — substrate relative permittivity
- `d` — physical length (m)
- `Cn` — Fourier coefficient vector (length `N`; the last `M` entries are sine terms)
- `M` — number of sine terms (0 = symmetric / even profile)

Key methods:

```cpp
ntl.Z(z)                          // impedance profile at position z
ntl.T_matrix(f, K)                // 2×2 ABCD matrix at frequency f with K segments
ntl.S_matrix(f, Zs, Zl, K)       // 2×2 S-matrix
ntl.Y_matrix(f, K)                // 2×2 Y-matrix
ntl.Zin(Zl, f, K)                 // input impedance
ntl.electrical_length(f, K)       // total electrical length (radians)
ntl.W_H(z)                        // microstrip width-to-height ratio at z
```

Effective permittivity dispersion is computed from the Hammerstad–Jensen closed-form microstrip model.

### T-Junction (TJ)

A `TJ::TJ` holds two NTL arms (arm 2 and arm 3) sharing a common input port (port 1). The 3×3 Y-matrix is assembled by stamping the two 2×2 NTL Y-matrices into a shared node, then converted to power-wave S-parameters via the Kurokawa formulation.

### Wilkinson Power Divider (WPD)

A `WPD::WPD` extends the TJ by adding an isolation resistor `R` between the two output ports. The 3×3 Y-matrix includes the resistor conductance term. A `system_S_matrix` method further attaches output matching transformers (additional NTLs) and performs a Schur complement reduction to present a clean 3-port from the external reference planes.

---

## Optimisation workflow

### NTL optimiser

Minimises the sum of squared reflection coefficients `|Γ(fᵢ)|²` across all target frequencies by finding the Fourier coefficients `Cn`.

Constraints:
- Equality: `Z(0) = Z_at_0`, `Z(d) = Z_at_d`
- Inequality: `Z_min ≤ Z(z) ≤ Z_max` at K sample points

Two-phase search per attempt:
1. **Global phase** — `GN_ISRES` stochastic search with `GBL_MAX` evaluations
2. **Local phase** — `LN_COBYLA` (gradient-free, default) or `LD_MMA` via `LD_AUGLAG` for gradient-based refinement

The impedance and effective permittivity profiles are pre-computed once per `Cn` evaluation (frequency-independent), then passed to a vectorised T-matrix function, with the inner frequency loop parallelised with OpenMP.

An `optimise_d` mode additionally trims the line length by `1 mm` steps until the accepted error threshold is no longer met, giving the shortest viable line.

### TJ optimiser

Decomposes the problem: each arm is independently optimised as an NTL impedance transformer using `NTL::opt`. Source and load impedances are derived from the desired power split ratio:

- Arm 2: `Zs = Zref · (1 + split)`, `Zl = Zref`
- Arm 3: `Zs = Zref · (1 + 1/split)`, `Zl = Zref`

### WPD optimiser

Three sequential stages:

1. **Arms** — optimise arm 2 and arm 3 jointly (2N coefficients) targeting conjugate match at each arm's source impedance, with an additional cost term penalising imbalance in the A-element of the T-matrix.
2. **Output transformers** — two independent `NTL::opt` runs to match `Zref·√split` and `Zref/√split` back to `Zref`.
3. **Isolation resistor** — golden-section search on `R ∈ [R_min, R_max]` minimising the reflection seen from the resistor looking back into both output nodes.

---

## Setup files

Each optimiser reads a JSON file at startup. Example fields for `tj_setup.json`:

```json
{
  "json_type": "setup",
  "setup_type": "TJ_opt",
  "N": 6,
  "M": 0,
  "GBL_MAX": 5000,
  "LCL_MAX": 10000,
  "accepted_error": 1e-4,
  "max_attempts": 10,
  "Z0": 50.0,
  "er": 4.4,
  "d": 0.03,
  "K": 50,
  "Zref": 50.0,
  "Z_min": 20.0,
  "Z_max": 120.0,
  "freqs": [1e9, 2e9, 3e9],
  "split": [1.0]
}
```

`N` is the total number of Fourier coefficients; `M` is the number of sine terms (must satisfy `M ≤ N`). Setting `M = 0` enforces a symmetric impedance profile and halves the T-matrix computation via even-symmetry exploitation.

---

## Simulation and output

All simulator executables produce Qt chart windows showing:

- **Impedance profile** `Z(z)` along the line length
- **Width-to-height profile** `W/H(z)` (microstrip geometry)
- **Input impedance** magnitude and phase vs frequency
- **S-parameters** (matching, split, isolation) in dB vs frequency

The `NTL_data_extract` tool exports the same data to CSV files compatible with HFSS and CST import formats, and writes coefficient tables suitable for inclusion in reports.

The `NTL_AutoCAD` tool exports the microstrip conductor outline as AutoCAD `.scr` scripts for both straight segments and bent arcs at arbitrary angles.

---

## Saving and loading results

Results are serialised to JSON:

```cpp
fh::ntl_to_file(ntl, "folder/name_ntl");
fh::tj_to_file(tj,   "folder/name_tj");
fh::wpd_to_file(wpd, "folder/name_wpd");
fh::setup_to_file<TJ::opt_setup>(setup, "folder/name_setup");
```

and reloaded by the simulator executables:

```cpp
NTL::NTL ntl = fh::file_to_ntl("folder/name_ntl");
TJ::TJ   tj  = fh::file_to_tj("folder/name_tj");
```

---

## Build

Open `NTL.slnx` in Visual Studio 2022. Build order is managed via solution dependencies:

```
NTL_core  →  NTL_sim  →  *_optimiser / *_simulator executables
```

Property sheets in `../../my_props/` must resolve correctly for include paths and library directories. Adjust `dir_props.props`, `optimisation_release.props`, `Link_NTL_core.props`, `Link_NTL_sim.props`, and `qt_plot.props` to match your local installation paths for Eigen, NLopt, nlohmann/json, and Qt.

Release builds enable OpenMP (`/openmp`) and whole-program optimisation (`/GL`). The TJ and WPD optimisers are single-threaded at the outer loop but parallelise the inner frequency evaluation across up to `min(F, hardware_concurrency)` threads.
