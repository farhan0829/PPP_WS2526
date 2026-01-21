# Reaction-Diffusion Simulator (MATLAB)

MATLAB reaction–diffusion simulator for **Schnakenberg**, **Brusselator**, and **Heat** (pure diffusion) models.  
Includes a full **App Designer GUI**, a **console runner**, and a **scripted runner**, all built on an included **FEM (OOPDELight)** toolbox.

This repo also contains example outputs (MAT/PNG/AVI), VisualPDE reference images, screenshots, and a packaged MATLAB Toolbox release.

---

## Features

- GUI (App Designer) with real-time controls and visualization
- Console runner for quick validation (no GUI)
- Script runner for reproducible batch runs
- FEM backend (Bilinear2D + RectangleR) bundled in `OOOPDELight/`
- Export data (MAT/CSV) and animations (AVI/MP4) from the GUI
- Preset pattern modes (Turing spots, stripes) and boundary conditions

---

## Requirements

- **MATLAB** (App Designer required only for `.mlapp`)
- No external dependencies beyond the included **OOPDELight** library

If you can’t open `.mlapp`, use the generated class file `Simulator_GUI.m` instead.

---

## Quick Start

### 1) Add repo to MATLAB path
```matlab
addpath(genpath(pwd))
```

### 2) Launch GUI
```matlab
run_rd("gui")
```

### 3) Launch Console Runner (interactive)
```matlab
run_rd("console")
```

---

## Console Usage (non-GUI)

Direct calls:
```matlab
run_rd_console('schnakenberg', 0.1, 0.9, 1, 8)
run_rd_console('brusselator',  1.0, 3.0, 1, 8)
run_rd_console('heat',         0,   0,   1, 1)
```

---

## Scripted Runs

Edit and run:
```matlab
run_rd_script.m
```

Key inputs inside the script:

- `modelName` = `"schnakenberg" | "brusselator" | "heat"`
- `a`, `b`, `Du`, `Dv`, `dt`, `tEnd`, mesh size `h`

---

## GUI Notes

The GUI is implemented in:

- `Simulator_GUI.mlapp` (App Designer project)
- `Simulator_GUI.m` (exported class)

Key features in the GUI:

- Model selection (Schnakenberg / Brusselator / Heat)
- Boundary conditions (Periodic / Dirichlet / Neumann)
- Plot types (2D heatmap, 3D surface, contour)
- History recording, data export, and animation export
- Preset pattern buttons (Turing spots / stripes)

---

## Data & Outputs Included

### Example Outputs (simulation results)
- `Heat.mat`, `Heat.png`, `Heat.avi`
- `Brusselator.mat`, `Brusselator.png`, `Brusselator.avi`
- `Stripes.mat`, `Stripes.png`, `Stripes.avi`
- `Turing_Spots.mat`, `Turing_Spots.png`, `Turing_Spots.avi`
- `rd_animation.avi`

### VisualPDE Reference Images
- `Heat_visualpde.png`
- `Stripes_visualpde.png`
- `Brusselator_visualpde.png`
- `turing_spots_visualpde.png`

### Console Screenshots
- `RD_console:heat.png`
- `RD_console:schnakenberg.png`
- `RD_console:brusselator.png`

### UI / Misc Images
- `Matlab_Appmaker.png`
- `Screenshot.png`

---

## Main Files

- `run_rd.m` — main entry point (GUI or console)
- `run_rd_console.m` — console runner
- `run_rd_script.m` — scripted runner
- `Simulator_GUI.mlapp` — GUI project
- `Simulator_GUI.m` — GUI class
- `setModel.m` — default model parameters
- `FEM.m` — FEM solver helper
- `OOOPDELight/` — FEM toolbox (Bilinear2D, RectangleR, mesh utilities)

---

## Utilities

- `Convert.m` — extracts `params` from all `.mat` files and writes to `output_params.xlsx`
- `testSymmetryBC.m` — periodic boundary condition demo (OOPDELight)
- `solverScript.m` — legacy example script
- `simulator.m` — older console/GUI wrapper (legacy)

---

## Releases / Build Artifacts

- `release/Farhan_app_2.mltbx` — packaged MATLAB Toolbox
- `release/build/` — compiled artifacts + logs

---

## References

- `OOOPDELight/testOOPDELight.pdf` — FEM toolbox tests

---

## License

MIT — see `LICENSE`
