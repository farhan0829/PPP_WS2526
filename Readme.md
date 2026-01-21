<img src="https://r2cdn.perplexity.ai/pplx-full-logo-primary-dark%402x.png" style="height:64px;margin-right:32px"/>

# README.md - Reaction-Diffusion Simulator Test Cases \& Examples

## Test Case Files Mapping

| **Test Case** | **Generated Files** | **Purpose** | **Location** |
| :-- | :-- | :-- | :-- |
| **Test Case 1: Heat Equation** | `Heat.mat`, `Heat.png` | Pure diffusion verification (Neumann BC) | `data/` |
| **Test Case 2: Turing Spots** | `Turing_spots.mat`, `Turing_spots.png` | Schnakenberg spots ($D_v=100$) | `data/` |
| **Test Case 3: Stripes** | `Stripes.mat`, `Stripes.png` | Schnakenberg stripes ($D_v=30$) | `data/` |
| **Test Case 4: Brusselator** | `Brusselator.mat`, `Brusselator.png` | Brusselator Turing pattern ($A=1.945,B=2.88$) | `data/` |
| **Test Case 5: Console Runner** | `RD_console-schnakenberg.jpg`<br>`RD_console-brusselator.jpg`<br>`RD_console-heat.jpg` | Non-GUI validation (3 models) | `screenshots/` or `data/` |

## Reference Images (Validation)

| **File** | **Test Case** | **Source** |
| :-- | :-- | :-- |
| `Heat_visualpde.png` | Test Case 1 | VisualPDE.com reference |
| `turing_spots_visualpde.png` | Test Case 2 | VisualPDE.com reference |
| `Stripes_visualpde.png` | Test Case 3 | VisualPDE.com reference |
| `Brusselator_visualpde.png` | Test Case 4 | VisualPDE.com reference |

## Core Files Structure

```
project/
├── run_rd.m                 # Main launcher (GUI + console)
├── run_rd_console.m         # Console-only runner (Test Case 5)
├── SimulationGUI.m          # GUI implementation
├── OOPDELight/              # FEM library (Bilinear2D.m, RectangleR.m)
├── data/                    # Test outputs (.mat, .png)
├── screenshots/             # Console screenshots
└── report.pdf               # Validation documentation
```


## Quick Test Reproduction

```matlab
% Test Case 5 (Console - all 3 models)
run_rd_console('schnakenberg', 0.1, 0.9, 1, 8)
run_rd_console('brusselator',  1.0, 3.0, 1, 8) 
run_rd_console('heat',         0,   0,   1, 1)
```

**All `.mat` files contain full simulation history** (`u,v,x,y,t,parameters`) for post-processing and reproducibility verification.
<span style="display:none">[^1][^10][^2][^3][^4][^5][^6][^7][^8][^9]</span>

<div align="center">⁂</div>

[^1]: paste.txt

[^2]: image.jpg

[^3]: image.jpg

[^4]: paste.txt

[^5]: Screenshot-2026-01-21-at-18.14.07.jpg

[^6]: paste.txt

[^7]: RD_console-brusselator.jpg

[^8]: RD_console-heat.jpg

[^9]: RD_console-schnakenberg.jpg

[^10]: paste.txt

