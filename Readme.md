# Calibrating earthquake impacts with satellite-mapped damage (code & demos)

**Repository:** `eq-impact-pipeline-with-calibration`  

This repository contains MATLAB code supporting the paper:

**_Calibrating earthquake impacts with satellite-mapped damage revises national risk hotspots_**

---

## Summary

Credible earthquake-impact estimation is essential for seismic risk reduction, yet remains difficult where models are weakly calibrated and rarely validated. We develop an end-to-end, modular earthquake-impact estimation pipeline that combines globally available hazard, exposure and fragility datasets to estimate fatalities, direct building economic losses and post-earthquake hospital accessibility. The pipeline explicitly represents shaking-triggered cascades (liquefaction, landslides and building-to-road debris blockage) and uses satellite-mapped building collapses and road closures to calibrate a small set of asset-level correction parameters that reconcile biases across modules and their integration. Calibrated to satellite-mapped damage in 25 areas of interest from the 2025 Mw 7.7 Myanmar earthquake, the pipeline reproduces national impacts for this and earlier events and reshapes national seismic risk maps.

---

## What’s included

- **Two demos** (MATLAB scripts) showing how to estimate:
  - fatalities, direct economic losses, and injuries  
    - `demo_loss_fatality.m`
  - injuries with post-earthquake travel time to hospital exceeding specified thresholds  
    - `demo_accessibility.m`
- **Supporting demo data** in `data/`
- **Module code**:
  - `asset_damage/` — asset damage modelling components
  - `socio_economic_impact/` — socio-economic impact components
  - `utils/` — shared helper functions/utilities
- **Figure reproduction scripts (Figures 2–5)** in `figures/`
  - `figure_2.m`, `figure_3.m`, `figure_4.m`, `figure_5.m`

Other related code and/or resources can be provided upon request.  
Contact: **chongyang_du@hust.edu.cn**

---


## Repository structure

.
├── asset_damage/ # Asset damage modelling
├── socio_economic_impact/ # Socio-economic impact modelling
├── utils/ # Shared utilities
├── data/ # Data required to run the demos
├── figures/ # Scripts to reproduce Figures 2–5
│ ├── results/ # Output directory for generated results
│ ├── figure_2.m
│ ├── figure_3.m
│ ├── figure_4.m
│ ├── figure_5.m
│ └── Figure 4c and 4e.xlsx
├── demo_loss_fatality.m # Demo: fatalities, losses, injuries
└── demo_accessibility.m # Demo: injuries beyond travel-time thresholds

---

## Requirements

- **MATLAB** (version R2025a)

---

## Quickstart

1. Clone the repository and open MATLAB.
2. Set MATLAB’s current folder to the repository root.

---

## Run the demos

### 1) Fatalities, losses, and injuries

demo_loss_fatality.m

### 2) Accessibility / injuries beyond travel-time thresholds

demo_accessibility.m

---

## Reproduce Figures 2–5

All scripts used to generate Figures 2–5 are located in figures/.

From the repository root:

figures/figure_2.m
figures/figure_3.m
figures/figure_4.m
figures/figure_5.m

---

## Contact

For questions, issues, or requests for additional code/resources:
chongyang_du@hust.edu.cn