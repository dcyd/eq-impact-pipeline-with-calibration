# Calibrating earthquake impacts with satellite-mapped damage (code & demos)

**Repository:** `eq-impact-pipeline-with-calibration`  
This repository contains MATLAB code supporting the paper:

**_Calibrating earthquake impacts with satellite-mapped damage revises national risk hotspots_**

---

## Summary

Credible earthquake-impact estimation is essential for seismic risk reduction but is often difficult where models are weakly calibrated and rarely validated. This repository provides an end-to-end, modular MATLAB pipeline that incorporates satellite-mapped damage to improve calibration of impact estimates (fatalities, economic loss, injuries) and reproduce figures from the paper.

---

## What’s included

- Two demo scripts showing example workflows:
  - [demo_loss_fatality.m](https://github.com/dcyd/eq-impact-pipeline-with-calibration/blob/main/demo_loss_fatality.m) — estimate fatalities, direct economic losses, and injuries
  - [demo_accessibility.m](https://github.com/dcyd/eq-impact-pipeline-with-calibration/blob/main/demo_accessibility.m) — compute injuries based on post-earthquake travel time to hospitals exceeding thresholds
- Supporting demo data: [data/](https://github.com/dcyd/eq-impact-pipeline-with-calibration/tree/main/data)
- Module code:
  - [asset_damage/](https://github.com/dcyd/eq-impact-pipeline-with-calibration/tree/main/asset_damage) — asset damage modelling components
  - [socio_economic_impact/](https://github.com/dcyd/eq-impact-pipeline-with-calibration/tree/main/socio_economic_impact) — socio-economic impact components
  - [utils/](https://github.com/dcyd/eq-impact-pipeline-with-calibration/tree/main/utils) — shared helper functions/utilities
- Figure reproduction scripts (Figures 2–5): [figures/](https://github.com/dcyd/eq-impact-pipeline-with-calibration/tree/main/figures)
  - [figures/figure_2.m](https://github.com/dcyd/eq-impact-pipeline-with-calibration/blob/main/figures/figure_2.m)
  - [figures/figure_3.m](https://github.com/dcyd/eq-impact-pipeline-with-calibration/blob/main/figures/figure_3.m)
  - [figures/figure_4.m](https://github.com/dcyd/eq-impact-pipeline-with-calibration/blob/main/figures/figure_4.m)
  - [figures/figure_5.m](https://github.com/dcyd/eq-impact-pipeline-with-calibration/blob/main/figures/figure_5.m)

Other related code and/or resources can be provided upon request. Contact: **chongyang_du@hust.edu.cn**

---

## Repository structure (high level)

.
- [asset_damage/](https://github.com/dcyd/eq-impact-pipeline-with-calibration/tree/main/asset_damage) — Asset damage modelling
- [socio_economic_impact/](https://github.com/dcyd/eq-impact-pipeline-with-calibration/tree/main/socio_economic_impact) — Socio-economic impact modelling
- [utils/](https://github.com/dcyd/eq-impact-pipeline-with-calibration/tree/main/utils) — Shared utilities
- [data/](https://github.com/dcyd/eq-impact-pipeline-with-calibration/tree/main/data) — Data required to run the demos
- [figures/](https://github.com/dcyd/eq-impact-pipeline-with-calibration/tree/main/figures) — Scripts to reproduce Figures 2–5
- [demo_loss_fatality.m](https://github.com/dcyd/eq-impact-pipeline-with-calibration/blob/main/demo_loss_fatality.m)
- [demo_accessibility.m](https://github.com/dcyd/eq-impact-pipeline-with-calibration/blob/main/demo_accessibility.m)

---

## Requirements

- MATLAB (tested with R2025a)

---

## Quickstart

1. Clone the repository and open MATLAB:
   - git clone https://github.com/dcyd/eq-impact-pipeline-with-calibration.git
2. In MATLAB, set the Current Folder to the repository root.
3. Add repo folders to the MATLAB path if you prefer:
   - addpath(genpath(pwd))

---

## Run the demos

From the repository root in MATLAB:

1. Fatalities, losses, and injuries:
   - open and run: [demo_loss_fatality.m](https://github.com/dcyd/eq-impact-pipeline-with-calibration/blob/main/demo_loss_fatality.m)
   - In MATLAB command window: run('demo_loss_fatality.m') or just type `demo_loss_fatality`

2. Accessibility / injuries beyond travel-time thresholds:
   - open and run: [demo_accessibility.m](https://github.com/dcyd/eq-impact-pipeline-with-calibration/blob/main/demo_accessibility.m)
   - In MATLAB command window: run('demo_accessibility.m') or type `demo_accessibility`

---

## Reproduce Figures 2–5

All figure scripts are in [figures/](https://github.com/dcyd/eq-impact-pipeline-with-calibration/tree/main/figures). From the repository root in MATLAB:

- run('figures/figure_2.m')
- run('figures/figure_3.m')
- run('figures/figure_4.m')
- run('figures/figure_5.m')

Figures will write outputs to [figures/results/](https://github.com/dcyd/eq-impact-pipeline-with-calibration/tree/main/figures/results).

---

## Contact

For questions, issues, or requests for additional code/resources:
chongyang_du@hust.edu.cn
