# modalNEB

Automated Python/Matlab Toolkit for Modal NEB Analysis and Visualization  
Author: Kiarash Gordiz  
Contact: [add your email or link, optional]

---

## Overview

This repository contains tools for running, analyzing, and visualizing **modal NEB (Nudged Elastic Band)** calculations, especially for studying ion migration and phonon mode coupling in solid-state materials.  
It provides a robust workflow that uses both Python and MATLAB for:

- Automated simulation setup and data generation (Python)
- Advanced post-processing and scientific visualization (MATLAB)

The theoretical framework and results are published in:

- [A modal analysis of lithium diffusion in solid electrolytes](https://www.sciencedirect.com/science/article/pii/S2666386421001260), *Cell Reports Physical Science* (2021)
- [A general atomistic framework for computing modal contributions to migration barriers](https://arxiv.org/abs/2305.01632), *arXiv preprint* (2023)

---

## Features

- **Batch Processing:** Automated analysis across multiple simulation folders.
- **Flexible Analysis:** Extracts, computes, and compares atomistic projections, energy contributions, and mode participation.
- **High-Quality Visualization:** Generates publication-quality figures, using MATLAB’s advanced plotting.
- **Hybrid Workflow:** Python for data generation and simulation, MATLAB for post-processing and visualization.
- **Fully Scriptable:** Easily adaptable for different systems and simulation types.

---

## Workflow

1. **Simulation Setup & Execution (Python)**
   - Follow instructions in the main package and scripts to generate output folders (see `/modalNEB/`).

2. **Post-Processing & Plotting (MATLAB)**
   - Navigate to `/analysis/matlab/`
   - Place all simulation output folders (from Python runs) in the working directory or update the script’s `foldernames` as needed.
   - Run the provided `.m` script to generate analysis and figures.

3. **Result Interpretation**
   - Analyze the outputs for trends in ion migration, the role of different phonon modes, and elemental contributions.
   - Visualize frequency-resolved and atom-resolved contributions to migration barriers.

4. **Example Results**
   - *Optional: Add 1–2 plot images or a link to an example figure here.*

---

## Directory Structure

modalNEB/
|-- modalNEB/ # Main Python package for running modal NEB calculations
|-- scripts/ # Example scripts and automation tools
|-- analysis/
|-- matlab/ # MATLAB code for post-processing and visualization
|-- README.md
|-- setup.py
|-- requirements.txt


---

## Publications

- **A modal analysis of lithium diffusion in solid electrolytes**  
  [Cell Reports Physical Science, 2021](https://www.sciencedirect.com/science/article/pii/S2666386421001260)
- **A general atomistic framework for computing modal contributions to migration barriers**  
  [arXiv preprint, 2023](https://arxiv.org/abs/2305.01632)

---

## License

[Add license info here if applicable.]

---

## Contact

For questions, suggestions, or collaboration inquiries, please contact:  
Kiarash Gordiz ([add email/link])

---

## Acknowledgments

Development supported by [your lab/PI/funding if you want to mention].

---

**Instructions for contributors:**  
- To add images, use the Markdown syntax:  
  `![Description](relative/path/to/figure.png)`  
  *(You can upload an image to the repo and reference it this way.)*

---
