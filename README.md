# JAK1-drug-discovery-pipeline
Open-source, automated pipeline for JAK1 drug discovery—includes data acquisition/curation, QSAR modeling, molecular docking, de novo generation, and ADMET prediction using Python and open tools such as RDKit, AutoDock Vina, and DeepChem. Fully reproducible and configurable.

JAK1 Drug Discovery Pipeline


This repository provides an open-source, automated pipeline for end-to-end drug discovery targeting JAK1.


Project Overview


A modular Python workflow for discovering and optimizing drug candidates against JAK1. Designed for full reproducibility, scalability, and ease of use.


Workflow Stages


Data Acquisition & Curation: Download, clean, and prepare datasets relevant for JAK1.
QSAR Modeling: Calculate descriptors and build predictive models for compound bioactivity.
Molecular Docking: Virtual screening of molecules against JAK1 using open-source docking tools.
De Novo/Generative Design: Generate and optimize new candidate molecules based on docking and QSAR results.
ADMET Prediction: Assess absorption, distribution, metabolism, excretion, and toxicity with open-source ADMET tools.

Technologies Used


Python (main language)
RDKit, AutoDock Vina, DeepChem, REINVENT, Pandas, scikit-learn
Tools for ADMET prediction and visualization
Getting Started


Clone the repository
Install dependencies from `requirements.txt` or `environment.yml`
Edit `config.yaml` to adjust parameters for each pipeline stage
Run the workflow via `main.py` or individual stage scripts

Repository Structure


`/data_preprocessing/` — Data acquisition and cleaning
`/qsar_modeling/` — Feature engineering and QSAR models
`/docking/` — Docking automation and analysis
`/generative_design/` — De novo molecule generation
`/admet/` — ADMET prediction scripts
`/notebooks/` — Jupyter exploratory notebooks
`/utils/` — Shared helper functions and modules

Contribution


PRs and issues welcome! Please see CONTRIBUTING.md for guidelines (to be added).
Contact: Jake Aquilina, jaquilinag@gmail.com
