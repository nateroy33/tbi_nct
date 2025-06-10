# TBI Network Control Analysis

This repository contains MATLAB scripts used in our network control theory analysis of traumatic brain injury (TBI) data, now published in NeuroImage: Clinical (https://pubmed.ncbi.nlm.nih.gov/40381376/)

The code implements and adapts methods from:
- [Eli Cornblath’s brain_states repository](https://github.com/ejcorn/brain_states) ([Cornblath et al., *Communications Biology*, 2020])
- [Parker Singleton’s energy_landscape repository](https://github.com/singlesp/energy_landscape). ([Singleton et al., *Nature Communications*, 2022])

Some functions from these repositories have been modified and consolidated for this project.

### Contents
The scripts here are organized to replicate the full pipeline used in the corresponding manuscript, including:
- Data preparation and input structuring
- Repeated k-means clustering and control energy calculations
- Entropy and transition analysis
- Visualization scripts (violin plots, energy matrices, etc.)

Each script was designed for a specific step in our analysis workflow and may not be fully generalized or cleaned.

### Notes
This code is shared for transparency and reproducibility. For questions or clarification, contact [Nate Roy] at nr284@cornell.edu.

