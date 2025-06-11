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

Each script was designed for a specific step in our analysis workflow and may not be fully generalized or cleaned:

[**repeatkmeans.m**](repeatkmeans.m) and elbow.m were used to find and visualize the optimal number of clusters for explaining variance in the data.

[**ami_calc.m**](ami_calc.m) was used to assess clustering stability and to choose the partition with the highest amount of adjusted mutual information shared with all other partitions.

[**systems_plot.m**](systems_plot.m) and [**plotcentroidsSPSnate.m**](plotcentroidsSPSnate.m) (via https://github.com/kjamison/atlasblobs) were used to visualize brain states.

[**comb_clusters.m**](comb_clusters.m) was used to reorder brain states for interpretative clarity.

[**transprobs5gDS1.m**](transprobs5gDS1.m) was used to calculate and visualize state transition probabilities.

[**countclusters5g.m**](countclusters5g.m) was used to calulcate state fractional occupancy.

[**violin5G.m**](violin5g.m) was used to visualize data from countclusters.m.

[**TEtbi4F.m**](TEtbi4F.m) was used for all energetic calculations and visualization.

[**entropmild.m**](entropmild.m) was used for entropy calculations and visualization.

### Notes
This code is shared for transparency and reproducibility. For questions or clarification, contact Nate Roy at nr284@cornell.edu.

