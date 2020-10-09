# Procedure used in the article of Zanatta et al 2020 in Nature Communications
This repository contains the source code and scripts that were used to run the dispersal
simulations presented in Zanatta et al 2020.

The R scripts prefixed `01_` to `04_` allow to run the dispersal simulations for all species. They
are intended to be run in their prefix's number order.  
The file `migclim.tar.gz` contains a custom version of the migclim R package used in the project.
The code must be compiled into an actual R package before it can be used.

### Brief description of each R script purpose:
* **01_SDMs_and_kernel_map_computation.R**: script to generate the habitat suitability maps for the
  different Bryophytes species under current and future climatic conditions. The habitat suitability
  maps are then used as inputs for the dispersal simulations.  
* **02_setup_input_data.R:** script that prepares input data for running the dispersal simulations
  across the multiple species and multiple climate change, wind speed, release height and
  long-distance dispersal scenarios.
* **03_run_migclim.R:** script to run the actual dispersal simulations.
* **04_merge_outputs.R:** script to aggregate all simulation outputs across different species and
  replicates into a single table.
* **05_plot_species_projections.R:** script that produces plots of species distributions under
  current and projected future climatic conditions. This script is not strictly necessary to
  produce the simulation outputs.
