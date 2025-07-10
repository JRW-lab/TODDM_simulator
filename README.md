# TODDM Simulator (MATLAB, SQL)
A MATLAB-based simulator for calculating error rates in OTFS, ODDM and TODDM.
Results can be stored either locally in an Excel file or in an SQL database.

## Introduction
To use this code, you must run it in MATLAB 2024b or higher. The parallelization toolbox is used in the current implementation but can be turned "off" in settings. Additionally, the database toolbox and several others are required, both to simulate and to upload simulation results to MySQL. SQL files are included to run in order to automatially create the needed tables for MATLAB to read/write from.

## Instructions
The code included here is lengthy and may be confusing so here is an overview of how it works:

1. MAIN_simulator_v2.m includes the configurations and when run, the user selects from a series of options.
2. If a sufficient number of frames is not already simulated, sim_save.m is run for a specific system with a set of defined parameters
3. Based on the system_name in parameters, a simulation file is selected and additional frames are run.
4. Steps 2 and 3 are repeated until all configurations have the sufficient number of frames for figure rendering.
5. gen_figure_v2 or gen_hex_layout.m is run to generate a figure and save, if specified.

NOTE: At of time of writing, gen_hex_layout.m does not support loading local simulation results.

## Configuration Setup
In MAIN_simulator_v2.m, there is a large section named Configurations. There you will see several pre-configured profiles to use as reference for your own profile. Profiles work by defining the primary variable for a parametric sweep and the corresponding range. If a figure is being rendered, this is the range of the plot, and each line of the parameter 'configs' specifies a line on the plot, and each line has its own custom parameters separate from those specified in default_parameters. Once all the configs are defined, the user can be specific in defining the appearance of plots using several customizable parameters.

Hexgrids offer a more robust way of viewing tradeoffs being variables with TODDM (see paper this code is a reference for). Hexgrid plots function a bit differently, where there is no parametric sweep and 'configs' is defined by a special function. There also aren't any visualization settings, as a simple colorbar is all you can really change and the current one feels more fitting and easiest to read.

## Further Questions
For any questions please contact jrwimer@uark.edu.
