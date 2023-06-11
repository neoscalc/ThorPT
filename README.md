# ThorPT - This is ThorPT
ThorPT is a modelling package for fluid coupled geolgical processes in the crust. Fluid production and migration in the crust is fundamental for the earth's geochemical cycling and ThorPT is specified for modelling this fluid production and fluid transfer. It includes the latest developments in petrogeochemical modelling and helps the user to test different scenarios:
- Single rock modelling:
    1. Petrological purposes to test changes in mineral assemblages and coupled fluid production, e.g., while prograde metamorphism
    2. Mechanical testing of the rock based on the Mohr-Coulomb theory
    3. Quantification of the production of fluid masses

- Multi rock modelling (*Work in progress):
    1. Stacks of rock sequences involving fluid transfer
    2. Outcrop modelling
    3. *Contact zones and fluid-rock ratios

# 1. Getting started and use the package
Information and prerequisits:
There are two ways to use the package at the moment. Prerequisit is in both cases to get and install a working theriak-domino version from
https://github.com/Theriak-Domino/theriak-domino

Getting started:

1.1 (Recommended for all users) Quick start: Get the jupyter working directory that includes all instructions and direct connection to ThorPT.

1.2 (Recommended to users familair with coding) Install the package manually into your python environment and develop your own script.

# 1.1 Use the package with the jupyter environment
Download the jupyter working directory that includes:
- Jupyter script
- Datafiles directory with
    - initial file
    - examples to model

Before starting open the "_init_.txt" file in the working directory and in the first line behind "Theriak:" write the path of the directory of the working "theriak" version. (It is the path to the folder where the theriak executable of theriak-domino is located. In the version from 06.06.2023 this is usually the "programs" folder from the software.)

# 1.2.1 Install the package
ThorPT is public on the test PyPi servers and will be published on the offical servers soon. You can get the package by using the following command to import it to your python environment.

pip install -i https://test.pypi.org/simple/ thorpt-thorPT

# 1.2.2 Import the package
from thorpt_thorPT import start_ThorPT


# 1.2.3 Run the main module
The software will ask you for two input files via a GUI:
1. init file that defines parameters for the modelling
2. depending if you want to generate or use a preexisting P-T path

start_ThorPT.run_routine()


# 1.2.4 Plotting module with ThorPT
Plotting is based on the taufrir module. This module comprises pre defined plots for phase assemblages, oxygen isotopes and several plots for visualising the fluid production and transfer. Options are possible to save the image files and generate a gif.


First steps:
Import taufrir module
Read variables from a hdf5 with the function
```python
from thorpt_thorPT.taufrir02 import *
data = ThorPT_hdf5_reader()
data.open_ThorPT_hdf5()
```


Activate the plotting module - Predefined plots for evaluation:
```python
compPlot = ThorPT_plots(data.filename, data.mainfolder, data.rock, data.compiledrock)
```
Default plotting functions are then:
```python
compPlot.boxplot_to_GIF(rock_tag='NAME')

compPlot.pt_path_plot(rock_tag='NAME')

compPlot.permeability_plot(rock_tag='NAME')

compPlot.time_int_flux_plot(rock_tag='NAME')

compPlot.porosity_plot(rock_tag='NAME')

compPlot.release_fluid_volume_plot(rock_tag='NAME')
```