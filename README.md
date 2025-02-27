
<img src="https://tenor.com/view/thor-love-and-thunder-mjolnir-thor-thor-hammer-the-mighty-thor-gif-25431540" width="500" />
# ThorPT - This is ThorPT
ThorPT is a modelling package for fluid coupled geolgical processes in the crust. Fluid production and migration in the crust is fundamental for the earth's geochemical cycling and ThorPT is specified for modelling this fluid production and fluid transfer. It includes the latest developments in petrogeochemical modelling and helps the user to test different scenarios:
- Single rock modelling:
    1. Petrological purposes to test changes in mineral assemblages and coupled fluid production, e.g., while prograde metamorphism
    2. Mechanical testing of the rock based on the Mohr-Coulomb theory
    3. Quantification of the production of fluid masses

- Multi rock modelling (Work in progress):
    1. Stacks of rock sequences involving fluid transfer
    2. Outcrop modelling
    3. *Contact zones and fluid-rock ratios

# 1. Getting started with the package
Information and requirements:
There are currently two ways to use the package. Prerequisite in both cases is to get and install a working version of theriak-domino, because we need access to its executables. Check the software at
https://github.com/Theriak-Domino/theriak-domino

The latest updates incorporate packages written in the programming language Julia. Therefore, it is now necessary to have a version of Julia (>=v1.9; https://julialang.org/) installed on your machine (see juypter installation notebook for more details).

Further, you need the ThorPT input files to pass informations to the software. This is a simple text file and can be found in the latest release of ThorPT jupyter.

1.1 Recommended for all users

Quick start: Get the jupyter working directory that includes all instructions and installs ThorPT with its dependencies. You can find it in the latest release of ThorPT jupyter.

1.2 Recommended to users familair with Python

Install the package manually into your python environment and develop your own script.

# 1.1 Use the package with the jupyter environment
`Download the jupyter release that includes:`
- Jupyter notebook script to install ThorPT and all dependencies
- Jupyter script to run ThorPT and containing functions
- Datafiles directory with
    - initial file
    - examples for modelling


<div class="alert alert-block alert-danger">

__Before you start, make sure that you have followed all the steps and executed all the cells in the "install_ThorPT" Jupyter notebook!__

</div>


`Instructions: Creating the link to the software theriak`

- Declare the working directory of the theriak executable in order to do minimizations with ThorPT.
- Open the file "_init_.txt" and edit the first line after "Theriak:".<br>(You find the file in the released package you downloaded from github, it is located in the DataFiles folder)

- Insert the file path of theriak. In the latest versions of TheriakDomino it is located in the "programs" folder.<br>Windows example:"Theriak:C:\TheriakDominoWIN.2023.06.11\Programs"

# 1.2 Manual installation to your python environment
# 1.2.1 Install the package
ThorPT is public on the test PyPi servers and will be published on the offical servers soon. You can get the package by using the following command to import it to your python environment.

``` python
pip install thorpt
```

# 1.2.2 Import the package
from thorpt import start_ThorPT


# 1.2.3 Run the main module
The software will prompt you for two input files:
1. init file that defines parameters for the modelling
2. usage of a P-T path (pressure unit has to be defined in the init file, using "kbar" or "GPa")

start_ThorPT.run_main_routine()


# 1.2.4 Plotting module with ThorPT
Plotting is based on the ξόρκι module. This module comprises pre defined plots for phase assemblages, oxygen isotopes and several plots for visualising the fluid production and transfer. Options are possible to save the image files and generate a gif.


First steps:
from thorpt.xorki import *
Read variables from a hdf5 with the function
```python
from thorpt.xorki import *
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


# Funding and author

The development of this package was supported by the European Research Council (ERC) under the European Union's Horizon 2020 research and innovation programme (grant agreement No 850530).

Author: Thorsten A. Markmann (thorsten.markmann@unibe.ch)
