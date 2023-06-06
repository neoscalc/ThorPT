# ThorPT - This is ThorPT

pip install -i https://test.pypi.org/simple/ thorpt-thorPT

# Import the package
from thorpt_thorPT import start_ThorPT


# Run the main module
"""The software will ask you for two input files via a GUI:"""
"""1. init file that defines parameters for the modelling"""
"""2. depending if you want to generate or use a preexisting P-T path"""

start_ThorPT.run_routine()


# Plotting module "taufrir"

"""Plotting default results"""
from thorpt_thorPT.taufrir02 import *
"""Read variables from a hdf5 output file from ThorPT"""
data = ThorPT_hdf5_reader()
data.open_ThorPT_hdf5()
"""Activate the plotting module - Predefined plots for evaluation"""
compPlot = ThorPT_plots(data.filename, data.mainfolder, data.rock, data.compiledrock)
"""Default plotting functions are"""
compPlot.boxplot_to_GIF(rock_tag='NAME')
compPlot.pt_path_plot(rock_tag='NAME')
compPlot.permeability_plot(rock_tag='NAME')
compPlot.time_int_flux_plot(rock_tag='NAME')
compPlot.porosity_plot(rock_tag='NAME')
compPlot.release_fluid_volume_plot(rock_tag='NAME')