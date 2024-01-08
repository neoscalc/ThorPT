"""
Written by
Thorsten Markmann
thorsten.markmann@geo.unibe.ch
status: 11.06.2023
"""

# Plotting module for ThorPT
# Reading hdf5 file
# Suit of plotting functions for petrological modelling


import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import h5py
import os

from pathlib import Path
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from tkinter import *
from tkinter import filedialog
import imageio.v2
from imageio.v2 import imread
import copy
from dataclasses import dataclass, field
from dataclasses import fields


def file_opener():
    root = Tk()
    root.withdraw()
    root.update()
    filein = filedialog.askopenfilename(
        title="Select h5-file to read",
        filetypes=(
            ("hdf5 file", "*.hdf5"),
            ("All files", "*.*"))
    )
    root.update()
    return filein


def remove_items(test_list, item):
    """[removes item in a list - used in code to remove blanks or tabs in a list]

    Args:
        test_list ([list]): [a list with items read from a txt file]
        item ([string or number ...]): [can be anything - used in this case for symbols, blanks or tabs]

    Returns:
        [list]: [list of items now cleaned from the item defined]
    """
    # using list comprehension to perform the task
    res = [i for i in test_list if i != item]
    return res


def phases_and_colors_XMT(database, phases):

    # XMapTools mineral names and colors
    script_folder = Path(__file__).parent.absolute()
    file_to_open = script_folder / "DataFiles" / "XMap_MinColors.txt"
    minerals_XMT = pd.read_csv(file_to_open, header=16, delimiter='\t', names=[
                               'Name', 'R', 'G', 'B'])

    # Translation file - database to XMT names
    dot_indi = database.index('.')
    file_to_open = script_folder / "DataFiles" / \
        f"MINERAL_NAMES_{database[:dot_indi]}_to_XMT.txt"
    min_translation = pd.read_csv(file_to_open, delimiter='\t')

    # Iterating through phases and select colors
    colpal = sns.color_palette("flare", 20)
    color_set = []
    phase_set = []
    z = 0
    for mini in phases:
        if mini == 'fluid' or mini == 'H2O.liq':
            color_set.append((84/255, 247/255, 242/255))
            phase_set.append("Water")
        elif mini in list(min_translation['Database-name']):
            xmt_name = min_translation[min_translation['Database-name']
                                       == mini]['XMT-name'].iloc[0]
            minerals_XMT[minerals_XMT['Name'] == xmt_name]
            r = minerals_XMT[minerals_XMT['Name'] == xmt_name]['R'].iloc[0]
            g = minerals_XMT[minerals_XMT['Name'] == xmt_name]['G'].iloc[0]
            b = minerals_XMT[minerals_XMT['Name'] == xmt_name]['B'].iloc[0]
            color_set.append((r/255, g/255, b/255))
            phase_set.append(xmt_name)
        elif '_' in mini:
            mini_inid = mini.index('_')
            base_name = mini[:mini_inid]
            if base_name in list(min_translation['Database-name']):
                xmt_name = min_translation[min_translation['Database-name']
                                           == base_name]['XMT-name'].iloc[0]
                minerals_XMT[minerals_XMT['Name'] == xmt_name]
                r = minerals_XMT[minerals_XMT['Name'] == xmt_name]['R'].iloc[0]
                g = minerals_XMT[minerals_XMT['Name'] == xmt_name]['G'].iloc[0]
                b = minerals_XMT[minerals_XMT['Name'] == xmt_name]['B'].iloc[0]
                color_set.append((r/255, g/255, b/255))
                phase_set.append(xmt_name)
            else:
                color_set.append(colpal[z])
                phase_set.append(mini.title())
                z += 1
        else:
            color_set.append(colpal[z])
            phase_set.append(mini.title())
            z += 1

    return phase_set, color_set


def mineral_translation(database):
    """
    read DB_database in DataFiles file and stores to list
    - is the translation between theriak nad the DBOxygen dataset

    Returns:
        Dictionary: Including 'DB_phases', 'Theriak_phases', 'SSolution_Check', 'SSolution_set', 'Oxygen_count'
    """

    # Checks for database used and reads the name, to be used for translating file
    name = "SSModel_" + database + ".txt"

    # selecting path and file depending on OS
    mainfolder = Path(__file__).parent.absolute()
    data_folder = mainfolder / "DataFiles/"

    file_to_open = data_folder / name

    # data_folder = Path("DataFiles/")
    # file_to_open = data_folder / name

    # saving lines read from file to list
    read_translation = []
    with open(file_to_open, encoding='utf8') as f:
        reading = f.readlines()
        for line in reading:
            read_translation.append(line)

    # Variables to store
    phases_DB = []
    phases_Ther = []
    phases_if_ssolution = []
    phases_of_solution = []
    phase_counts = []
    # Extract mineral names from mineral data
    mineral_data_line = read_translation.index('*** MINERAL DATA ***\n')
    end_line = len(read_translation) - 16
    comp_list = read_translation[mineral_data_line+1:end_line]
    for num, item in enumerate(comp_list):
        line = item.split('\t')
        # function to remove all occurences of '' in list
        line = remove_items(line, '')
        # Store in comprehensive list
        phases_DB.append(line[0])
        line[1] = line[1].replace(' ', '')
        phases_Ther.append(line[1])
        phase_counts.append(line[2].rstrip())
        phases_if_ssolution.append(False)
        phases_of_solution.append(False)

    phase_counts = [float(ele) for ele in phase_counts]
    # extrqacting solid solution phases
    solid_solution_line = read_translation.index(
        '*** SOLID SOLUTIONS ***\n')
    end_line = read_translation.index('*** MINERAL DATA ***\n')
    comp_list = read_translation[solid_solution_line+1:end_line]
    comp_list = remove_items(comp_list, '\n')
    for i, line in enumerate(comp_list):
        line = line.rstrip()
        line = line.split()
        for el in line:
            if el == '>>':
                phases_DB.append(line[1])
                phases_Ther.append(line[2])
                phases_if_ssolution.append(True)
                cache1 = comp_list[i+1].rstrip()
                cache1 = cache1.split()
                cache2 = comp_list[i+2].rstrip()
                cache2 = cache2.split()
                phases_of_solution.append([cache1, cache2])
                ox_num_set = comp_list[i+3].rstrip()
                ox_num_set = ox_num_set.split()
                ox_num_set = [float(ele) for ele in ox_num_set]
                phase_counts.append(ox_num_set)

    translation_dic = {'DB_phases': phases_DB, 'Theriak_phases': phases_Ther,
                       'SSolution_Check': phases_if_ssolution, 'SSolution_set': phases_of_solution,
                       'Oxygen_count': phase_counts}

    return translation_dic


def Merge_phase_group(data):
    frame = data
    phases = list(frame.columns)
    li_BIO = []
    li_CCDO = []
    li_FSP = []
    li_CHLR = []
    li_PHNG = []
    li_ClAMP = []
    li_OMPH = []
    li_GARNET = []
    li_OLIVINE = []
    li_OPX = []
    li_SERP = []
    li_all = []
    rev_data = pd.DataFrame()

    for item in phases:
        if 'CCDO' in item:
            li_CCDO.append(item)
            li_all.append(item)
        if 'FSP' in item:
            li_FSP.append(item)
            li_all.append(item)
        if 'CHLR' in item:
            li_CHLR.append(item)
            li_all.append(item)
        if 'PHNG' in item:
            li_PHNG.append(item)
            li_all.append(item)
        if 'ClAMP' in item:
            li_ClAMP.append(item)
            li_all.append(item)
        if 'OMPH' in item:
            li_OMPH.append(item)
            li_all.append(item)
        if 'GARNET' in item:
            li_GARNET.append(item)
            li_all.append(item)
        if 'BIO' in item:
            li_BIO.append(item)
            li_all.append(item)
        if 'OLIVINE' in item:
            li_OLIVINE.append(item)
            li_all.append(item)
        if 'OPX' in item:
            li_OPX.append(item)
            li_all.append(item)
        if 'SERP' in item:
            li_SERP.append(item)
            li_all.append(item)

    selection = [li_BIO, li_CCDO, li_FSP, li_CHLR,
                 li_PHNG, li_ClAMP, li_OMPH, li_GARNET, li_OLIVINE, li_OPX, li_SERP]

    for collection in selection:
        chache = pd.DataFrame()
        phase_frame = pd.DataFrame()
        # Check each dataframe in collection and combines multiple
        for phase in collection:
            if chache.empty is True:
                chache = frame[phase]
            else:
                chache = chache + frame[phase]

        if not collection:
            pass
        else:
            name = collection[0]
            name = name[:name.index('_')]
            phase_frame[name] = chache

        rev_data = pd.concat([rev_data, phase_frame], axis=1)

    for phase in li_all:
        frame = frame.drop(phase, axis=1)

    frame = pd.concat([frame, rev_data], axis=1)

    frame = frame.T

    if 'LIQtc_h2oL' in frame.index and 'fluid' in frame.index:
        water_1 = frame.loc['LIQtc_h2oL']
        water_2 = frame.loc['fluid']
        water_2[water_2.isna()] = water_1[water_2.isna()]
        frame.loc['fluid'] = water_2
        frame = frame.drop('LIQtc_h2oL', axis=0)

    return frame


def Merge_phase_group_oxy(data):
    frame = data
    phases = list(frame.columns)
    li_BIO = []
    li_CCDO = []
    li_FSP = []
    li_CHLR = []
    li_PHNG = []
    li_ClAMP = []
    li_OMPH = []
    li_GARNET = []
    li_OLIVINE = []
    li_OPX = []
    li_SERP = []
    li_all = []
    rev_data = pd.DataFrame()

    for item in phases:
        if 'CCDO' in item:
            li_CCDO.append(item)
            li_all.append(item)
        if 'FSP' in item:
            li_FSP.append(item)
            li_all.append(item)
        if 'CHLR' in item:
            li_CHLR.append(item)
            li_all.append(item)
        if 'PHNG' in item:
            li_PHNG.append(item)
            li_all.append(item)
        if 'ClAMP' in item:
            li_ClAMP.append(item)
            li_all.append(item)
        if 'OMPH' in item:
            li_OMPH.append(item)
            li_all.append(item)
        if 'GARNET' in item:
            li_GARNET.append(item)
            li_all.append(item)
        if 'BIO' in item:
            li_BIO.append(item)
            li_all.append(item)
        if 'OLIVINE' in item:
            li_OLIVINE.append(item)
            li_all.append(item)
        if 'OPX' in item:
            li_OPX.append(item)
            li_all.append(item)
        if 'SERP' in item:
            li_SERP.append(item)
            li_all.append(item)

    selection = [li_BIO, li_CCDO, li_FSP, li_CHLR,
                 li_PHNG, li_ClAMP, li_OMPH, li_GARNET, li_OLIVINE, li_OPX, li_SERP]

    for collection in selection:
        chache = pd.DataFrame()
        phase_frame = pd.DataFrame()
        # Check each dataframe in collection and combines multiple
        for phase in collection:
            if chache.empty is True:
                chache = frame[phase].fillna(0.00)
            else:
                chache = chache + frame[phase].fillna(0.00)

        if not collection:
            pass
        else:
            name = collection[0]
            name = name[:name.index('_')]
            phase_frame[name] = chache

        rev_data = pd.concat([rev_data, phase_frame], axis=1)
    rev_data = rev_data.replace(0.00, np.nan)
    for phase in li_all:
        frame = frame.drop(phase, axis=1)

    frame = pd.concat([frame, rev_data], axis=1)

    frame = frame.T

    if 'LIQtc_h2oL' in frame.index and 'fluid' in frame.index:
        water_1 = frame.loc['LIQtc_h2oL']
        water_2 = frame.loc['fluid']
        water_2[water_2.isna()] = water_1[water_2.isna()]
        frame.loc['fluid'] = water_2
        frame = frame.drop('LIQtc_h2oL', axis=0)

    return frame


def create_gif(phase_data, mainfolder, filename, group_key, subfolder='default'):
    # reading images and put it to GIF
    frames = []
    for i in range(len(phase_data.index)):
        image = imread(
            f'{mainfolder}/img_{filename}/{subfolder}/{group_key}/img_{i}.png')
        frames.append(image)
    imageio.mimsave(
        f'{mainfolder}/img_{filename}/{subfolder}/{group_key}/output.gif', frames, duration=400)


# Progressbar init
def progress(percent=0, width=40):
    left = width * percent // 100
    right = width - left
    tags = "#" * left
    spaces = " " * right
    percents = f"{percent:.0f}%"
    print("\r[", tags, spaces, "]", percents, sep="", end="", flush=True)

def clean_frame(dataframe_to_clean, legend_phases, color_set):
            dataframe = dataframe_to_clean.copy()
            # copy information before manipulation
            colorframe = pd.DataFrame(color_set, index=legend_phases)
            new_color_set2 = []
            new_legend_phases = []
            new_dataframe = pd.DataFrame()

            # routine to filter multiple assigned phase columns
            for phase in legend_phases:
                # Merge dataframe rows of multiple entries with name of phase
                pcount = legend_phases.count(phase)
                if pcount > 1:
                    if phase in new_legend_phases:
                        pass
                    else:
                        # Sum the dataframe rows of multiple entries with name of phase
                        dat = dataframe.loc[phase].sum()
                        # Merge dat to new_dataframe
                        new_dataframe = pd.concat([new_dataframe, dat], axis=1)
                        new_legend_phases.append(phase)
                        new_color_set2.append(tuple(colorframe.loc[phase].mean()))
                else:
                    dat = dataframe.loc[phase]
                    # Merge dat to new_dataframe
                    new_dataframe = pd.concat([new_dataframe, dat], axis=1)
                    new_legend_phases.append(phase)
                    new_color_set2.append(tuple(colorframe.loc[phase]))

            new_dataframe = new_dataframe.T
            new_dataframe.index = new_legend_phases
            print('Cleaning done!')


            return new_dataframe, new_legend_phases, new_color_set2


def depth_porosity(data_box, num_rocks=3, rock_colors=["#0592c1", "#f74d53", '#10e6ad'], rock_symbol = ['d', 's']):

    for j, subdata in enumerate(data_box):
        step = int(len(subdata.compiledrock.all_porosity)/num_rocks)

        for i in range(num_rocks):
            frac = np.concatenate(subdata.compiledrock.all_frac_bool[step*i:step*(i+1)])
            x = np.concatenate(subdata.compiledrock.all_porosity[step*i:step*(i+1)])*100
            y = np.concatenate(subdata.compiledrock.all_depth[step*i:step*(i+1)])/1000

            y = y[frac>0]
            x = x[frac>0]
            plt.plot(x,y, rock_symbol[j], color=rock_colors[i], markeredgecolor='black')

    plt.ylabel("Depth [km]")
    plt.xlabel("Fluid-filled porosity @ extraction [Vol.%]")
    plt.legend(['Basalt', 'Sediment', 'Serpentinite'])
    plt.xscale('log')
    plt.ylim(100,0)
    plt.xlim(10e-4,10e2)
    plt.show()

def density_porosity(data_box, num_rocks=3, rock_colors=["#0592c1", "#f74d53", '#10e6ad'], rock_symbol = ['d', 's']):


    for j, subdata in enumerate(data_box):
        step = int(len(subdata.compiledrock.all_porosity)/num_rocks)

        for i in range(num_rocks):
            frac = np.concatenate(subdata.compiledrock.all_frac_bool[step*i:step*(i+1)])
            x = np.concatenate(subdata.compiledrock.all_porosity[step*i:step*(i+1)])*100
            y1 = np.concatenate(subdata.compiledrock.all_system_density_post[step*i:step*(i+1)])
            y2 = np.concatenate(subdata.compiledrock.all_system_density_pre[step*i:step*(i+1)])
            y = y1-y2
            y = y[frac>0]
            x = x[frac>0]
            plt.plot(x,y, rock_symbol[j], color=rock_colors[i], markeredgecolor='black')

    plt.ylabel(r"$\Delta$ density")
    plt.xlabel("Fluid-filled porosity at extraction")
    plt.legend(['Basalt', 'Sediment', 'Serpentinite'])
    plt.yscale('log')
    plt.xscale('log')
    plt.ylim(0,0.5)
    plt.xlim(0,25)
    plt.show()

def extraction_stress(data_box, num_rocks=3, rock_colors=["#0592c1", "#f74d53", '#10e6ad'], rock_symbol = ['d', 's'], rock_names=False):


    for j, subdata in enumerate(data_box):
        step = int(len(subdata.compiledrock.all_porosity)/num_rocks)

        for i in range(num_rocks):
            track_diff = []
            track_shear = []
            track_num_ext = []
            for k, item in enumerate(subdata.compiledrock.all_diffs[step*i:step*(i+1)]):
                diff = np.unique(item)[-1]
                fracture_bool = subdata.compiledrock.all_frac_bool[k+step*i]
                num_ext = len(fracture_bool[fracture_bool>0])
                track_diff.append(diff)
                # track_shear.append()
                track_num_ext.append(num_ext)
            if rock_names is False:
                plt.plot(track_diff, track_num_ext, rock_symbol[j], color=rock_colors[i], markersize=14, markeredgecolor='black')
            else:
                plt.plot(track_diff, track_num_ext, rock_symbol[i], color=rock_colors[j], markersize=14, markeredgecolor='black')
    plt.ylabel("# of extractions", fontsize=18)
    plt.xlabel("Differential stress [$\it{MPa}$]", fontsize=18)

    plt.yticks(fontsize=16)
    plt.xticks(fontsize=16)

    plt.show()


def resort_frame(y, legend_phases, color_set):
            # read the position of "Water" in legen_phases and move "Water" to the end of the legend_phases list
            # move the color of "Water" to the end of the color_set list
            color_set.append(color_set.pop(legend_phases.index('Water')))
            legend_phases.append(legend_phases.pop(legend_phases.index('Water')))
            # reindex the dataframe and move "Water" to the end of the dataframe
            y = y.reindex(legend_phases)
            return y, legend_phases, color_set

class Simple_binary_plot():

    def __init__(self, x, y):
        self.x = x
        self.y = y

    def path_plot(self, save=False):

        # TODO - Add exception when user aborts input or gives no input - Needs a default mode (x,y)
        x_ax_label = input("Please provide the x-axis label...")
        y_ax_label = input("Please provide the y-axis label...")
        plt.scatter(self.x, self.y, s=18, c='#7fffd4', markeredgecolor='black')
        plt.xlabel(x_ax_label)
        plt.ylabel(y_ax_label)

        # saving option from user input
        if save is True:
            plt.savefig(Path("plotting_results/PTt_path.png"), dpi=400)

@dataclass
class Phasecomp:
    name: str
    temperature: float = field(metadata={'unit': 'degree C'})
    pressure: float = field(metadata={'unit': 'bar'})
    moles: float
    volume: float = field(metadata={'unit': 'ccm'})
    volp: float = field(metadata={'unit': 'V%'})
    mass: float = field(metadata={'unit': 'g'})
    massp: float = field(metadata={'unit': 'wt%'})
    density: float = field(metadata={'unit': 'g/ccm'})
    elements: any
    volPmole: float = field(metadata={'unit': 'ccm/mol'})


# single rock data as dataclass
@dataclass
class Rock:
    name: str
    # temperature: float = field(metadata={'unit': 'degree C'})
    temperature: any
    pressure: any
    depth: any
    time: any

    systemVolPre: any
    systemVolPost: any
    fluidbefore: any

    tensile_strength: any
    differentialstress: any
    frac_bool: any
    permeability: any
    timeintflux: any

    arr_line: any
    geometry: any
    database: any
    phases: any
    phases2: any
    phase_data: any
    td_data: any
    oxygen_data: any
    bulk_deltao_pre: any
    bulk_deltao_post: any

    group_key: str
    extracted_fluid_volume: any
    extraction_boolean: any
    porosity: any
    time_int_flux: any
    time_int_flux2: any

    v_permeability: any
    v_timeintflux: any

    garnet: any
    garnets_bools: any


# Compiled data as dataclass
@dataclass
class CompiledData:
    all_ps: any
    all_tensile: any
    all_diffs: any
    all_frac_bool: any
    arr_line: any
    all_permea: any
    all_t_flux: any
    all_depth: any
    all_virtual_perm: any
    all_virtual_flux: any
    all_geometry: any
    all_system_vol_pre: any
    all_system_vol_post: any
    all_system_density_pre: any
    all_system_density_post: any
    all_fluid_before: any
    all_phases: any
    all_phases2: any
    all_database: any
    all_porosity: any


# HDF5 reader
class ThorPT_hdf5_reader():

    def __init__(self):
        self.rock = {}
        self.rock_names = []
        self.compiledrock = 0
        self.mainfolder = False
        self.filename = False

    def open_ThorPT_hdf5(self):

        # lists for compiling data
        all_ps = []
        all_tensile = []
        all_diffs = []
        all_frac_bool = []
        arr_line = []
        all_permea = []
        all_t_flux = []
        all_depth = []
        all_virtual_perm = []
        all_virtual_flux = []
        all_geometry = []
        all_system_vol_pre = []
        all_system_vol_post = []
        all_system_density_pre = []
        all_system_density_post = []
        all_fluid_before = []
        all_phases = []
        all_phases2 = []
        all_database = []
        all_td_data = []
        all_porosity = []
        garnets = {}
        sysv = {}
        sysv_post = {}
        garnets_bools = {}
        volume_data = {}

        o_file = file_opener()

        if len(o_file) < 1:
            print("No data file selected. Quitting routine.")
            quit()

        self.mainfolder = Path(o_file).parent.absolute()
        self.filename = o_file.split('/')[-1].split('_')[0].split('.')[0]

        # open the thorPT output file
        with h5py.File(o_file, 'r') as f:
            rocks = list(f.keys())
            for rock in rocks:
                self.rock[rock] = 0

            for group_key in rocks:

                # physical input data from modelling
                ts = np.array(f[group_key]['Parameters']['temperatures'])
                ps = np.array(f[group_key]['Parameters']['pressures'])
                vtime = np.array(f[group_key]['Parameters']['time_frame'])
                vtime = vtime.T
                extr_time = np.array(f[group_key]['FluidData']['extr_time'])
                frac_bool = np.array(f[group_key]['MechanicsData']['fracture bool'])
                depth = np.array(f[group_key]['Parameters']['depth'])
                geometry = list(f[group_key]['Parameters']['geometry'].asstr())

                tensile = np.float64(f[group_key]['Parameters']['tensile strength'])

                if 'line' in f[group_key].keys():
                    line = np.array(f[group_key]['Parameters']['line'])
                else:
                    line = False

                # Modelled params
                permea = np.array(f[group_key]['permeability'])
                t_flux = np.array(f[group_key]['time-int flux'])
                virtual_perm = np.array(f[group_key]['virtual permeability'])
                virtual_flux = np.array(f[group_key]['virtual time-int flux'])

                # Fluid extraction
                extr_d = pd.DataFrame(f[group_key]['FluidData']['extracted_fluid_data'])
                fluid_before = np.array(f[group_key]['MechanicsData']['st_fluid_before'])

                # get some physics arrays
                sys_physc = {}
                for item in f[group_key]['sys_physicals']:
                    sys_physc[item] = np.array(
                        f[group_key]['sys_physicals'][item])

                system_vol_pre = sys_physc['system_vol_pre']
                system_vol_post = sys_physc['system_vol_post']
                system_density_pre = sys_physc['system_density_pre']
                system_density_post = sys_physc['system_density_post']

                # ////////////////////////////////////////////////////////

                # Thermodynamic data

                # phases
                phases = list(f[group_key].attrs['Phases'])
                phases2 = (f[group_key].attrs['Phases'])
                database = str(f[group_key].attrs['database'])

                # retrieve all the thermodyamica physical data
                df_var_d = {}
                for item in f[group_key]['df_var_dictionary'].keys():
                    df_var_d[item] = pd.DataFrame(
                        f[group_key]['df_var_dictionary'][item])
                    if len(phases) > len(df_var_d[item].columns):
                        pass
                    else:
                        df_var_d[item].columns = phases
                # element data
                td_data = pd.DataFrame(f[group_key]['df_element_total'])
                td_data_tag = list(f[group_key].attrs['el_index'])
                td_data = td_data.T
                td_data.columns = td_data_tag
                # Get the chemical potential data
                pot_tag = list(f[group_key].attrs['pot_tag'])
                p_data = pd.DataFrame(f[group_key]['pot_data'])
                p_data.index = pot_tag

                # Oxygen fractionation data
                # oxygen isotope data
                d_oxy = pd.DataFrame(f[group_key]['IsotopeData']['save_oxygen'])
                d_oxy_phases = list(f[group_key].attrs['oxy_data_phases'])
                d_oxy.columns = d_oxy_phases
                save_bulk_oxygen_pre = np.array(
                    f[group_key]['IsotopeData']['save_bulk_oxygen_pre'])
                save_bulk_oxygen_post = np.array(
                    f[group_key]['IsotopeData']['save_bulk_oxygen_post'])

                # geometry factor
                bloc_a = np.float64(geometry[0])
                bloc_b = np.float64(geometry[1])
                bloc_c = np.float64(geometry[2])
                area = bloc_b*bloc_c
                xxx = bloc_a
                size = bloc_a * bloc_b * bloc_c

                v_pre = sys_physc['system_vol_pre']
                v_post = sys_physc['system_vol_post']
                v_fluid_extr = v_pre - v_post
                scaling_factor = size*1_000_000/system_vol_pre
                v_fluid_extr = v_fluid_extr*scaling_factor

                unfiltered_porosity = fluid_before/v_pre
                unfiltered_int_flux = v_fluid_extr/area

                # Regional time int flux calculation from Ague 2013 (Fluid flow in deep crust)
                mass_data = copy.deepcopy(df_var_d['df_wt%'])
                mass_data.columns = phases
                mass_data[np.isnan(mass_data) == True] = 0

                volume_data = copy.deepcopy(df_var_d['df_volume'])
                volume_data.columns = phases
                volume_data[np.isnan(volume_data) == True] = 0

                mass_abs_data = copy.deepcopy(df_var_d['df_wt'])
                mass_abs_data.columns = phases
                mass_abs_data[np.isnan(mass_abs_data) == True] = 0

                if 'fluid' in volume_data.columns:
                    solid_volumes = volume_data.T.sum()-volume_data['fluid']
                    solid_weight = mass_abs_data.T.sum(
                    )-mass_abs_data['fluid']
                else:
                    solid_volumes = volume_data.T.sum()
                    solid_weight = mass_abs_data.T.sum()

                solid_density = np.array(
                    solid_weight)/np.array(solid_volumes)*1000

                # massperc_fluid = mass_data['fluid']
                if 'fluid' in mass_abs_data.columns:
                    massperc_fluid = mass_abs_data['fluid'] / \
                        (mass_abs_data.T.sum()-mass_abs_data['fluid'])
                else:
                    massperc_fluid = 0

                q_ti = massperc_fluid*solid_density * \
                    (1-unfiltered_porosity)*bloc_a/area

                # the new permeability
                # ANCHOR
                mü_water = 0.0001
                if 'fluid' in mass_abs_data.columns:
                    density_cont = solid_density - \
                        mass_abs_data['fluid']/volume_data['fluid']
                else:
                    density_cont = 0
                permeability2 = q_ti/(151000*365*24*60*60) * \
                    mü_water / 9.81 / density_cont

                # extraction marker
                # NOTE - new masking method
                extraction_boolean = np.array(frac_bool, dtype=bool)
                extraction_boolean = np.invert(extraction_boolean)



                #######################################################
                el = list(f[group_key].attrs['garnet'])
                params = {}
                for item in f[group_key]['garnet']:
                    if item == 'elements':
                        params[item] = pd.DataFrame(f[group_key]['garnet'][item], index=el)
                    elif item == 'name':
                        params[item] = list(f[group_key]['garnet'][item])
                    else:
                        params[item] = np.array(f[group_key]['garnet'][item])

                phase = Phasecomp(params['name'],
                                params['temperature'],
                                params['pressure'],
                                params['moles'],
                                params['volume'],
                                params['volp'],
                                params['mass'],
                                params['massp'],
                                params['density'],
                                params['elements'],
                                params['VolumePMole']
                                )
                garnets[group_key] = phase

                garnets_bools[group_key] = np.array(f[group_key]['garnet_check'])

                #######################################################
                # Write the dataclass
                self.rock[group_key] = Rock(
                    name=None,
                    temperature=ts,
                    pressure=ps,
                    depth=depth,
                    time=vtime[1],
                    systemVolPost=v_post,
                    systemVolPre=v_pre,
                    fluidbefore=fluid_before,
                    tensile_strength=tensile,
                    differentialstress=0,
                    frac_bool=frac_bool,
                    timeintflux=0,
                    arr_line=line,
                    geometry=geometry,
                    database=database,
                    phases=phases,
                    phases2=phases2,
                    phase_data=df_var_d,
                    td_data=td_data,
                    oxygen_data=d_oxy,
                    bulk_deltao_pre=save_bulk_oxygen_pre,
                    bulk_deltao_post=save_bulk_oxygen_post,
                    group_key=group_key,
                    extraction_boolean=extraction_boolean,
                    extracted_fluid_volume=v_fluid_extr,
                    permeability=permea,
                    porosity=unfiltered_porosity,
                    time_int_flux=unfiltered_int_flux,
                    time_int_flux2=q_ti,
                    v_permeability=0,
                    v_timeintflux=0,
                    garnet=garnets,
                    garnets_bools=garnets_bools)

                #######################################################
                # Compiling section
                all_ps.append(ps)
                all_tensile.append(tensile)
                all_diffs.append(np.array(f[group_key]['Parameters']['diff. stress']))
                all_frac_bool.append(frac_bool)
                arr_line.append(line)
                all_permea.append(permea)
                all_t_flux.append(np.array(q_ti))
                all_depth.append(depth)
                all_virtual_perm.append(virtual_perm)
                all_virtual_flux.append(virtual_flux)
                all_geometry.append(geometry)
                all_system_vol_pre.append(system_vol_pre)
                all_system_vol_post.append(system_vol_post)
                all_system_density_pre.append(system_density_pre)
                all_system_density_post.append(system_density_post)
                all_fluid_before.append(fluid_before)
                all_phases.append(phases)
                all_phases2.append(phases2)
                all_database.append(database)
                all_td_data.append(td_data)
                all_porosity.append(unfiltered_porosity)

                self.compiledrock = CompiledData(
                    all_ps,
                    all_tensile,
                    all_diffs,
                    all_frac_bool,
                    arr_line,
                    all_permea,
                    all_t_flux,
                    all_depth,
                    all_virtual_perm,
                    all_virtual_flux,
                    all_geometry,
                    all_system_vol_pre,
                    all_system_vol_post,
                    all_system_density_pre,
                    all_system_density_post,
                    all_fluid_before,
                    all_phases,
                    all_phases2,
                    all_database,
                    all_porosity)


# Module of plotting functions for reducing the data from modelling with ThorPT
class ThorPT_plots():

    def __init__(self, filename, mainfolder, rockdata, compiledrockdata):

        self.filename = filename
        self.mainfolder = mainfolder
        self.rockdic = rockdata
        self.comprock = compiledrockdata

    def boxplot_to_GIF(self, rock_tag, img_save=False, gif_save=False):

        if gif_save is True:
            img_save = True

        ts = self.rockdic[rock_tag].temperature
        database = self.rockdic[rock_tag].database
        phases = self.rockdic[rock_tag].phases
        phases2 = self.rockdic[rock_tag].phases2
        phase_data = self.rockdic[rock_tag].phase_data['df_N']
        group_key = self.rockdic[rock_tag].group_key

        # Clean up dataframe to be used by mineral names and no NaN
        phase_data.columns = phases
        phase_data[np.isnan(phase_data) == True] = 0

        # saving folder
        subfolder = 'phase_modes'

        # XMT naming and coloring
        legend_phases, color_set = phases_and_colors_XMT(database, phases2)

        # when image save is active
        if img_save is True:
            print("\n\nGenerating boxplots images. Please wait...")
            # Start progress bar
            k = 0
            kk = len(ts)
            progress(int(k/kk)*100)

            os.makedirs(
                f'{self.mainfolder}/img_{self.filename}', exist_ok=True)

            for i in range(len(phase_data.index)):
                mole_fractions = phase_data.iloc[i, :] / \
                    np.sum(phase_data.iloc[i, :])*100

                os.makedirs(
                    f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}', exist_ok=True)
                # -------------------------------------------------------------------
                fig = plt.figure(constrained_layout=False,
                                 facecolor='0.9', figsize=(9, 9))
                gs = fig.add_gridspec(nrows=3, ncols=3, left=0.15,
                                      right=0.95, hspace=0.6, wspace=0.5)
                # All the single plots
                # plot 1
                y_offset = 0
                ax2 = fig.add_subplot(gs[:, :2])
                for t, val in enumerate(np.array(mole_fractions)):
                    if val == np.float64(0):
                        ax2.bar(
                            1, val, bottom=0, color=color_set[t], edgecolor='black', linewidth=0.1, label=legend_phases[t])
                    else:
                        ax2.bar(1, val, bottom=y_offset,
                                color=color_set[t], edgecolor='black', linewidth=0.1, label=legend_phases[t])
                    y_offset = y_offset + val
                # legend
                handles, labels = ax2.get_legend_handles_labels()
                legend_properties = {'weight': 'bold'}
                ax2.legend(handles[::-1], labels[::-1],
                           title='Phases', bbox_to_anchor=(1.2, 0.8), fontsize=18)
                ax2.get_xaxis().set_visible(False)
                ax2.get_yaxis().set_visible(False)
                ax2.set_aspect(0.02)

                # Finishing plot - save
                plt.savefig(f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}/img_{i}.png',
                            transparent=False, facecolor='white')

                # plt.pause(2)
                plt.close()

                k += 1
                ic = int(np.ceil(k/kk*100))
                print("=====Progress=====")
                progress(ic)

        # reading images and put it to GIF
        if gif_save is True:
            # call gif function
            create_gif(phase_data, self.mainfolder, self.filename,
                       group_key, subfolder=subfolder)

    def pt_path_plot(self, rock_tag, img_save=False, gif_save=False):

        if gif_save is True:
            img_save = True

        # Variables to be used
        ts = self.rockdic[rock_tag].temperature
        ps = self.rockdic[rock_tag].pressure
        phases = self.rockdic[rock_tag].phases
        phase_data = self.rockdic[rock_tag].phase_data['df_N']
        group_key = self.rockdic[rock_tag].group_key
        # Clean up dataframe to be used by mineral names and no NaN
        phase_data.columns = phases
        phase_data[np.isnan(phase_data) == True] = 0

        # saving folder
        subfolder = 'PT_path'

        if img_save is True:

            print("\n\nGenerating PT-path images. Please wait...")
            # Start progress bar
            k = 0
            kk = len(ts)
            progress(int(k/kk)*100)

            os.makedirs(
                f'{self.mainfolder}/img_{self.filename}', exist_ok=True)

            for i in range(len(phase_data.index)):
                os.makedirs(
                    f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}', exist_ok=True)
                # -------------------------------------------------------------------
                fig = plt.figure(constrained_layout=False,
                                 facecolor='0.9', figsize=(9, 9))
                gs = fig.add_gridspec(nrows=3, ncols=3, left=0.15,
                                      right=0.95, hspace=0.6, wspace=0.5)

                # plot 6
                ax5 = fig.add_subplot(gs[:, :])
                ax5.plot(ts[:i+1], ps[:i+1]/10_000, 'd--',
                         color='black', markersize=4)
                ax5.plot(ts[i:i+1], ps[i:i+1]/10_000, 'd--',
                         color='#7fffd4', markersize=8, markeredgecolor='black')
                # ax5 figure features
                ax5.set_xlim(200, 750)
                ax5.set_ylim(0, 3)
                ax5.set_ylabel("Pressure [GPa]", fontsize=18)
                ax5.set_xlabel("Temperature [°C]", fontsize=18)
                ax5.tick_params(axis='both', labelsize=18)

            # Finishing plot - save
                plt.savefig(f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}/img_{i}.png',
                            transparent=False, facecolor='white')
                # plt.pause(2)
                plt.close()

                # tail of progress bar
                k += 1
                ic = int(np.ceil(k/kk*100))
                print("=====Progress=====")
                progress(ic)
        else:
            # image framework
            fig = plt.figure(constrained_layout=False,
                             facecolor='0.9', figsize=(9, 9))
            gs = fig.add_gridspec(nrows=3, ncols=3, left=0.15,
                                  right=0.95, hspace=0.6, wspace=0.5)
            ax5 = fig.add_subplot(gs[:, :])
            for i in range(len(phase_data.index)):
                ax5.plot(ts[:i+1], ps[:i+1]/10_000, 'd--',
                         color='black', markersize=4)
                ax5.plot(ts[i:i+1], ps[i:i+1]/10_000, 'd--',
                         color='#7fffd4', markersize=8, markeredgecolor='black')
            # ax5 figure features
            ax5.set_xlim(200, 750)
            ax5.set_ylim(0, 3)
            ax5.set_ylabel("Pressure [GPa]", fontsize=18)
            ax5.set_xlabel("Temperature [°C]", fontsize=18)
            ax5.tick_params(axis='both', labelsize=18)
            plt.show()

        if gif_save is True:
            # call gif function
            create_gif(phase_data, self.mainfolder, self.filename,
                       group_key, subfolder=subfolder)

    def permeability_plot(self, rock_tag, img_save=False, gif_save=False):

        if gif_save is True:
            img_save = True

        group_key = self.rockdic[rock_tag].group_key

        phases = self.rockdic[rock_tag].phases
        phase_data = self.rockdic[rock_tag].phase_data['df_N']
        # Clean up dataframe to be used by mineral names and no NaN
        phase_data.columns = phases
        phase_data[np.isnan(phase_data) == True] = 0

        ts = self.rockdic[rock_tag].temperature
        ps = self.rockdic[rock_tag].pressure
        depth = self.rockdic[rock_tag].depth

        # Filter method of data array
        permea = self.rockdic[rock_tag].permeability
        extraction_boolean = self.rockdic[rock_tag].extraction_boolean

        if False in extraction_boolean:
            filtered_permeability = np.ma.masked_array(
                permea, extraction_boolean)
        else:
            filtered_permeability = np.ma.masked_array(
                np.zeros(len(extraction_boolean)), extraction_boolean)

        if len(np.unique(self.rockdic[rock_tag].frac_bool)) == 1 and np.unique(self.rockdic[rock_tag].frac_bool)[-1] == 0:
            filtered_permeability = np.zeros(
                len(self.rockdic[rock_tag].frac_bool))

        # saving folder
        subfolder = 'permeability'

        if img_save is True:

            print("\n\nGenerating permeability plot images. Please wait...")
            # Start progress bar
            k = 0
            kk = len(ts)
            progress(int(k/kk)*100)

            os.makedirs(
                f'{self.mainfolder}/img_{self.filename}', exist_ok=True)

            for i in range(len(phase_data.index)):
                os.makedirs(
                    f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}', exist_ok=True)
                # -------------------------------------------------------------------
                fig = plt.figure(constrained_layout=False,
                                 facecolor='0.9', figsize=(9, 9))
                gs = fig.add_gridspec(nrows=3, ncols=3, left=0.15,
                                      right=0.95, hspace=0.6, wspace=0.5)

                ax5 = fig.add_subplot(gs[:, :])
                if False in filtered_permeability.mask[:i+1]:
                    # plot 4
                    ax5.plot(
                        filtered_permeability[:i+1], depth[:i+1]/1000, 'd--', color='black', alpha=0.7)
                    ax5.plot(filtered_permeability[i:i+1], depth[i:i+1]/1000, 'd',
                             color='#7fffd4', markersize=8, markeredgecolor='black')
                else:
                    # plot 4
                    ax5.plot(np.ones(i+1)*1e-30, depth[:i+1]/1000, 'd')

                # add Manning & Ingebritsen to plot 3
                z = np.arange(1, 110, 1)
                c_k = 10**(-14 - 3.2 * np.log10(z))
                ax5.plot(c_k, z, '--', color='#b18e4e', linewidth=3)

                # ax5 figure features
                ax5.set_title('T:{:.2f} °C P:{:.2f} GPa'.format(
                    ts[i], ps[i]/10000), fontsize=18)
                ax5.set_xlim(1e-26, 1e-14)
                ax5.set_xscale('log')
                ax5.set_ylim(100, 0)
                ax5.set_xlabel("permeability [$m^{2}$]", fontsize=18)
                ax5.set_ylabel("Depth [km]", fontsize=18)
                ax5.fill_between([1e-30, 1e-21], [100, 100],
                                 color="#c5c5c5", edgecolor='black', hatch="/", alpha=0.4)
                ax5.tick_params(axis='both', labelsize=20)
                ax5.annotate("Manning &\nIngebritsen\n(1999)",
                             (4e-17, 15), fontsize=16)

            # Finishing plot - save
                plt.savefig(f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}/img_{i}.png',
                            transparent=False, facecolor='white')
                # plt.pause(2)
                plt.close()

                # tail of progress bar
                k += 1
                ic = int(np.ceil(k/kk*100))
                print("=====Progress=====")
                progress(ic)

        else:
            fig = plt.figure(constrained_layout=False,
                             facecolor='0.9', figsize=(9, 9))
            gs = fig.add_gridspec(nrows=3, ncols=3, left=0.15,
                                  right=0.95, hspace=0.6, wspace=0.5)
            ax5 = fig.add_subplot(gs[:, :])
            for i in range(len(phase_data.index)):
                if False in filtered_permeability.mask[:i+1]:
                    # plot 4
                    ax5.plot(
                        filtered_permeability[:i+1], depth[:i+1]/1000, 'd--', color='black', alpha=0.7)
                    ax5.plot(filtered_permeability[i:i+1], depth[i:i+1]/1000, 'd',
                             color='#7fffd4', markersize=8, markeredgecolor='black')
                else:
                    # plot 4
                    ax5.plot(np.ones(i+1)*1e-30, depth[:i+1]/1000, 'd')

            # add Manning & Ingebritsen to plot 3
            z = np.arange(1, 110, 1)
            c_k = 10**(-14 - 3.2 * np.log10(z))
            ax5.plot(c_k, z, '--', color='#b18e4e', linewidth=3)

            # ax5 figure features
            ax5.set_title('T:{:.2f} °C P:{:.2f} GPa'.format(
                ts[i], ps[i]/10000), fontsize=18)
            ax5.set_xlim(1e-26, 1e-14)
            ax5.set_xscale('log')
            ax5.set_ylim(100, 0)
            ax5.set_xlabel("permeability [$m^{2}$]", fontsize=18)
            ax5.set_ylabel("Depth [km]", fontsize=18)
            ax5.fill_between([1e-30, 1e-21], [100, 100],
                             color="#c5c5c5", edgecolor='black', hatch="/", alpha=0.4)
            ax5.tick_params(axis='both', labelsize=20)
            ax5.annotate("Manning &\nIngebritsen\n(1999)",
                         (4e-17, 15), fontsize=16)
            plt.show()

        # reading images and put it to GIF
        if gif_save is True:
            # call gif function
            create_gif(phase_data, self.mainfolder, self.filename,
                       group_key, subfolder=subfolder)

    def time_int_flux_plot(self, rock_tag, img_save=False, gif_save=False):

        group_key = self.rockdic[rock_tag].group_key

        phases = self.rockdic[rock_tag].phases
        phase_data = self.rockdic[rock_tag].phase_data['df_N']
        # Clean up dataframe to be used by mineral names and no NaN
        phase_data.columns = phases
        phase_data[np.isnan(phase_data) == True] = 0

        ts = self.rockdic[rock_tag].temperature
        ps = self.rockdic[rock_tag].pressure
        depth = self.rockdic[rock_tag].depth

        # Filter method of data array
        unfiltered_int_flux = self.rockdic[rock_tag].time_int_flux
        q_ti = self.rockdic[rock_tag].time_int_flux2
        extraction_boolean = self.rockdic[rock_tag].extraction_boolean

        if False in extraction_boolean:
            filtered_int_flux = np.ma.masked_array(
                unfiltered_int_flux, extraction_boolean)
            filtered_q_ti = np.ma.masked_array(
                q_ti, extraction_boolean)
            regional_filtered_flux = filtered_q_ti


        if len(np.unique(self.rockdic[rock_tag].frac_bool)) == 1 and np.unique(self.rockdic[rock_tag].frac_bool)[-1] == 0:
            filtered_q_ti = np.zeros(len(self.rockdic[rock_tag].frac_bool))
            filtered_int_flux = np.zeros(len(self.rockdic[rock_tag].frac_bool))
            regional_filtered_flux = filtered_q_ti


        if len(extraction_boolean) == 0:
            regional_filtered_flux = unfiltered_int_flux
            arr = np.invert(np.array(unfiltered_int_flux, dtype=bool))
            filtered_int_flux = np.ma.masked_array(
                regional_filtered_flux, arr)
        # saving subfolder
        subfolder = 'time_int_flux'

        if gif_save is True:

            print("\n\nGenerating time-int. fluid flux images. Please wait...")
            # Start progress bar
            k = 0
            kk = len(ts)
            progress(int(k/kk)*100)

            os.makedirs(
                f'{self.mainfolder}/img_{self.filename}', exist_ok=True)

            for i in range(len(phase_data.index)):
                os.makedirs(
                    f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}', exist_ok=True)
                # -------------------------------------------------------------------
                fig = plt.figure(constrained_layout=False,
                                 facecolor='0.9', figsize=(9, 9))
                gs = fig.add_gridspec(nrows=3, ncols=3, left=0.15,
                                      right=0.95, hspace=0.6, wspace=0.5)

                ax5 = fig.add_subplot(gs[:, :])
                # ax5.plot(unfiltered_int_flux[:i+1],
                #          depth[:i+1]/1000, 'd--', color='green')
                # ax5.plot(regional_flux[:i+1],
                #          depth[:i+1]/1000, 'd--', color='green', markersize=3, alpha=0.5)
                if False in filtered_int_flux.mask[:i+1]:
                    # plot 4
                    # ax5.plot(filtered_int_flux[:i+1], depth[:i+1]/1000, 'd--', color='black', alpha=0.7)
                    # ax5.plot(filtered_int_flux[i:i+1], depth[i:i+1]/1000, 'd',
                    #          color='#7fffd4', markersize=8, markeredgecolor='black')
                    ax5.plot(
                        regional_filtered_flux[:i+1], depth[:i+1]/1000, 'd--', color='black', alpha=0.7, markersize=10)
                    ax5.plot(regional_filtered_flux[i:i+1], depth[i:i+1]/1000, 'd',
                             color='#7fffd4', markersize=15, markeredgecolor='black')
                else:
                    # plot 4
                    ax5.plot(np.ones(i+1)*1e-30, depth[:i+1]/1000, 'd')

                # ax5 figure features
                ax5.set_title('T:{:.2f} °C P:{:.2f} GPa'.format(
                    ts[i], ps[i]/10000), fontsize=18)
                ax5.set_xscale('log')
                ax5.set_xlim(10**0, 10**6)
                ax5.set_ylim(100, 0)
                ax5.set_xlabel("Time int.-flux [$m^{3}/m^{2}$]", fontsize=18)
                ax5.set_ylabel("Depth [km]", fontsize=18)

                # pervasive range
                """ax5.fill_between([10**(2), 10**(4)], [100, 100],
                                color="#a9d5b2", edgecolor='black', alpha=0.1)"""
                ax5.fill_between([10**(1), 10**(4)], [100, 100],
                                 cmap="Purples", alpha=0.2, step='mid')
                # channelized range
                ax5.fill_between([10**(4), 10**(6)], [100, 100],
                                 color="#ca638f", alpha=0.1)
                # regional range
                """ax5.fill_between([10**(2.7-0.5), 10**(2.7+0.5)], [100, 100],
                                color="#c5c5c5", edgecolor='black', alpha=0.5)"""

                ax5.tick_params(axis='both', labelsize=20)
                ax5.vlines(10**4, 0, 100, linestyles='--',
                           color='black', linewidth=4)

                props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

                # ax5.annotate("Regional range", (10**(2.7+0.5-0.3), 25), fontsize=16, bbox=props, rotation=90)
                ax5.annotate("Channelized", (10**(4+0.2), 5),
                             fontsize=16, bbox=props)
                ax5.annotate("Pervasive", (10**(1+0.2), 5),
                             fontsize=16, bbox=props)
            # Finishing plot - save
                plt.savefig(f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}/img_{i}.png',
                            transparent=False, facecolor='white')
                # plt.pause(2)
                plt.close()

                # tail of progress bar
                k += 1
                ic = int(np.ceil(k/kk*100))
                print("=====Progress=====")
                progress(ic)

        else:
            fig = plt.figure(constrained_layout=False,
                             facecolor='0.9', figsize=(9, 9))
            gs = fig.add_gridspec(nrows=3, ncols=3, left=0.15,
                                  right=0.95, hspace=0.6, wspace=0.5)
            ax5 = fig.add_subplot(gs[:, :])
            for i in range(len(phase_data.index)):
                # ax5.plot(unfiltered_int_flux[:i+1],
                #          depth[:i+1]/1000, 'd--', color='green')
                # ax5.plot(regional_flux[:i+1],
                #          depth[:i+1]/1000, 'd--', color='green', markersize=3, alpha=0.5)
                # Test whether it is a masked array
                if np.ma.isMaskedArray(filtered_int_flux) is True:
                    if False in filtered_int_flux.mask[:i+1]:
                        # plot 4
                        # ax5.plot(filtered_int_flux[:i+1], depth[:i+1]/1000, 'd--', color='black', alpha=0.7)
                        # ax5.plot(filtered_int_flux[i:i+1], depth[i:i+1]/1000, 'd',
                        #          color='#7fffd4', markersize=8, markeredgecolor='black')
                        ax5.plot(
                            regional_filtered_flux[:i+1], depth[:i+1]/1000, 'd', color='black', alpha=0.7, markersize=10)
                        ax5.plot(regional_filtered_flux[i:i+1], depth[i:i+1]/1000, 'd',
                                color='#7fffd4', markersize=15, markeredgecolor='black')
                    else:
                        # plot 4
                        ax5.plot(np.ones(i+1)*1e-30, depth[:i+1]/1000, 'd')
                else:
                    # plot 4
                    ax5.plot(np.ones(i+1)*1e-30, depth[:i+1]/1000, 'd')

            # ax5 figure features
            ax5.set_xscale('log')
            ax5.set_xlim(10**0, 10**6)
            ax5.set_ylim(100, 0)
            ax5.set_xlabel("Time int.-flux [$m^{3}/m^{2}$]", fontsize=18)
            ax5.set_ylabel("Depth [km]", fontsize=18)
            # pervasive range
            """ax5.fill_between([10**(2), 10**(4)], [100, 100],
                            color="#a9d5b2", edgecolor='black', alpha=0.1)"""
            ax5.fill_between([10**(1), 10**(4)], [100, 100],
                             cmap="Purples", alpha=0.2, step='mid')
            # channelized range
            ax5.fill_between([10**(4), 10**(6)], [100, 100],
                             color="#ca638f", alpha=0.1)
            # regional range
            """ax5.fill_between([10**(2.7-0.5), 10**(2.7+0.5)], [100, 100],
                            color="#c5c5c5", edgecolor='black', alpha=0.5)"""
            ax5.tick_params(axis='both', labelsize=20)
            ax5.vlines(10**4, 0, 100, linestyles='--',
                       color='black', linewidth=4)
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            # ax5.annotate("Regional range", (10**(2.7+0.5-0.3), 25), fontsize=16, bbox=props, rotation=90)
            ax5.annotate("Channelized", (10**(4+0.2), 5),
                         fontsize=16, bbox=props)
            ax5.annotate("Pervasive", (10**(1+0.2), 5),
                         fontsize=16, bbox=props)

            if img_save is True:
                os.makedirs(
                    f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}', exist_ok=True)
                plt.savefig(f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}/fluid_flux_summary.png',
                            transparent=False, facecolor='white')
                # plt.pause(2)
                plt.close()
            else:
                plt.show()

        # reading images and put it to GIF
        if gif_save is True:
            # call gif function
            create_gif(phase_data, self.mainfolder, self.filename,
                       group_key, subfolder=subfolder)

    def time_int_flux_summary(self, img_save=False, rock_colors=["#0592c1", "#f74d53", '#10e6ad']):

        # saving subfolder
        subfolder = 'time_int_flux'

        fig = plt.figure(constrained_layout=False,
                        facecolor='0.9', figsize=(9, 9))
        gs = fig.add_gridspec(nrows=3, ncols=3, left=0.15,
                            right=0.95, hspace=0.6, wspace=0.5)
        ax5 = fig.add_subplot(gs[:, :])

        for i, rock_tag in enumerate(self.rockdic.keys()):

            if len(self.rockdic.keys()) == 3:
                rock_color = rock_colors[i]
            else:
                rock_color = "#7fffd4"

            group_key = self.rockdic[rock_tag].group_key
            phases = self.rockdic[rock_tag].phases
            phase_data = self.rockdic[rock_tag].phase_data['df_N']
            # Clean up dataframe to be used by mineral names and no NaN
            phase_data.columns = phases
            phase_data[np.isnan(phase_data) == True] = 0
            ts = self.rockdic[rock_tag].temperature
            ps = self.rockdic[rock_tag].pressure
            depth = self.rockdic[rock_tag].depth
            # Filter method of data array
            unfiltered_int_flux = self.rockdic[rock_tag].time_int_flux
            q_ti = self.rockdic[rock_tag].time_int_flux2
            extraction_boolean = self.rockdic[rock_tag].extraction_boolean
            if False in extraction_boolean:
                filtered_int_flux = np.ma.masked_array(
                    unfiltered_int_flux, extraction_boolean)
                filtered_q_ti = np.ma.masked_array(
                    q_ti, extraction_boolean)
                regional_filtered_flux = filtered_q_ti
            if len(np.unique(self.rockdic[rock_tag].frac_bool)) == 1 and np.unique(self.rockdic[rock_tag].frac_bool)[-1] == 0:
                filtered_q_ti = np.zeros(len(self.rockdic[rock_tag].frac_bool))
                filtered_int_flux = np.zeros(len(self.rockdic[rock_tag].frac_bool))
                regional_filtered_flux = filtered_q_ti
            if len(extraction_boolean) == 0:
                regional_filtered_flux = unfiltered_int_flux
                arr = np.invert(np.array(unfiltered_int_flux, dtype=bool))
                filtered_int_flux = np.ma.masked_array(
                    regional_filtered_flux, arr)

            if np.ma.isMaskedArray(filtered_int_flux) is True:
                if False in filtered_int_flux.mask[:]:
                    # plot 4
                    # ax5.plot(filtered_int_flux[:i+1], depth[:i+1]/1000, 'd--', color='black', alpha=0.7)
                    # ax5.plot(filtered_int_flux[i:i+1], depth[i:i+1]/1000, 'd',
                    #          color='#7fffd4', markersize=8, markeredgecolor='black')
                    ax5.plot(regional_filtered_flux, depth/1000, 'd',
                            color=rock_color, markersize=15, markeredgecolor='black')
                else:
                    # plot 4
                    ax5.plot(np.ones(len(filtered_int_flux))*1e-30, depth/1000, 'd')
            else:
                # plot 4
                ax5.plot(np.ones(len(filtered_int_flux))*1e-30, depth/1000, 'd')

        # General plot edits
        # ax5 figure features
        ax5.set_xscale('log')
        ax5.set_xlim(10**0, 10**6)
        ax5.set_ylim(100, 0)
        ax5.set_xlabel("Time int.-flux [$m^{3}/m^{2}$]", fontsize=18)
        ax5.set_ylabel("Depth [km]", fontsize=18)
        # pervasive range
        """ax5.fill_between([10**(2), 10**(4)], [100, 100],
                        color="#a9d5b2", edgecolor='black', alpha=0.1)"""
        ax5.fill_between([10**(1), 10**(4)], [100, 100],
                        cmap="Purples", alpha=0.2, step='mid')
        # channelized range
        ax5.fill_between([10**(4), 10**(6)], [100, 100],
                        color="#ca638f", alpha=0.1)
        # regional range
        """ax5.fill_between([10**(2.7-0.5), 10**(2.7+0.5)], [100, 100],
                        color="#c5c5c5", edgecolor='black', alpha=0.5)"""
        ax5.tick_params(axis='both', labelsize=20)
        ax5.vlines(10**4, 0, 100, linestyles='--',
                color='black', linewidth=4)
        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
        # ax5.annotate("Regional range", (10**(2.7+0.5-0.3), 25), fontsize=16, bbox=props, rotation=90)
        ax5.annotate("Channelized", (10**(4+0.2), 5),
                    fontsize=16, bbox=props)
        ax5.annotate("Pervasive", (10**(1+0.2), 5),
                    fontsize=16, bbox=props)

        if img_save is True:
            os.makedirs(
                f'{self.mainfolder}/img_{self.filename}/{subfolder}', exist_ok=True)
            plt.savefig(f'{self.mainfolder}/img_{self.filename}/{subfolder}/fluid_flux_summary.png',
                        transparent=False, facecolor='white')
            # plt.pause(2)
            plt.close()
        else:
            plt.show()

    def porosity_plot(self, rock_tag, img_save=False, gif_save=False):

        if gif_save is True:
            img_save = True

        group_key = self.rockdic[rock_tag].group_key

        phases = self.rockdic[rock_tag].phases
        phase_data = self.rockdic[rock_tag].phase_data['df_N']
        # Clean up dataframe to be used by mineral names and no NaN
        phase_data.columns = phases
        phase_data[np.isnan(phase_data) == True] = 0

        ts = self.rockdic[rock_tag].temperature
        ps = self.rockdic[rock_tag].pressure
        depth = self.rockdic[rock_tag].depth

        # Filter method of data array
        unfiltered_porosity = self.rockdic[rock_tag].porosity
        extraction_boolean = self.rockdic[rock_tag].extraction_boolean

        if False in extraction_boolean:
            filtered_porosity = np.ma.masked_array(
                unfiltered_porosity, extraction_boolean)
        else:
            filtered_porosity = np.ma.masked_array(
                np.zeros(len(unfiltered_porosity)), extraction_boolean)

        if len(np.unique(self.rockdic[rock_tag].frac_bool)) == 1 and np.unique(self.rockdic[rock_tag].frac_bool)[-1] == 0:
            filtered_porosity = np.zeros(
                len(self.rockdic[rock_tag].frac_bool))

        subfolder = 'porosity'

        if img_save is True:

            print("\n\nGenerating fluid-filled porosity images. Please wait...")
            # Start progress bar
            k = 0
            kk = len(ts)
            progress(int(k/kk)*100)
            os.makedirs(
                f'{self.mainfolder}/img_{self.filename}', exist_ok=True)

            for i in range(len(phase_data.index)):
                os.makedirs(
                    f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}', exist_ok=True)
                # -------------------------------------------------------------------
                fig = plt.figure(constrained_layout=False,
                                 facecolor='0.9', figsize=(9, 9))
                gs = fig.add_gridspec(nrows=3, ncols=3, left=0.15,
                                      right=0.95, hspace=0.6, wspace=0.5)

                ax5 = fig.add_subplot(gs[:, :])
                ax5.plot(unfiltered_porosity[:i+1]*100,
                         depth[:i+1]/1000, 'd--', color='green')
                if False in filtered_porosity.mask[:i+1]:
                    # plot 4
                    ax5.plot(
                        filtered_porosity[:i+1]*100, depth[:i+1]/1000, 'd--', color='black', alpha=0.7)
                    ax5.plot(filtered_porosity[i:i+1]*100, depth[i:i+1]/1000, 'd',
                             color='#7fffd4', markersize=8, markeredgecolor='black')
                else:
                    # plot 4
                    ax5.plot(np.ones(i+1)*1e-30, depth[:i+1]/1000, 'd')

                # ax5 figure features
                ax5.set_title('T:{:.2f} °C P:{:.2f} GPa'.format(
                    ts[i], ps[i]/10000), fontsize=18)
                ax5.set_xlim(0.01, 10)
                ax5.set_xscale('log')
                ax5.set_ylim(100, 0)
                ax5.set_xlabel("Fluid-filled porosity [$vol.%$]", fontsize=18)
                ax5.set_ylabel("Depth [km]", fontsize=18)
                # ax5.set_xlim(10**0, 10**6)
                # ax5.fill_between([10**2.7 - 0.5, 10**2.7 + 0.5], [100, 100],
                #                  color="#c5c5c5", edgecolor='black', hatch="/", alpha=0.4)
                ax5.tick_params(axis='both', labelsize=20)
                # ax5.annotate("Ague 2003 Comp", (10**2.7, 15), fontsize=16)

            # Finishing plot - save
                plt.savefig(f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}/img_{i}.png',
                            transparent=False, facecolor='white')
                # plt.pause(2)
                plt.close()

                # tail of progress bar
                k += 1
                ic = int(np.ceil(k/kk*100))
                print("=====Progress=====")
                progress(ic)

        else:
            fig = plt.figure(constrained_layout=False,
                             facecolor='0.9', figsize=(9, 9))
            gs = fig.add_gridspec(nrows=3, ncols=3, left=0.15,
                                  right=0.95, hspace=0.6, wspace=0.5)
            ax5 = fig.add_subplot(gs[:, :])

            # looping the array
            for i in range(len(phase_data.index)):
                ax5.plot(unfiltered_porosity[:i+1]*100,
                         depth[:i+1]/1000, 'd--', color='green')
                if False in filtered_porosity.mask[:i+1]:
                    # plot 4
                    ax5.plot(
                        filtered_porosity[:i+1]*100, depth[:i+1]/1000, 'd--', color='black', alpha=0.7)
                    ax5.plot(filtered_porosity[i:i+1]*100, depth[i:i+1]/1000, 'd',
                             color='#7fffd4', markersize=8, markeredgecolor='black')
                else:
                    # plot 4
                    ax5.plot(np.ones(i+1)*1e-30, depth[:i+1]/1000, 'd')

            # ax5 figure features
            ax5.set_title('T:{:.2f} °C P:{:.2f} GPa'.format(
                ts[i], ps[i]/10000), fontsize=18)
            ax5.set_xlim(0.01, 10)
            ax5.set_xscale('log')
            ax5.set_ylim(100, 0)
            ax5.set_xlabel("Fluid-filled porosity [$vol.%$]", fontsize=18)
            ax5.set_ylabel("Depth [km]", fontsize=18)
            # ax5.set_xlim(10**0, 10**6)
            # ax5.fill_between([10**2.7 - 0.5, 10**2.7 + 0.5], [100, 100],
            #                  color="#c5c5c5", edgecolor='black', hatch="/", alpha=0.4)
            ax5.tick_params(axis='both', labelsize=20)
            # ax5.annotate("Ague 2003 Comp", (10**2.7, 15), fontsize=16)
            plt.show()

        # reading images and put it to GIF
        # reading images and put it to GIF
        if gif_save is True:
            # call gif function
            create_gif(phase_data, self.mainfolder, self.filename,
                       group_key, subfolder=subfolder)

    def release_fluid_volume_plot(self, rock_tag, img_save=False, gif_save=False):

        if gif_save is True:
            img_save = True

        ts = self.rockdic[rock_tag].temperature
        group_key = self.rockdic[rock_tag].group_key
        phases = self.rockdic[rock_tag].phases
        phase_data = self.rockdic[rock_tag].phase_data['df_N']
        phase_data.columns = phases
        phase_data[np.isnan(phase_data) == True] = 0

        depth = self.rockdic[rock_tag].depth

        # Filter method of data array
        unfiltered_v_fluid_extr = self.rockdic[rock_tag].extracted_fluid_volume
        permea = self.rockdic[rock_tag].permeability
        extraction_boolean = self.rockdic[rock_tag].extraction_boolean

        if False in extraction_boolean:
            filtered_v_fluid_extr = np.ma.masked_array(
                unfiltered_v_fluid_extr, extraction_boolean)
            filtered_permeability = np.ma.masked_array(
                permea, extraction_boolean)
        else:
            filtered_v_fluid_extr = np.ma.masked_array(
                np.zeros(len(unfiltered_v_fluid_extr)), extraction_boolean)
            filtered_permeability = np.ma.masked_array(
                np.zeros(len(extraction_boolean)), extraction_boolean)

        if len(np.unique(self.rockdic[rock_tag].frac_bool)) == 1 and np.unique(self.rockdic[rock_tag].frac_bool)[-1] == 0:
            filtered_v_fluid_extr = np.zeros(
                len(self.rockdic[rock_tag].frac_bool))
            filtered_permeability = np.zeros(
                len(self.rockdic[rock_tag].frac_bool))

        subfolder = 'released_volume'

        if img_save is True:

            print("\n\nGenerating release fluid volume images. Please wait...")
            # Start progress bar
            k = 0
            kk = len(ts)
            progress(int(k/kk)*100)

            os.makedirs(
                f'{self.mainfolder}/img_{self.filename}', exist_ok=True)

            for i in range(len(phase_data.index)):
                os.makedirs(
                    f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}', exist_ok=True)
                # -------------------------------------------------------------------
                fig = plt.figure(constrained_layout=False,
                                 facecolor='0.9', figsize=(9, 9))
                gs = fig.add_gridspec(nrows=3, ncols=3, left=0.15,
                                      right=0.95, hspace=0.6, wspace=0.5)

                ax5 = fig.add_subplot(gs[:, :])
                ax5.plot(unfiltered_v_fluid_extr[:i+1],
                         depth[:i+1]/1000, 'd--', color='green')
                if False in filtered_permeability.mask[:i+1]:
                    ax5.plot(filtered_v_fluid_extr[:i+1] /
                             1_000_000, depth[:i+1]/1000, 'd--')
                    ax5.plot(filtered_v_fluid_extr[:i+1]/1_000_000,
                             depth[:i+1]/1000, 'd--', color='black', alpha=0.5)
                    ax5.plot(filtered_v_fluid_extr[i:i+1]/1_000_000, depth[i:i+1]/1000,
                             'd--', color='#7fffd4', markersize=8, markeredgecolor='black')

                # ax5 figure features
                ax5.set_title("Extracted fluid volume")
                ax5.set_xlim(0.001, 100)
                ax5.set_xscale('log')
                ax5.set_ylim(100, 0)
                # TODO ccm because it is not the geometry scale
                ax5.set_xlabel("Volume extracted [$m^{3}$]", fontsize=12)
                ax5.set_ylabel("Depth [km]", fontsize=12)
                ax5.tick_params(axis='both', labelsize=14)

            # Finishing plot - save
                plt.savefig(f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}/img_{i}.png',
                            transparent=False, facecolor='white')
                # plt.pause(2)
                plt.close()

                # tail of progress bar
                k += 1
                ic = int(np.ceil(k/kk*100))
                print("=====Progress=====")
                progress(ic)

        else:
            fig = plt.figure(constrained_layout=False,
                             facecolor='0.9', figsize=(9, 9))
            gs = fig.add_gridspec(nrows=3, ncols=3, left=0.15,
                                  right=0.95, hspace=0.6, wspace=0.5)
            ax5 = fig.add_subplot(gs[:, :])
            for i in range(len(phase_data.index)):
                ax5.plot(unfiltered_v_fluid_extr[:i+1],
                         depth[:i+1]/1000, 'd--', color='green')
                if False in filtered_permeability.mask[:i+1]:
                    ax5.plot(filtered_v_fluid_extr[:i+1] /
                             1_000_000, depth[:i+1]/1000, 'd--')
                    ax5.plot(filtered_v_fluid_extr[:i+1]/1_000_000,
                             depth[:i+1]/1000, 'd--', color='black', alpha=0.5)
                    ax5.plot(filtered_v_fluid_extr[i:i+1]/1_000_000, depth[i:i+1]/1000,
                             'd--', color='#7fffd4', markersize=8, markeredgecolor='black')

            # ax5 figure features
            ax5.set_title("Extracted fluid volume")
            ax5.set_xlim(0.001, 100)
            ax5.set_xscale('log')
            ax5.set_ylim(100, 0)
            # TODO ccm because it is not the geometry scale
            ax5.set_xlabel("Volume extracted [$m^{3}$]", fontsize=12)
            ax5.set_ylabel("Depth [km]", fontsize=12)
            ax5.tick_params(axis='both', labelsize=14)
            plt.show()

        # reading images and put it to GIF
        if gif_save is True:
            # call gif function
            create_gif(phase_data, self.mainfolder, self.filename,
                       group_key, subfolder=subfolder)

    def phases_stack_plot(self, rock_tag, img_save=False, val_tag=False, transparent=False, fluid_porosity=True):
        group_key = self.rockdic[rock_tag].group_key
        subfolder = 'stack_plot'

        # XMT naming and coloring
        database = self.rockdic[rock_tag].database
        phases2 = self.rockdic[rock_tag].phases2
        legend_phases, color_set = phases_and_colors_XMT(database, phases2)

        phases = self.rockdic[rock_tag].phases
        phase_data = self.rockdic[rock_tag].phase_data
        frac_bool = self.rockdic[rock_tag].frac_bool
        garnet = self.rockdic[rock_tag].garnet[rock_tag]
        garnet_bool = self.rockdic[rock_tag].garnets_bools[rock_tag]
        # Input for the variable of interest
        if val_tag is False:
            tag_in = input(
                "Please provide what you want to convert to a stack. ['vol%', 'volume', 'wt%', 'wt']")
            tag = 'df_'+tag_in
        else:
            tag_in = val_tag
            if val_tag in ['vol%', 'volume', 'wt%', 'wt']:
                tag = 'df_'+tag_in
                pass
            else:
                print("Try again and select a proper value for stack plot input")
                quit()

        # compile data for plotting
        system_vol_pre = self.rockdic[rock_tag].systemVolPre
        system_vol_post = self.rockdic[rock_tag].systemVolPost
        st_fluid_before = self.rockdic[rock_tag].fluidbefore

        # max_vol = max(system_vol_pre)
        temperatures = self.rockdic[rock_tag].temperature
        pressures = self.rockdic[rock_tag].pressure
        line = np.arange(1, len(temperatures)+1, 1)

        y = phase_data[tag].fillna(value=0)
        y.columns = legend_phases
        y = y.T

        # cleaning dataframe from multiple phase names - combine rows into one
        if len(legend_phases) == len(np.unique(legend_phases)):
            pass
        else:
            y, legend_phases, color_set = clean_frame(y, legend_phases, color_set)

        # resort the dataframe, legend_phases and color_set so that "Water" is always last
        if 'Water' in legend_phases:
            y, legend_phases, color_set = resort_frame(y, legend_phases, color_set)
        else:
            pass

        # drop rows with all zeros and delete the corresponding legend_phases and color_set
        new_list = y.loc[(y != 0).any(axis=1)].index.tolist()
        new_color_set = []
        for phase in y.index:
            if phase in new_list:
                indi = y.index.tolist().index(phase)
                new_color_set.append(color_set[indi])
        y = y.loc[(y != 0).any(axis=1)]
        legend_phases = new_list
        color_set = new_color_set

        # "Serp_Atg', 'Br_Br', 'Olivine', 'Chlorite', 'Magnetite', 'Water"
        # color_set = ["#91D7EA", "#DB7012", "#4E459B", "#B32026", "#FAEA64", "#2E83D0"]
        # color_set = ["#91D7EA", "#DB7012", "#FAEA64", "#4E459B", "#E724C5",  "#969696", "#2E83D0"]
        # plt.rcParams['font.family'] = 'sans-serif'

        # add metastable garnet to the dataframe
        if 'Garnet' in y.index:
            y.loc['Garnet'][np.isnan(y.loc['Garnet'])] = 0
            kk = 0
            garnet_in = False
            for k, xval in enumerate(y.loc['Garnet']):
                if xval == 0 and garnet_in is True:
                    y.loc['Garnet',k] = np.cumsum(garnet.volume[:1+kk])[-1]
                elif xval == 0:
                    pass
                else:
                    garnet_in = True
                    y.loc['Garnet',k] = np.cumsum(garnet.volume[:1+kk])[-1]
                    kk += 1
            y.loc['Garnet'][y.loc['Garnet']==0] = np.nan

        # nomrmalize the volume data to starting volume
        y = y/y.sum()[0]

        # Test if temperature has a negative slope
        if np.any(np.diff(temperatures)[np.diff(temperatures)<0]) is True:


            # plot
            plt.rc('axes', labelsize=16)
            plt.rc('xtick', labelsize=12)  # fontsize of the x tick labels
            plt.rc('ytick', labelsize=12)  # fontsize of the y tick labels
            plt.figure(111, figsize=(10,6), facecolor=None, )

            ax2 = host_subplot(111)
            #fig.suptitle('Phase changes for P-T-t')
            # main plot
            if transparent is True:
                ax2.stackplot(line, y, labels=legend_phases,
                        colors=color_set, alpha=.35, edgecolor='black')
            else:
                ax2.stackplot(line, y, labels=legend_phases,
                                    colors=color_set, alpha=.7, edgecolor='black')

            # define the legend with its handles, labels and its position
            handles, labels = ax2.get_legend_handles_labels()
            legend = ax2.legend(
                handles[::-1], labels[::-1], loc='right',
                borderaxespad=0.1, title="Stable phases", fontsize=14
                )

            # tick stepping
            con = round(10/len(line), 2)
            con = round(len(line)*con)
            step = round(len(line)/con)
            if step == 0:
                step = 1
            if step == 4:
                step=5

            # if np.any(np.diff(temperatures)[np.diff(temperatures)<0]) is True:
            # define x-axis bottom (temperature) tick labels
            ax2.set_xticks(line[::step])
            if len(line) < len(temperatures):
                selected_temp = temperatures[list(line-1)]
                f_arr = np.around(selected_temp[::step], 0)
                ax2.set_xticklabels(f_arr.astype(int))
            else:
                f_arr = np.around(temperatures[::step], 0)
                ax2.set_xticklabels(f_arr.astype(int))
            ax2.xaxis.set_minor_locator(AutoMinorLocator())

            # second y axis
            # free fluid content
            fluid_porosity_color = "#4750d4"
            # fluid content as vol% of total system volume
            y2 = (st_fluid_before)/system_vol_pre*100
            # extraction steps
            # NOTE extraction marker boolean
            mark_extr = extraction_boolean = self.rockdic[rock_tag].extraction_boolean
            mark_extr = np.array(system_vol_pre-system_vol_post, dtype='bool')
            if fluid_porosity is True:
                twin1 = ax2.twinx()
                twin1.plot(line, y2, 'o--', c=fluid_porosity_color, linewidth=2, markeredgecolor='black')
                twin1.set_ylabel("Vol% of fluid-filled porosity", color=fluid_porosity_color, fontsize=15, weight='bold')
                twin1.set_ymargin(0)

                if len(frac_bool) > 0:
                    if 1 in frac_bool or 2 in frac_bool or 3 in frac_bool or 10 in frac_bool:
                        extension_bool = np.isin(frac_bool, 1)
                        extend_shear_bool = np.isin(frac_bool, 2)
                        compress_shear_bool = np.isin(frac_bool, 3)
                        ten_bool = np.isin(frac_bool, 10)
                        twin1.plot(line[extension_bool], y2[extension_bool], 'Dr')
                        twin1.plot(line[extend_shear_bool], y2[extend_shear_bool], 'Dg')
                        twin1.plot(line[compress_shear_bool], y2[compress_shear_bool], 'Db')
                        twin1.plot(line[ten_bool], y2[ten_bool], 'D', c='violet')
                else:
                    pass


            # define x-axis top (pressure) tick labels
            ax3 = ax2.twin()
            ax3.set_xticks(line[::step])

            if len(line) < len(pressures):
                selected_pres = pressures[list(line-1)]
                ax2.set_xticklabels(np.around(selected_pres[::step], 1))
            else:
                ax3.set_xticklabels(np.around(np.array(pressures[::step])/10000, 1))

            # labeling and style adjustments
            plt.subplots_adjust(right=0.75)
            pressures = list(pressures)
            temperatures = list(temperatures)
            peakp = pressures.index(max(pressures))
            peakt = temperatures.index(max(temperatures))
            if len(temperatures) > peakp+1:
                ax2.axvline(x=line[peakp+1], linewidth=1.5,
                            linestyle='--', color='black')
                ax2.axvline(x=line[peakp], linewidth=1.5,
                            linestyle='--', color='black')
                ax2.axvline(x=line[peakt], linewidth=1.5,
                            linestyle='--', color='red')
            if tag[3:] == 'volume[ccm]':
                ax2.set_ylabel("Relative volume")
            else:
                ax2.set_ylabel(tag[3:])
            ax2.set_xlabel('Temperature [°C]')
            ax3.set_xlabel('Pressure [GPa]')
            ax3.axis["right"].major_ticklabels.set_visible(False)
            ax3.axis["right"].major_ticks.set_visible(False)

            if fluid_porosity is True:
                twin1.tick_params(colors=fluid_porosity_color)
                twin1.set_yticklabels(twin1.get_yticks(), weight='heavy', size=14)
                ymax= np.round(max(y2),2)
                twin1.set_ylim(0, np.round(ymax+ymax*0.05,2))

            ax2.hlines(1, -2, max(line)+5, 'black', linewidth=1.5, linestyle='--')
            ax2.set_xlim(0,max(line)+1)

            ax2.set_facecolor("#fcfcfc")
            ax2.set_alpha(0.5)

            # Fix the first y axis to a range from 0 to 1.1
            ax2.set_ylim(0, 1.1)

            # fit figure and legend in figsize
            plt.tight_layout(rect=[0, 0, 1, 1])
            legend = ax2.legend(bbox_to_anchor=(1.4, 0.9))

        else:
            # plot
            plt.rc('axes', labelsize=16)
            plt.rc('xtick', labelsize=12)  # fontsize of the x tick labels
            plt.rc('ytick', labelsize=12)  # fontsize of the y tick labels
            plt.figure(111, figsize=(10,6), facecolor=None, )

            ax2 = host_subplot(111)
            #fig.suptitle('Phase changes for P-T-t')
            # main plot
            if transparent is True:
                ax2.stackplot(temperatures, y, labels=legend_phases,
                        colors=color_set, alpha=.35, edgecolor='black')
            else:
                ax2.stackplot(temperatures, y, labels=legend_phases,
                                    colors=color_set, alpha=.7, edgecolor='black')

            # second y axis
            # free fluid content
            fluid_porosity_color = "#4750d4"
            # fluid content as vol% of total system volume
            y2 = (st_fluid_before)/system_vol_pre*100
            y2 = np.round(y2,1)
            # extraction steps
            if fluid_porosity is True:
                twin1 = ax2.twinx()
                twin1.plot(temperatures, y2, 'o--', c=fluid_porosity_color, linewidth=2, markeredgecolor='black')
                twin1.set_ylabel("Vol% of fluid-filled porosity", color=fluid_porosity_color, fontsize=15, weight='bold')
                twin1.set_ymargin(0)

                if len(frac_bool) > 0:
                    if 1 in frac_bool or 2 in frac_bool or 3 in frac_bool or 10 in frac_bool:
                        extension_bool = np.isin(frac_bool, 1)
                        extend_shear_bool = np.isin(frac_bool, 2)
                        compress_shear_bool = np.isin(frac_bool, 3)
                        ten_bool = np.isin(frac_bool, 10)
                        twin1.plot(temperatures[extension_bool], y2[extension_bool], 'Dr')
                        twin1.plot(temperatures[extend_shear_bool], y2[extend_shear_bool], 'Dg')
                        twin1.plot(temperatures[compress_shear_bool], y2[compress_shear_bool], 'Db')
                        twin1.plot(temperatures[ten_bool], y2[ten_bool], 'D', c='violet')
                else:
                    pass

            # define x-axis top (pressure) tick labels
            ax3 = ax2.twiny()
            # plot zeros of length of pressure array on ax3 which will be invisible
            ax3.plot(np.zeros(len(pressures)), pressures, alpha=0)

            # labeling and style adjustments
            if tag[3:] == 'volume[ccm]':
                ax2.set_ylabel("Relative volume")
            else:
                ax2.set_ylabel(tag[3:])
            ax2.set_xlabel('Temperature [°C]')
            ax3.set_xlabel('Pressure [GPa]')
            ax3.axis["right"].major_ticklabels.set_visible(False)
            ax3.axis["right"].major_ticks.set_visible(False)

            if fluid_porosity is True:
                twin1.tick_params(colors=fluid_porosity_color)
                tick_list = list(np.round(twin1.get_yticks(),2))
                twin1.set_yticks(ticks=tick_list, labels=tick_list, fontweight='heavy', fontsize=14)
                ymax= max(y2)
                twin1.set_ylim(0, ymax+np.round(ymax*0.01, 1))

            ax2.set_facecolor("#fcfcfc")
            ax2.set_alpha(0.5)
            ax2.hlines(1, min(temperatures), max(temperatures), 'black', linewidth=1.5, linestyle='--')

            # Fix the first y axis to a range from 0 to 1.1
            ax2.set_ylim(0, 1.1)

            # fit figure and legend in figsize
            plt.tight_layout(rect=[0, 0, 0.79, 1])
            # define the legend with its handles, labels and its position
            handles, labels = ax2.get_legend_handles_labels()
            legend = ax2.legend(
                handles[::-1], labels[::-1], bbox_to_anchor=(1.12, 1.01),
                borderaxespad=0.1, title="Stable phases", fontsize=14
                )

        # save image
        if img_save is True:
            os.makedirs(
                    f'{self.mainfolder}/img_{self.filename}/{subfolder}', exist_ok=True)
            if transparent is True:
                plt.savefig(f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}_stack_plot_transparent.png',
                                                transparent=True)
            else:
                plt.savefig(f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}_stack_plot.png',
                                transparent=False)

        else:
            plt.show()
        plt.clf()
        plt.close()

    def binary_plot(self, rock_tag, img_save=False):

        data = self.rockdic[rock_tag]
        attributes = [field.name for field in fields(data)]

        # TODO - Add exception when user aborts input or gives no input - Needs a default mode (x,y)
        x_ax_label = input("Please provide the x-axis data...")
        y_ax_label = input("Please provide the y-axis data...")

        whitelist_attributes = ['temperature', 'pressure', 'depth', 'time', 'systemVolPre',
                                'systemVolPost', 'permeability', 'extracted_fluid_volume', 'porosity', 'time_int_flux2']
        whitelist_phase_data = ['N', 'vol%', 'volume[ccm]', 'wt%', 'wt[g]']

        if x_ax_label in attributes and x_ax_label in whitelist_attributes:
            if x_ax_label == 'temperature':
                xdata = data.temperature
            if x_ax_label == 'pressure':
                xdata = data.pressure
            if x_ax_label == 'depth':
                xdata = data.depth
            if x_ax_label == 'time':
                xdata = data.time
            if x_ax_label == 'systemVolPost':
                xdata = data.systemVolPost
            if x_ax_label == 'permeability':
                xdata = data.permeability
            if x_ax_label == 'extracted_fluid_volume':
                xdata = data.extracted_fluid_volume
            if x_ax_label == 'porosity':
                xdata = data.porosity
            if x_ax_label == 'time_int_flux2':
                xdata = data.time_int_flux2
        elif x_ax_label in whitelist_phase_data:
            sel_d = data.phase_data[f'df_{x_ax_label}']
            xdata = input(
                f"Your x data query requires a phase input from the following list:\n{sel_d.columns}")
            xdata = data.phase_data[f'df_{x_ax_label}'][xdata]
        else:
            print("Selected data for x label is not available.")
            x_ax_label = False

        if y_ax_label in attributes and y_ax_label in whitelist_attributes:
            if y_ax_label == 'temperature':
                ydata = data.temperature
            if y_ax_label == 'pressure':
                ydata = data.pressure
            if y_ax_label == 'depth':
                ydata = data.depth
            if y_ax_label == 'systemVolPost':
                ydata = data.systemVolPost
            if y_ax_label == 'permeability':
                ydata = data.permeability
            if y_ax_label == 'extracted_fluid_volume':
                ydata = data.extracted_fluid_volume
            if y_ax_label == 'porosity':
                ydata = data.porosity
            if y_ax_label == 'time_int_flux2':
                ydata = data.time_int_flux2
        elif y_ax_label in whitelist_phase_data:
            sel_d = data.phase_data[f'df_{y_ax_label}']
            ydata = input(
                f"Your x data query requires a phase input from the following list:\n{sel_d.columns}")
            ydata = data.phase_data[f'df_{y_ax_label}'][ydata]
        else:
            print("Selected data for x label is not available.")
            y_ax_label = False

        if x_ax_label is False or y_ax_label is False:
            print("No binary plot generated because no legit user input")
            pass
        else:
            plt.scatter(xdata, ydata, s=24, c='#7fffd4', edgecolor='black')
            plt.xlabel(x_ax_label)
            plt.ylabel(y_ax_label)
            # saving option from user input
            if img_save is True:
                plt.savefig(Path(self.mainfolder /
                            f"{self.rock_key}_binary_plot.png"), dpi=300)
            else:
                plt.show()

    def oxygen_isotopes(self, rock_tag, img_save=False):
        """Plotting function for the modelled oxygen isotope data

        Args:
            rock_tag (str): the name of the rock such as "rock0", "rock1", ...
            img_save (bool, optional): Optional argument to save the plot as an image to the directory. Defaults to False.
        """

        group_key = self.rockdic[rock_tag].group_key
        subfolder = 'oxygen_plot'

        # Read the oxygen data from the dictionary
        summary = self.rockdic[rock_tag].oxygen_data.T.describe()
        # oxyframe = Merge_phase_group(self.rockdic[rock_tag].oxygen_data)
        oxyframe = self.rockdic[rock_tag].oxygen_data.T
        oxyframe.columns = self.rockdic[rock_tag].temperature

        # XMT naming and coloring
        database = self.rockdic[rock_tag].database
        phases2 = oxyframe.index
        legend_phases, color_set = phases_and_colors_XMT(database, phases2)

        # Plotting routine
        for item in ['"H"', '"Si"', '"Ti"', '"Al"']:
            if item in legend_phases:
                indi = legend_phases.index(item)
                del legend_phases[indi]

        oxyframe.index = legend_phases

        # cleaning dataframe from multiple phase names - combine rows into one
        if len(legend_phases) == len(np.unique(legend_phases)):
            pass
        else:
            oxyframe, legend_phases, color_set = clean_frame(oxyframe, legend_phases, color_set)

        # manual colorset for the plot
        # "Serp_Atg', 'Br_Br', 'Olivine', 'Chlorite', 'Magnetite', 'Water"
        # color_set = ["#91D7EA", "#DB7012", "#4E459B", "#B32026", "#FAEA64", "#2E83D0"]
        # color_set = ["#91D7EA", "#DB7012", "#FAEA64", "#2E83D0", "#4E459B", "#E724C5",  "#969696"]
        # plt.rcParams['font.family'] = 'arial'
        fig, ax211 = plt.subplots(1, 1, figsize=(8, 5))
        for t, phase in enumerate(list(oxyframe.index)):
            print(t)
            ax211.plot(oxyframe.columns, oxyframe.loc[phase], '--d',
                       color=color_set[t], linewidth=0.7, markeredgecolor='black', markersize=9)

        ax211.plot(self.rockdic[rock_tag].temperature,
                   self.rockdic[rock_tag].bulk_deltao_post, '-', c='black')
        # ax211.plot(self.rockdic[rock_tag].temperature,
        #            self.rockdic[rock_tag].bulk_deltao_pre, '-.', c='black')
        # legend_list = list(oxyframe.index) + ["pre bulk", "post bulk"]
        legend_list = list(oxyframe.index) + ["bulk"]
        ax211.legend(legend_list, bbox_to_anchor=(1.28, 0.9))
        min_min = min(summary.loc['min'])
        max_max = max(summary.loc['max'])
        min_val = min_min - (max_max-min_min)*0.05
        max_val = max_max + (max_max-min_min)*0.05
        ax211.set_ylim(min_val, max_val)
        ax211.set_ylabel("$\delta^{18}$O [‰ vs. VSMOW]")
        ax211.set_xlabel("Temperature [°C]")
        plt.subplots_adjust(right=0.8)

        if img_save is True:
            os.makedirs(
                    f'{self.mainfolder}/img_{self.filename}/{subfolder}', exist_ok=True)
            plt.savefig(f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}_oxygen_isotope_plot.png',
                            transparent=False, facecolor='white')
        else:
            plt.show()
        plt.clf()
        plt.close()

    # function to plot a column consisting of all rocks showing the fluid content
    def fluid_content(self, img_save=False):
        """Plotting function for the fluid content of all rocks in the model and the glaucophane content in two subplots.

        Args:
            img_save (bool, optional): Optional argument to save the plot as an image to the directory. Defaults to False.
        """
        # read the temperature data from the first rock in the dictionary
        ts = self.rockdic['rock000'].temperature
        # set the color palette for the plot based on the number of rocks
        color_palette = sns.color_palette("viridis", len(self.rockdic.keys()))

        # plotting routine
        fig, (ax1, ax2, ax3, ax4, ax5, ax6) = plt.subplots(6, 1, figsize=(10, 9))
        for i, rock in enumerate(self.rockdic.keys()):
            # skip the first rock when i == 0
            if i == 0:
                pass
            else:
                # read the database from the rockdic and store it in a variable
                database = self.rockdic[rock].database
                phases2 = list(self.rockdic[rock].phase_data['df_vol%'].columns)
                legend_phases, color_set = phases_and_colors_XMT(database, phases2)

                # read the vol% data to a variable and replace the NaN values with 0 and transpose the dataframe
                y = self.rockdic[rock].phase_data['df_vol%'].fillna(value=0)
                y.columns = legend_phases
                y = y.T

                # clean the dataframe from multiple phase names - combine rows into one
                if len(legend_phases) == len(np.unique(legend_phases)):
                    pass
                else:
                    y, legend_phases, color_set = clean_frame(y, legend_phases, color_set)

                # drop rows with all zeros
                y = y.loc[(y != 0).any(axis=1)]
                legend_phases = list(y.index)

                # define the fluid content from the dataframe based on 'Water' and plot it
                if 'Water' in legend_phases:
                    fluid_content = y.loc['Water']
                    ax2.plot(ts, fluid_content,
                            '-', c=color_palette[i], linewidth=1.2, markeredgecolor='black')
                    # y axes for ax2 subplot with 10 ticks
                    start, end = ax2.get_ylim()
                    ax2.yaxis.set_ticks(np.round(np.linspace(start, end+end*0.02, 5),1))
                    # plot horizontal line at y = 0
                    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1.0, alpha=0.6)

                # if -Glaucophane- is in y plot it as a dashed line
                if 'Glaucophane' in legend_phases:
                    glaucophane_content = y.loc['Glaucophane']
                    ax3.plot(ts, glaucophane_content,
                            '--', c=color_palette[i], linewidth=1.2, markeredgecolor='black')
                    # y axes for ax3 subplot with 10 ticks
                    start, end = ax3.get_ylim()
                    ax3.yaxis.set_ticks(np.round(np.linspace(start, end+end*0.02, 5),1))
                    # plot horizontal line at y = 0
                    ax3.axhline(y=0, color='black', linestyle='-', linewidth=1.0, alpha=0.6)


                # if -Amphibole- is in y plot it as another dashed line with a different color
                if 'Amphibole' in legend_phases:
                    amphibole_content = y.loc['Amphibole']
                    ax4.plot(ts, amphibole_content,
                            '-.', c=color_palette[i], linewidth=1.2, markeredgecolor='black')
                    # y axes for ax4 subplot with 10 ticks
                    start, end = ax4.get_ylim()
                    ax4.yaxis.set_ticks(np.round(np.linspace(start, end+end*0.02, 5),1))
                    # plot horizontal line at y = 0
                    ax4.axhline(y=0, color='black', linestyle='-', linewidth=1.0, alpha=0.6)

                    # lopp through hydrous phases, test if they are in the legend_phases, get the indices and sum their values from y
                    hydrous_phases = ['Muscovite', 'Amphibole', 'Glaucophane', 'Lawsonite', 'Talc', 'Phengite']
                    inid_list = []
                    for item in hydrous_phases:
                        if item in legend_phases:
                            inid_list.append(legend_phases.index(item))
                    hydrous_content = y.iloc[inid_list].sum(axis=0)

                    # plot the hydrous content as a dashed line with a different color
                    ax4.plot(ts, hydrous_content,
                            '--', c=color_palette[i], linewidth=1.2, markeredgecolor='black')

                """# if -Omphacite- is in y plot it as a dashed line
                if 'Omphacite' in legend_phases:
                    omphacite_content = y.loc['Omphacite']
                    ax5.plot(ts, omphacite_content,
                            '--', c=color_palette[i], linewidth=1.2, markeredgecolor='black')
                    # y ticks for the omphacite subplot with 10 ticks based on max and min of omphacite content
                    # y axes for ax5 subplot with 10 ticks
                    start, end = ax5.get_ylim()
                    ax5.yaxis.set_ticks(np.round(np.linspace(start, end+end*0.02, 5),1))
                    # plot horizontal line at y = 0
                    ax5.axhline(y=0, color='black', linestyle='-', linewidth=1.0, alpha=0.6)"""

                # if -Lawsonite- is in y plot it as a dashed line
                if 'Lawsonite' in legend_phases:
                    omphacite_content = y.loc['Lawsonite']
                    ax5.plot(ts, omphacite_content,
                            '--', c=color_palette[i], linewidth=1.2, markeredgecolor='black')
                    # y ticks for the omphacite subplot with 10 ticks based on max and min of omphacite content
                    # y axes for ax5 subplot with 10 ticks
                    start, end = ax5.get_ylim()
                    ax5.yaxis.set_ticks(np.round(np.linspace(start, end+end*0.02, 5),1))
                    # plot horizontal line at y = 0
                    ax5.axhline(y=0, color='black', linestyle='-', linewidth=1.0, alpha=0.6)


        # plot the water content of rock000 in ax4
        if 'fluid' in self.rockdic['rock000'].phase_data['df_vol%'].columns:
            water_fluid_content = self.rockdic['rock000'].phase_data['df_vol%']['fluid']
            ax6.plot(ts, water_fluid_content,
                    '-', c='black', linewidth=2, markeredgecolor='black')
            # y axes for ax6 subplot with 10 ticks
            start, end = ax6.get_ylim()
            ax6.yaxis.set_ticks(np.round(np.linspace(start, end+end*0.02, 5),1))
            # plot horizontal line at y = 0
            ax6.axhline(y=0, color='black', linestyle='-', linewidth=1.0, alpha=0.6)

        # plot the water that is release from the least rock in rockdic.keys() in ax1
        top_rock_name = list(self.rockdic.keys())[-1]
        if 'fluid' in self.rockdic[top_rock_name].phase_data['df_vol%'].columns:
            water_release = self.rockdic[top_rock_name].extracted_fluid_volume/1000000 # in m3
            # water release as cumulative sum
            water_release = np.cumsum(water_release)
            ax1.plot(ts, water_release,
                    '-', c='black', linewidth=2, markeredgecolor='black')


        # Set the x-label at lower most subplot represetning all plots
        ax6.set_xlabel("Prograde P-T")

        # Title each subplot
        ax1.set_title("Water release at top [$m^{3}$]")
        ax2.set_title("Fluid content [vol.%]")
        ax3.set_title("Glaucophane [vol.%]")
        ax4.set_title("Hydrous phases -- and Amphibole .- content [vol.%]")
        ax5.set_title("lawsonite content [vol.%]")
        ax6.set_title("Baserock free fluid [vol.%]")

        # ax2, ax3, ax4 die not show a x axis label
        ax1.set_xticklabels([])
        ax2.set_xticklabels([])
        ax3.set_xticklabels([])
        ax4.set_xticklabels([])
        ax5.set_xticklabels([])

        # define the range of the ax2 to ax5 axis and assign the xticks from 300 to 700 with an increment of 50.
        ax1.set_xlim([300, 700])
        ax1.set_xticks(np.arange(300, 700, 50))
        ax2.set_xlim([300, 700])
        ax2.set_xticks(np.arange(300, 700, 50))
        ax3.set_xlim([300, 700])
        ax3.set_xticks(np.arange(300, 700, 50))
        ax4.set_xlim([300, 700])
        ax4.set_xticks(np.arange(300, 700, 50))
        ax5.set_xlim([300, 700])
        ax5.set_xticks(np.arange(300, 700, 50))
        ax6.set_xlim([300, 700])
        ax6.set_xticks(np.arange(300, 700, 50))

        # add a legend to the figure on the right side
        import matplotlib.patches as mpatches
        # create a list of patches for each rock
        patches = []
        for i, rock in enumerate(self.rockdic.keys()):
            patches.append(mpatches.Patch(color=color_palette[i], label=rock))
        # add the patches to the legend and set the legend position
        plt.legend(handles=patches, bbox_to_anchor=(1.0, 1.4), title="Rocks", fontsize=14)

        # adjust the subplot spacing
        plt.subplots_adjust(hspace=0.5)
        # and size to fit the legend
        plt.subplots_adjust(right=0.8)
        plt.tight_layout()

        # save the plot to the directory if img_save is True
        if img_save is True:
            os.makedirs(
                    f'{self.mainfolder}/img_{self.filename}/fluid_content', exist_ok=True)
            plt.savefig(f'{self.mainfolder}/img_{self.filename}/fluid_content/fluid_content.png',
                            transparent=True)
        else:
            plt.show()
        plt.clf()
        plt.close()

    # function to plot a column consisting of all rocks showing the fluid content
    def fluid_content_msg2023(self, img_save=False):
        """Plotting function for the fluid content of all rocks in the model and the glaucophane content in two subplots.

        Args:
            img_save (bool, optional): Optional argument to save the plot as an image to the directory. Defaults to False.
        """
        # read the temperature data from the first rock in the dictionary
        ts = self.rockdic['rock000'].temperature
        # set the color palette for the plot based on the number of rocks
        color_palette = sns.color_palette("viridis", len(self.rockdic.keys()))

        # plotting routine
        fig, (ax1, ax2, ax3, ax6) = plt.subplots(4, 1, figsize=(10, 7))
        for i, rock in enumerate(self.rockdic.keys()):
            # skip the first rock when i == 0
            if i == 0:
                pass
            else:
                # read the database from the rockdic and store it in a variable
                database = self.rockdic[rock].database
                phases2 = list(self.rockdic[rock].phase_data['df_vol%'].columns)
                legend_phases, color_set = phases_and_colors_XMT(database, phases2)

                # read the vol% data to a variable and replace the NaN values with 0 and transpose the dataframe
                y = self.rockdic[rock].phase_data['df_vol%'].fillna(value=0)
                y.columns = legend_phases
                y = y.T

                # clean the dataframe from multiple phase names - combine rows into one
                if len(legend_phases) == len(np.unique(legend_phases)):
                    pass
                else:
                    y, legend_phases, color_set = clean_frame(y, legend_phases, color_set)

                # drop rows with all zeros
                y = y.loc[(y != 0).any(axis=1)]
                legend_phases = list(y.index)

                # define the fluid content from the dataframe based on 'Water' and plot it
                if 'Water' in legend_phases:
                    fluid_content = y.loc['Water']
                    ax2.plot(ts, fluid_content,
                            '-', c=color_palette[i], linewidth=2, markeredgecolor='black')
                    # y axes for ax2 subplot with 10 ticks
                    start, end = ax2.get_ylim()
                    ax2.yaxis.set_ticks(np.round(np.linspace(start, end+end*0.02, 5),1))
                    # plot horizontal line at y = 0
                    ax2.axhline(y=0, color='black', linestyle='-', linewidth=1.0, alpha=0.6)

        # plot the water content of rock000 in ax4
        if 'fluid' in self.rockdic['rock000'].phase_data['df_vol%'].columns:
            water_fluid_content = self.rockdic['rock000'].phase_data['df_vol%']['fluid']
            ax6.plot(ts, water_fluid_content,
                    '-', c='black', linewidth=5, markeredgecolor='black')
            # y axes for ax6 subplot with 10 ticks
            start, end = ax6.get_ylim()
            ax6.yaxis.set_ticks(np.round(np.linspace(start, end+end*0.02, 5),1))
            # plot horizontal line at y = 0
            ax6.axhline(y=0, color='black', linestyle='-', linewidth=1.0, alpha=0.6)

        # plot the water that is release from the least rock in rockdic.keys() in ax1
        top_rock_name = list(self.rockdic.keys())[-1]
        if 'fluid' in self.rockdic[top_rock_name].phase_data['df_vol%'].columns:
            water_release = self.rockdic[top_rock_name].extracted_fluid_volume/1000000 # in m3
            # water release as cumulative sum
            water_release = np.cumsum(water_release)
            ax1.plot(ts, water_release,
                    '-', c='black', linewidth=5, markeredgecolor='black')


        # Set the x-label at lower most subplot represetning all plots
        ax1.set_xlabel("Stack top - water release", fontsize=16)
        ax2.set_xlabel("Aq. fluid distribution", fontsize=16)
        ax6.set_xlabel("Stack base", fontsize=16)
        # set y label
        ax1.set_ylabel("Volume [$m^{3}$]", fontsize=16)
        ax2.set_ylabel("Content [vol.%]", fontsize=16)
        ax6.set_ylabel("Content [vol.%]", fontsize=16)


        """
        # Title each subplot
        ax1.set_title("Water release at top [$m^{3}$]", fontsize=16)
        ax2.set_title("Fluid content [vol.%]", fontsize=16)
        ax6.set_title("Baserock free fluid [vol.%]", fontsize=16)
        """

        # ax2, ax3, ax4 die not show a x axis label
        ax1.set_xticklabels([])
        # ax2.set_xticklabels([])

        # define the range of the ax2 to ax5 axis and assign the xticks from 300 to 700 with an increment of 50.
        ax1.set_xlim([300, 700])
        ax1.set_xticks(np.arange(300, 700, 50))
        ax2.set_xlim([300, 700])
        ax2.set_xticks(np.arange(300, 700, 50))
        ax3.set_xlim([300, 700])
        ax3.set_xticks(np.arange(300, 700, 50))
        ax6.set_xlim([300, 700])
        ax6.set_xticks(np.arange(300, 700, 50))

        # define y axis limits
        ax1.set_ylim([0, np.round(np.max(water_release) + 0.1*np.max(water_release),5)])
        ax2.set_yticks(np.arange(0, np.round(np.max(water_release)+0.1*np.max(water_release)), 5))
        ax2.set_ylim([0, 30])
        ax2.set_yticks(np.arange(0, 30, 5))
        ax6.set_ylim([0, np.round(np.max(water_fluid_content)+0.1*np.max(water_fluid_content),1)])
        ax6.set_yticks(np.arange(0, np.round(np.max(water_fluid_content)+0.1*np.max(water_fluid_content),1), 0.5))
        # tick labels
        ax1.tick_params(axis='y', labelsize=14)
        ax2.tick_params(axis='y', labelsize=14)
        ax6.tick_params(axis='y', labelsize=14)
        ax6.tick_params(axis='x', labelsize=14)
        ax2.tick_params(axis='x', labelsize=14)
        ax3.tick_params(axis='x', labelsize=14)

        # add a legend to the figure on the right side
        import matplotlib.patches as mpatches
        # create a list of patches for each rock
        patches = []
        for i, rock in enumerate(self.rockdic.keys()):
            patches.append(mpatches.Patch(color=color_palette[i], label=rock))
        # add the patches to the legend and set the legend position
        plt.legend(handles=patches, bbox_to_anchor=(1.0, 1.4), title="Rocks", fontsize=14)

        # adjust the subplot spacing
        plt.subplots_adjust(hspace=0.5)
        # and size to fit the legend
        plt.subplots_adjust(right=0.8)
        plt.tight_layout()

        # save the plot to the directory if img_save is True
        if img_save is True:
            os.makedirs(f'{self.mainfolder}/img_{self.filename}/fluid_content', exist_ok=True)
            plt.savefig(f'{self.mainfolder}/img_{self.filename}/fluid_content/fluid_content_msg23.png',
                            transparent=True)
        else:
            plt.show()
        plt.clf()
        plt.close()

    def lawsonite_surface(self, img_save=False):
        """Plotting function for the fluid content of all rocks in the model and the glaucophane content in two subplots.

        Args:
            img_save (bool, optional): Optional argument to save the plot as an image to the directory. Defaults to False.
        """
        # read the temperature data from the first rock in the dictionary
        ts = self.rockdic['rock000'].temperature
        ps = self.rockdic['rock000'].pressure

        # Create a dataframe with columns for the temperature, pressure and columns for the length of self.rockdic
        frame = np.zeros((len(ts), len(self.rockdic)+2))
        data = pd.DataFrame(frame)

        # assigning temperatures and pressures to first two columns
        data.iloc[:,0] = ts
        data.iloc[:,1] = ps

        # looping over each rock to create the surface plot for lawsonite abundance
        for i, rock in enumerate(self.rockdic.keys()):
            # read the database from the rockdic and store it in a variable
            database = self.rockdic[rock].database
            phases2 = list(self.rockdic[rock].phase_data['df_vol%'].columns)
            legend_phases, color_set = phases_and_colors_XMT(database, phases2)

            # read the vol% data to a variable and replace the NaN values with 0 and transpose the dataframe
            y = self.rockdic[rock].phase_data['df_vol%'].fillna(value=0)
            y.columns = legend_phases
            y = y.T

            # clean the dataframe from multiple phase names - combine rows into one
            if len(legend_phases) == len(np.unique(legend_phases)):
                pass
            else:
                y, legend_phases, color_set = clean_frame(y, legend_phases, color_set)
            # drop rows with all zeros
            y = y.loc[(y != 0).any(axis=1)]
            legend_phases = list(y.index)

            # if -Lawsonite- is in y plot it as a dashed line
            if 'Lawsonite' in legend_phases:
                z = y.loc['Lawsonite']
                data.iloc[:,i+2] = z

        # write data to a csv file
        data.to_csv(f'{self.mainfolder}/lawsonite_surface.csv', index=False)

    def fluid_distribution_sgm23(self, img_save=False, gif_save=False, x_axis_log=False):
        """Plotting function for the fluid content of all rocks in the model and the glaucophane content in two subplots.

        Args:
            img_save (bool, optional): Optional argument to save the plot as an image to the directory. Defaults to False.
        """
        # read the temperature data from the first rock in the dictionary
        ts = self.rockdic['rock000'].temperature
        # set the color palette for the plot based on the number of rocks
        color_palette = sns.color_palette("viridis", len(self.rockdic.keys()))

        # z position dependent on the number of rocks in the model (rockdic)
        column_size = len(self.rockdic.keys())*0.1
        z = np.linspace(0, column_size, len(self.rockdic.keys())) + 0.05

        # record the phase assemblage for each rock in the model and store it in a list
        phase_assemblages = []
        # store the formatted data for each rock in the model in a list
        phase_data = []
        # store the amount of amphibole
        amphibole_content = []
        # store the amount of glaucophane
        glaucophane_content = []
        # store the amount of omphacite
        omphacite_content = []
        # store the amount of hydrous phases
        hydrous_content = []
        # store the amount of lawsonite
        lawsonite_content = []
        # store the amount of garnet
        garnet_content = []
        # store the amount of water
        water_content = []
        # store the bulk d18O
        bulk_d18O = []

        # loop over each rock in the model receiving the data for the rock and mineral phases to reformat into matrices for the plot
        for i, rock in enumerate(self.rockdic.keys()):
            # read the database from the rockdic and store it in a variable
            database = self.rockdic[rock].database
            phases2 = list(self.rockdic[rock].phase_data['df_vol%'].columns)
            legend_phases, color_set = phases_and_colors_XMT(database, phases2)

            # read the vol% data to a variable and replace the NaN values with 0 and transpose the dataframe
            y = self.rockdic[rock].phase_data['df_vol%'].fillna(value=0)
            y.columns = legend_phases
            y = y.T

            # clean the dataframe from multiple phase names - combine rows into one
            if len(legend_phases) == len(np.unique(legend_phases)):
                pass
            else:
                y, legend_phases, color_set = clean_frame(y, legend_phases, color_set)

            # drop rows with all zeros
            y = y.loc[(y != 0).any(axis=1)]
            legend_phases = list(y.index)

            # add the legend_phases to the phase_assemblages list
            phase_assemblages.append(legend_phases)
            # add the y dataframe to the phase_data list
            phase_data.append(y)

            # define the fluid content from the dataframe based on 'Water' and plot it
            if 'Water' in legend_phases:
                fluid_val = np.array(y.loc['Water'])
            else:
                fluid_val = np.zeros(len(ts))

            # if -Glaucophane-
            if 'Glaucophane' in legend_phases:
                glaucophane_val = np.array(y.loc['Glaucophane'])
            else:
                glaucophane_val = np.zeros(len(ts))

            # if -Omphacite-
            if 'Omphacite' in legend_phases:
                omphacite_val = np.array(y.loc['Omphacite'])
            else:
                omphacite_val = np.zeros(len(ts))

            # if -Amphibole-
            if 'Amphibole' in legend_phases:
                amphibole_val = np.array(y.loc['Amphibole'])

                # lopp through hydrous phases, test if they are in the legend_phases, get the indices and sum their values from y
                hydrous_phases = ['Muscovite', 'Amphibole', 'Glaucophane', 'Lawsonite', 'Talc', 'Phengite']
                inid_list = []
                for item in hydrous_phases:
                    if item in legend_phases:
                        inid_list.append(legend_phases.index(item))
                hydrous_val = np.array(y.iloc[inid_list].sum(axis=0))
            else:
                amphibole_val = np.zeros(len(ts))

            # if -Lawsonite- 
            if 'Lawsonite' in legend_phases:
                lawsonite_val = np.array(y.loc['Lawsonite'])
            else:
                lawsonite_val = np.zeros(len(ts))

            # get the data for garnet
            if 'Garnet' in legend_phases:
                garnet_val = np.array(y.loc['Garnet'])
                # cumulative sum of garnet content
                garnet_val = np.cumsum(garnet_val)
            else:
                garnet_val = np.zeros(len(ts))

            # get oxygen isotope signature of bulk
            bulk_oxygen_rock = self.rockdic[rock].bulk_deltao_post

            # write the data to the empty list
            amphibole_content.append(amphibole_val)
            glaucophane_content.append(glaucophane_val)
            omphacite_content.append(omphacite_val)
            hydrous_content.append(hydrous_val)
            lawsonite_content.append(lawsonite_val)
            garnet_content.append(garnet_val)
            water_content.append(fluid_val)
            bulk_d18O.append(bulk_oxygen_rock)

        # transforming data to matrix and transpose
        amphibole_content = np.array(amphibole_content).T
        glaucophane_content = np.array(glaucophane_content).T
        omphacite_content = np.array(omphacite_content).T
        hydrous_content = np.array(hydrous_content).T
        lawsonite_content = np.array(lawsonite_content).T
        garnet_content = np.array(garnet_content).T
        water_content = np.array(water_content).T
        bulk_d18O = np.array(bulk_d18O).T

        # creating a bar plot for the glaucophane, omphacite and garnet content at step 22 of the arrays with the rock slice on x axis
        m = glaucophane_content + omphacite_content + garnet_content + lawsonite_content
        m_max = np.max(m)
        # array from 1 to len(self.rockdic.keys()) for the x axis
        for i in range(len(ts)):
            x = np.arange(1, len(self.rockdic.keys())+1)
            # create the figure and add the bar plots for galucophane, omphacite and garnet
            legend_phases, color_set = phases_and_colors_XMT(database, phases2)
            # create data frame of legend phases and color set
            color_phase_frame = pd.DataFrame(color_set, index=legend_phases)
            pt_position = i

            # figure
            fig2, box_ax = plt.subplots(dpi=300)

            # glaucophane section
            if 'Glaucophane' in legend_phases:
                color_position = legend_phases.index('Glaucophane')
                box_ax.bar(x,glaucophane_content[pt_position],
                        bottom=0,
                        edgecolor='black', linewidth=0.5,
                        color=color_set[color_position],
                        label="Glaucophane")

            # omphacite section
            if 'Omphacite' in legend_phases:
                color_position = legend_phases.index('Omphacite')
            elif 'Clinopyroxene' in legend_phases:
                color_position = legend_phases.index('Clinopyroxene')
            else:
                print("Bar plot issue:\n Error, no color assigned for omphacite")
            box_ax.bar(x,omphacite_content[pt_position],
                    bottom=glaucophane_content[pt_position],
                    edgecolor='black', linewidth=0.5,
                    color=color_set[color_position],
                    label="Omphacite")

            # garnet section
            color_position = legend_phases.index('Garnet')
            box_ax.bar(x,garnet_content[pt_position],
                    bottom=omphacite_content[pt_position]+glaucophane_content[pt_position],
                    edgecolor='black', linewidth=0.5,
                    color=color_set[color_position],
                    label="Garnet")

            # plot the lawsonite content
            """color_position = legend_phases.index('Lawsonite')
            box_ax.bar(x,lawsonite_content[pt_position],
                    bottom=omphacite_content[pt_position]+glaucophane_content[pt_position]+garnet_content[pt_position],
                    edgecolor='black', linewidth=0.5,
                    color=color_set[color_position],
                    label="Lawsonite")"""

            box_ax.set_aspect(aspect=0.2)
            # define x label
            box_ax.set_xlabel("Rock stack (bottom-top in cm)", fontsize=12)
            # define y label
            box_ax.set_title(f"Temperature {int(np.round(ts[i],0))} °C\nMineral content [vol.%]", fontsize=12)
            # set tick params
            box_ax.tick_params(axis='both', labelsize=12)
            # set y ticks
            box_ax.set_yticks(np.arange(0, 110, 10))
            # set legend
            handles, labels = box_ax.get_legend_handles_labels()
            box_ax.legend(handles[::-1], labels[::-1],
                           title='Phases', fontsize=8, ncols=3, loc="lower center", bbox_to_anchor=(0.5, -0.77), framealpha=0.0)
            if img_save is True:
                os.makedirs(f'{self.mainfolder}/img_{self.filename}/msg_redistribution_bar', exist_ok=True)
                plt.savefig(f'{self.mainfolder}/img_{self.filename}/msg_redistribution_bar/distribution_{i}.png',
                                transparent=False)
            plt.clf()
            plt.close()
        # condition to save a gif file from the plots
        if gif_save is True:
            frames = []
            for i in range(len(ts)):
                if i == 0:
                    pass
                else:
                    image = imread(f'{self.mainfolder}/img_{self.filename}/msg_redistribution_bar/distribution_{i}.png')
                    frames.append(image)
            imageio.mimsave(f'{self.mainfolder}/img_{self.filename}/msg_redistribution_bar/output.gif', frames, duration=400)


        # save the fluid data to a txt file
        fluid_to_txt_to_julia = True
        if fluid_to_txt_to_julia is True:
            # save water content to matrix txt file
            np.savetxt(f"{self.mainfolder}/fluid_matrix.txt", water_content)
            # save garnet content to matrix txt file
            np.savetxt(f"{self.mainfolder}/garnet_matrix.txt", garnet_content)
            # save glaucophane content to matrix txt file
            np.savetxt(f"{self.mainfolder}/glaucophane_matrix.txt", glaucophane_content)
            # save omphacity content to matrix txt file
            np.savetxt(f"{self.mainfolder}/omphacite_matrix.txt", omphacite_content)
            # save d18O content to matrix txt file
            np.savetxt(f"{self.mainfolder}/d18O_matrix.txt", bulk_d18O)
            # run julia script to plot the fluid distribution - not working yet
            """
            import julia
            j = julia.Julia()
            xxx = j.include(f"julia {self.mainfolder}/fluid_surface.jl")
            """
            # os.system(f"julia {self.mainfolder}/fluid_surface.jl")


        # plotting routine for evolving system
        # looping over length of temperature array to plot each model increment
        for i in range(len(ts)):
            # root figure
            fig, (ax1, ax2, ax3, ax4, ax5, ax6, ax7) = plt.subplots(1, 7, figsize=(11, 9), sharey=True)
            # Title each subplot
            ax1.set_title("Omphacite")
            ax2.set_title("Glaucophane")
            ax3.set_title("Hydrouse ph.")
            ax4.set_title("Lawsonite")
            ax5.set_title("garnet")
            ax6.set_title("Water")
            ax7.set_title("Bulk $\delta^{18}$O \n[‰ vs. VSMOW]")

            # set y axis label for the first subplot
            ax1.set_ylabel("Thickness [m] (rock increment = 0.1 m)")

            # define x lims for each subplot
            ax1.set_xlim([0, omphacite_content.max()+0.1*omphacite_content.max()])
            ax2.set_xlim([0, glaucophane_content.max()+0.1*glaucophane_content.max()])
            ax3.set_xlim([0, hydrous_content.max()+0.1*hydrous_content.max()])
            ax4.set_xlim([0, lawsonite_content.max()+0.1*lawsonite_content.max()])
            ax5.set_xlim([0, garnet_content.max()+0.1*garnet_content.max()])
            ax6.set_xlim([0, water_content.max()+0.1*water_content.max()])
            ax7.set_xlim([0, 25])

            # activate log scaled x axis if argument x_axis_log is True
            if x_axis_log is True:
                ax1.set_xscale('log')
                ax3.set_xscale('log')
                ax5.set_xscale('log')

            # define y lims for each subplot
            ax1.set_ylim([0, column_size])
            ax2.set_ylim([0, column_size])
            ax3.set_ylim([0, column_size])
            ax4.set_ylim([0, column_size])
            ax5.set_ylim([0, column_size])
            ax6.set_ylim([0, column_size])
            ax7.set_ylim([0, column_size])

            # plotting
            ax1.plot(omphacite_content[i], z, '-', c='black', linewidth=1.2, markeredgecolor='black')
            ax2.plot(glaucophane_content[i], z, '-', c='black', linewidth=1.2, markeredgecolor='black')
            ax3.plot(hydrous_content[i], z, '-', c='black', linewidth=1.2, markeredgecolor='black')
            ax4.plot(lawsonite_content[i], z, '-', c='black', linewidth=1.2, markeredgecolor='black')
            ax5.plot(garnet_content[i], z, '-', c='black', linewidth=1.2, markeredgecolor='black')
            ax6.plot(water_content[i], z, '-', c='black', linewidth=1.2, markeredgecolor='black')
            ax7.plot(bulk_d18O[i], z, '-', c='black', linewidth=1.2, markeredgecolor='black')

            # Figure title writing the temperature step
            fig.suptitle(f"Model phases [vol.%] at {np.round(ts[i], 1)}°C", fontsize=16)

            """# adjust the subplot spacing
            plt.subplots_adjust(hspace=0.5)
            # and size to fit the legend
            plt.subplots_adjust(right=0.8)
            plt.tight_layout()"""

            # save the plot to the directory if img_save is True
            if img_save is True:
                os.makedirs(
                        f'{self.mainfolder}/img_{self.filename}/msg_redistribution', exist_ok=True)
                plt.savefig(f'{self.mainfolder}/img_{self.filename}/msg_redistribution/redistribution_{i}.png',
                                transparent=False, facecolor='white')
            else:
                plt.show()
            plt.clf()
            plt.close()
            print(f"plot {i} done")

        # condition to save a gif file from the plots
        if gif_save is True:
            frames = []
            for i in range(len(ts)):
                if i == 0:
                    pass
                else:
                    image = imread(f'{self.mainfolder}/img_{self.filename}/msg_redistribution/redistribution_{i}.png')
                    frames.append(image)
            imageio.mimsave(f'{self.mainfolder}/img_{self.filename}/msg_redistribution/output.gif', frames, duration=400)



if __name__ == '__main__':

    # Read variables from a hdf5 output file from ThorPT
    data = ThorPT_hdf5_reader()
    data.open_ThorPT_hdf5()

    # Activate the plotting module - Predefined plots for evaluation
    compPlot = ThorPT_plots(
        data.filename, data.mainfolder, data.rock, data.compiledrock)

    """
    track_diff = []
    track_shear = []
    track_num_ext = []
    for i, item in enumerate(data.compiledrock.all_diffs):
        diff = np.unique(item)[-1]
        fracture_bool = data.compiledrock.all_frac_bool[i]
        num_ext = len(fracture_bool[fracture_bool>0])
        track_diff.append(diff)
        # track_shear.append()
        track_num_ext.append(num_ext)


    for i, item in enumerate(data.rock.keys()):
        # y = data.compiledrock.all_depth[i]/1000
        y = data.compiledrock.all_ps[i]
        x = data.compiledrock.all_porosity[i]*100
        x[x==0] = np.nan

        plt.plot(x,y,'d')
        plt.ylabel("Depth [km]")
        plt.xlabel("Fluid-filled porosity [Vol.%]")
        plt.ylim(100, 0)

        num_rocks = 3
        step = int(len(data.compiledrock.all_porosity)/num_rocks)
        for i in range(num_rocks):
            x = np.concatenate(data.compiledrock.all_porosity[step*i:step*(i+1)])*100
            y = np.concatenate(data.compiledrock.all_ps[step*i:step*(i+1)])
            plt.plot(x,y, 'd')
        plt.xscale('log')
        plt.legend(['Basalt', 'Sediment', 'Serpentinite'])

        from matplotlib import colors
        i = 0
        dist1 = np.concatenate(data.compiledrock.all_porosity[step*i:step*(i+1)])*100
        dist2 = np.concatenate(data.compiledrock.all_ps[step*i:step*(i+1)])
        fig, ax = plt.subplots(tight_layout=True)
        hist = ax.hist2d(dist1, dist2, bins=128, norm=colors.LogNorm())

    for i, item in enumerate(data.rock.keys()):
        y = data.compiledrock.all_depth[i]/1000
        fracture_bool = data.compiledrock.all_frac_bool[i]
        num_ext = len(fracture_bool[fracture_bool>0])
        # plt.plot(num_ext,y,'d')"""


    """
    plt.plot(track_diff, track_num_ext, 'd-')
    plt.ylabel("# of extractions")
    plt.xlabel("Differential stress [$\it{MPa}$]")
    """

    # Protocol for lawsonite surface plot
    lawsonite_surface = False
    if lawsonite_surface is True:
        # Plotting tool to investigate the lawsonite occurence of a grid
        compPlot.lawsonite_surface(True)

    # Protocol for standard plotting (Stack, oxygen isotopes, release fluid volume)
    standard_plots = True
    if standard_plots is True:
            for key in data.rock.keys():
                print(key)
                compPlot.phases_stack_plot(rock_tag=key, img_save=True,
                            val_tag='volume', transparent=True, fluid_porosity=True)
                compPlot.oxygen_isotopes(rock_tag=key, img_save=True)
                # compPlot.release_fluid_volume_plot(rock_tag=key, img_save=False)

            # compPlot.pt_path_plot(key, img_save=True, gif_save=True)

    # modelling for msg23 - internal redistribution of fluids in an outcrop
    fluid_distribution_sgm23 = True
    if fluid_distribution_sgm23 is True:
        compPlot.fluid_distribution_sgm23(img_save=True, gif_save=True, x_axis_log=False)

    fluid_recycling_single = False
    if fluid_recycling_single is True:
        # Routine for Fluid recycling plotting
        print("The 'data.rock[key]' keys are:")
        compPlot.fluid_content(img_save=False)

    plotting_routine_for_sgm23 = False
    if plotting_routine_for_sgm23 is True:
        # Routine for Fluid recycling plotting
        compPlot.fluid_content_msg2023(img_save=True)
        # compPlot.fluid_distribution_sgm23(img_save=True, gif_save=True, x_axis_log=False)



    # routine to plot the fluid content of all rocks and internal recycling
    fluid_recycling = False
    if fluid_recycling is True:
        # Routine for Fluid recycling plotting
        compPlot.time_int_flux_summary(img_save=True, rock_colors=['#10e6ad', "#0592c1", "#f74d53"])

        file_combination = False
        num_files = 1
        if file_combination is True:
            data2 = ThorPT_hdf5_reader()
            data_box = []
            for i in range(num_files):
                data2.open_ThorPT_hdf5()
                data_box.append(copy.deepcopy(data2))

            depth_porosity(
                data_box, num_rocks=3,
                rock_colors=["#0592c1", "#f74d53", '#10e6ad'],
                rock_symbol = ['d', 'd', 'd', 's'])

            density_porosity(
                data_box, num_rocks=3,
                rock_colors=["#0592c1", "#f74d53", '#10e6ad'],
                rock_symbol = ['d', 'd', 'd','s'])

            extraction_stress(
                data_box, num_rocks=3,
                rock_colors=["#0592c1", "#f74d53", '#10e6ad'],
                rock_symbol = ['d--', 'd--', 'd--', 's--'], rock_names=True)



        ##########################################
        print("The 'data.rock[key]' keys are:")
        compPlot.fluid_content(img_save=True)
        evaluate_loop = False
        if evaluate_loop is True:
            for key in data.rock.keys():
                print(key)
                compPlot.phases_stack_plot(rock_tag=key, img_save=True,
                            val_tag='volume', transparent=True, fluid_porosity=True)

                #compPlot.boxplot_to_GIF(rock_tag='rock004', img_save=True, gif_save=True)
                compPlot.time_int_flux_plot(rock_tag=key, img_save=True)

                compPlot.oxygen_isotopes(rock_tag=key, img_save=True)
                # compPlot.release_fluid_volume_plot(rock_tag=key, img_save=False)

            # compPlot.pt_path_plot(key, img_save=True, gif_save=True)


        """
        compPlot.time_int_flux_plot(
                rock_tag='rock1', img_save=False, gif_save=False)
        compPlot.boxplot_to_GIF(rock_tag='rock1', img_save=True, gif_save=True)

        compPlot.binary_plot(rock_tag='rock0')
        compPlot.pt_path_plot(rock_tag='rock0', img_save=False, gif_save=False)
        compPlot.permeability_plot(
            rock_tag='rock0', img_save=False, gif_save=False)
        compPlot.time_int_flux_plot(
            rock_tag='rock0', img_save=False, gif_save=False)
        compPlot.porosity_plot(rock_tag='rock0', img_save=False, gif_save=False)
        compPlot.release_fluid_volume_plot(
            rock_tag='rock0', img_save=False, gif_save=False)"""

        """dist1 = np.concatenate(data_box[0].compiledrock.all_porosity)*100
        dist2 = np.concatenate(data_box[1].compiledrock.all_porosity)*100
        plt.figure(102)
        plt.hist(dist2, bins=128)
        plt.yscale('log')
        plt.figure(103)
        plt.hist(dist1, bins=128)
        plt.yscale('log')
        """
        """from matplotlib import colors
        dist1 = np.concatenate(data_box[0].compiledrock.all_porosity)*100
        dist2 = np.concatenate(data_box[1].compiledrock.all_porosity)*100
        fig, ax = plt.subplots(tight_layout=True)
        hist = ax.hist2d(dist1, dist2, bins=128, norm=colors.LogNorm())
        fig, (ax2,ax3) = plt.subplots((1,2))
        ax2.hist(dist1, bins=128)
        ax3.hist(dist2, bins=128)
        plt.show()"""

        """from matplotlib import colors
        dist1 = np.concatenate(data_box.compiledrock.all_porosity)*100
        dist2 = np.concatenate(data_box.compiledrock.all_system_density_post)-np.concatenate(data.compiledrock.all_system_density_pre)
        fig, ax = plt.subplots(tight_layout=True)
        hist = ax.hist2d(dist1, dist2, bins=128, norm=colors.LogNorm())

        frac = np.concatenate(data.compiledrock.all_frac_bool)
        x = np.concatenate(data.compiledrock.all_porosity)*100
        y1 = np.concatenate(data.compiledrock.all_system_density_post)
        y2 = np.concatenate(data.compiledrock.all_system_density_pre)
        y = y1-y2
        y = y[frac>0]
        x = x[frac>0]
        plt.plot(x,y, 'd')"""

    # Multiple file option for comparing
    # Read variables from a hdf5 output file from ThorPT
    multi_file_routine = True
    # Starting loop over multiple files
    if multi_file_routine is True:
        num_files = 1
        # reading each file and storing it in a list
        data_box = []
        for i in range(num_files):
            # Function call, open hdf5 and append to data_box
            data2 = ThorPT_hdf5_reader()
            data2.open_ThorPT_hdf5()
            data_box.append(copy.deepcopy(data2))

        print("The 'data.rock[key]' keys are:")

        # read all porosity values for each rock in the model
        all_porosity = np.array(data_box[0].compiledrock.all_porosity)

        color_palette = sns.color_palette("viridis", len(all_porosity.T))

        plt.figure(101, figsize=(6,6))
        # figure with aspect ratio of 1:1
        plt.gca().set_aspect('equal', adjustable='box')
        plt.rc('axes', labelsize=16)
        plt.rc('ytick', labelsize=12)  # fontsize of the y tick labels
        plt.rc('xtick', labelsize=12)  # fontsize of the x tick labels
        # plot the porosity of each rock in the model vs the first rock in the model
        for i, item in enumerate(data_box[0].rock.keys()):
            plt.scatter(all_porosity[i], all_porosity[0], color=color_palette[i])
        plt.plot(all_porosity[0], all_porosity[0], color="black", linestyle='--', linewidth=1.0, alpha=0.6)
        plt.xlabel("Fluid-filled porosity [Vol.%]")
        plt.ylabel("Fluid filled porosity\n continous release")
        plt.show()

        plt.figure(102)
        plt.rc('axes', labelsize=16)
        plt.rc('ytick', labelsize=12)  # fontsize of the y tick labels
        plt.rc('xtick', labelsize=12)  # fontsize of the x tick labels
        # normalized fluid-filled porosity prograde plot (extraction steps/always extraction = 1:1 if same)
        porosity = all_porosity/all_porosity[0]
        for i, item in enumerate(porosity):
            plt.plot(porosity[i], c=color_palette[i], marker='o', linestyle='--', linewidth=1.0, alpha=0.6)
        plt.xlabel("Prograde met. aka. time")
        plt.ylabel("Fluid filled porosity norm.")
        plt.show()

        applied_diff_stress = np.array(data_box[0].compiledrock.all_diffs).T[-1]
        applied_diff_stress = np.around(applied_diff_stress,1)
        used_tensile_strengths = data_box[0].compiledrock.all_tensile

        # do the colorbar in a seperate plot
        plt.figure(103)
        plt.gca().set_aspect(aspect=0.2, adjustable='box')
        # plotting a bar for each rock along the x axis in the model using color_palette
        plt.barh(np.arange(len(data_box[0].rock.keys())), 1, color=color_palette)
        # annotate the bar plot with the applied differential stress and tensile strength
        for i, item in enumerate(data_box[0].rock.keys()):
            plt.annotate(f"{applied_diff_stress[i]} MPa, {used_tensile_strengths[i]} MPa",
                         (1, i), color='black', fontsize=5)
        plt.xlim([0, 2])
        plt.title("Colorbar for rocks")
        # remove axes visibility
        plt.gca().axes.get_yaxis().set_visible(False)
        plt.gca().axes.get_xaxis().set_visible(False)
        # remove the frame
        plt.gca().set_frame_on(False)
        plt.show()

        # number of extraction vs differential stress
        track_diff = []
        track_shear = []
        track_num_ext = []
        for i, item in enumerate(data.compiledrock.all_diffs):
            diff = np.unique(item)[-1]
            fracture_bool = data.compiledrock.all_frac_bool[i]
            num_ext = len(fracture_bool[fracture_bool>0])
            track_diff.append(diff)
            # track_shear.append()
            track_num_ext.append(num_ext)
        plt.plot(track_diff, track_num_ext, 'd-')
        plt.ylabel("# of extractions")
        plt.xlabel("Differential stress [$\it{MPa}$]")

        # get colour for each value in track_diff depending on the float value
        color_palette = sns.color_palette("viridis", int(max(track_diff)))
        # create a list of colours for each value in track_diff
        color_list = []
        for i, item in enumerate(track_diff):
            color_list.append(color_palette[int(item)-1])

        # import mpatches
        import matplotlib.patches as mpatches
        #plot the number of extractions vs the used tensile strength with the color from color_list
        plt.figure(105)
        for i, item in enumerate(track_num_ext):
            plt.plot(used_tensile_strengths[i], track_num_ext[i], 'd', color=color_list[i])
        plt.xlabel("Tensile strength [MPa]")
        plt.ylabel("# of extractions")
        # get the legend of unique values in track_diff with the color from color_palette
        patches = []
        for i, item in enumerate(np.unique(track_diff)):
            patches.append(mpatches.Patch(color=color_palette[int(item)-1], label=item))
        plt.legend(handles=patches, bbox_to_anchor=(1.0, 1.0), title="Differential stress [MPa]", fontsize=12)

        plt.show()



