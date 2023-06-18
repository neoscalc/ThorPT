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
import imageio
from imageio import imread
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
        f"MINERAL_NAMES_{database[:-dot_indi]}_to_XMT.txt"
    min_translation = pd.read_csv(file_to_open, delimiter='\t')

    # Iterating through phases and select colors
    colpal = sns.color_palette("flare", 20)
    color_set = []
    phase_set = []
    z = 0
    for mini in phases:
        if mini == 'water.fluid' or mini == 'H2O.liq':
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

    if 'LIQtc_h2oL' in frame.index and 'water.fluid' in frame.index:
        water_1 = frame.loc['LIQtc_h2oL']
        water_2 = frame.loc['water.fluid']
        water_2[water_2.isna()] = water_1[water_2.isna()]
        frame.loc['water.fluid'] = water_2
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

    if 'LIQtc_h2oL' in frame.index and 'water.fluid' in frame.index:
        water_1 = frame.loc['LIQtc_h2oL']
        water_2 = frame.loc['water.fluid']
        water_2[water_2.isna()] = water_1[water_2.isna()]
        frame.loc['water.fluid'] = water_2
        frame = frame.drop('LIQtc_h2oL', axis=0)

    return frame


def create_gif(phase_data, mainfolder, filename, group_key, subfolder='default'):
    # reading images and put it to GIF
    frames = []
    for i in range(len(phase_data.index)):
        image = imread(
            f'{mainfolder}/img_{filename}/{subfolder}/{group_key}/img_{i}.png')
        frames.append(image)
    imageio.mimsave(f'{mainfolder}/img_{filename}/{subfolder}/{group_key}/output.gif', frames, duration=1)


# Progressbar init
def progress(percent=0, width=40):
    left = width * percent // 100
    right = width - left
    tags = "#" * left
    spaces = " " * right
    percents = f"{percent:.0f}%"
    print("\r[", tags, spaces, "]", percents, sep="", end="", flush=True)


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
class Rock:
    name: str
    # temperature: float = field(metadata={'unit': 'degree C'})
    temperature: any
    pressure: any
    depth: any

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

@dataclass
class CompiledData:
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
    all_fluid_before: any
    all_phases: any
    all_phases2: any
    all_database: any


class ThorPT_hdf5_reader():

    def __init__(self):
        self.rock = {}
        self.rock_names = []
        self.compiledrock = 0
        self.mainfolder = False
        self.filename = False

    def open_ThorPT_hdf5(self):

        # lists for compiling data
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
        all_fluid_before = []
        all_phases = []
        all_phases2 = []
        all_database = []
        all_td_data = []

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
                ts = np.array(f[group_key]['temperatures'])
                ps = np.array(f[group_key]['pressures'])
                vtime = np.array(f[group_key]['time_frame'])
                vtime = vtime.T
                extr_time = np.array(f[group_key]['extr_time'])
                frac_bool = np.array(f[group_key]['fracture bool'])
                depth = np.array(f[group_key]['depth'])
                geometry = list(f[group_key]['geometry'].asstr())

                tensile = np.float64(f[group_key]['tensile strength'])

                if 'line' in f[group_key].keys():
                    line = np.array(f[group_key]['line'])
                else:
                    line = False

                # Modelled params
                permea = np.array(f[group_key]['permeability'])
                t_flux = np.array(f[group_key]['time-int flux'])
                virtual_perm = np.array(f[group_key]['virtual permeability'])
                virtual_flux = np.array(f[group_key]['virtual time-int flux'])

                # Fluid extraction
                extr_d = pd.DataFrame(f[group_key]['extracted_fluid_data'])
                fluid_before = np.array(f[group_key]['st_fluid_before'])

                # get some physics arrays
                sys_physc = {}
                for item in f[group_key]['sys_physicals']:
                    sys_physc[item] = np.array(
                        f[group_key]['sys_physicals'][item])

                system_vol_pre = sys_physc['system_vol_pre']
                system_vol_post = sys_physc['system_vol_post']

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
                d_oxy = pd.DataFrame(f[group_key]['save_oxygen'])
                d_oxy_phases = list(f[group_key].attrs['oxy_data_phases'])
                d_oxy.columns = d_oxy_phases
                save_bulk_oxygen_pre = np.array(
                    f[group_key]['save_bulk_oxygen_pre'])
                save_bulk_oxygen_post = np.array(
                    f[group_key]['save_bulk_oxygen_post'])

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

                volume_data = copy.deepcopy(df_var_d['df_volume[ccm]'])
                volume_data.columns = phases
                volume_data[np.isnan(volume_data) == True] = 0

                mass_abs_data = copy.deepcopy(df_var_d['df_wt[g]'])
                mass_abs_data.columns = phases
                mass_abs_data[np.isnan(mass_abs_data) == True] = 0

                solid_volumes = volume_data.T.sum()-volume_data['water.fluid']
                solid_weight = mass_abs_data.T.sum(
                )-mass_abs_data['water.fluid']
                solid_density = np.array(
                    solid_weight)/np.array(solid_volumes)*1000

                # massperc_fluid = mass_data['water.fluid']
                massperc_fluid = mass_abs_data['water.fluid'] / \
                    (mass_abs_data.T.sum()-mass_abs_data['water.fluid'])
                q_ti = massperc_fluid*solid_density * \
                    (1-unfiltered_porosity)*bloc_a/area

                # the new permeability
                # ANCHOR
                mü_water = 0.0001
                density_cont = solid_density - \
                    mass_abs_data['water.fluid']/volume_data['water.fluid']
                permeability2 = q_ti/(151000*365*24*60*60) * \
                    mü_water / 9.81 / density_cont

                # extraction marker
                # NOTE - new masking method
                extraction_boolean = np.array(frac_bool, dtype=bool)
                extraction_boolean = np.invert(extraction_boolean)

                # Compiling section
                all_tensile.append(tensile)
                all_diffs.append(np.array(f[group_key]['diff. stress']))
                all_frac_bool.append(frac_bool)
                arr_line.append(line)
                all_permea.append(permea)
                all_t_flux.append(t_flux)
                all_depth.append(depth)
                all_virtual_perm.append(virtual_perm)
                all_virtual_flux.append(virtual_flux)
                all_geometry.append(geometry)
                all_system_vol_pre.append(system_vol_pre)
                all_system_vol_post.append(system_vol_post)
                all_fluid_before.append(fluid_before)
                all_phases.append(phases)
                all_phases2.append(phases2)
                all_database.append(database)
                all_td_data.append(td_data)

                # Write the dataclass
                self.rock[group_key] = Rock(
                    name=None,
                    temperature=ts,
                    pressure=ps,
                    depth=depth,
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
                    oxygen_data = d_oxy,
                    bulk_deltao_pre= save_bulk_oxygen_pre,
                    bulk_deltao_post= save_bulk_oxygen_post,
                    group_key=group_key,
                    extraction_boolean=extraction_boolean,
                    extracted_fluid_volume=v_fluid_extr,
                    permeability=permea,
                    porosity=unfiltered_porosity,
                    time_int_flux=unfiltered_int_flux,
                    time_int_flux2=q_ti,
                    v_permeability=0,
                    v_timeintflux=0)

                self.compiledrock = CompiledData(
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
                    all_fluid_before,
                    all_phases,
                    all_phases2,
                    all_database)


class ThorPT_plots():

    def __init__(self, filename, mainfolder, rockdata, compiledrockdata):

        self.filename = filename
        self.mainfolder = mainfolder
        self.rockdic = rockdata
        self.comprock = compiledrockdata

    def boxplot_to_GIF(self, rock_tag, img_save=False, gif_save=False):

        if gif_save is True:
            img_save=True

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
                ax1 = fig.add_subplot(gs[:, :2])
                for t, val in enumerate(np.array(mole_fractions)):
                    if val == np.float64(0):
                        ax1.bar(
                            1, val, bottom=0, color=color_set[t], edgecolor='black', linewidth=0.1, label=legend_phases[t])
                    else:
                        ax1.bar(1, val, bottom=y_offset,
                                color=color_set[t], edgecolor='black', linewidth=0.1, label=legend_phases[t])
                    y_offset = y_offset + val
                # legend
                handles, labels = ax1.get_legend_handles_labels()
                legend_properties = {'weight': 'bold'}
                ax1.legend(handles[::-1], labels[::-1],
                           title='Phases', bbox_to_anchor=(1.2, 0.8), fontsize=18)
                ax1.get_xaxis().set_visible(False)
                ax1.get_yaxis().set_visible(False)
                ax1.set_aspect(0.02)

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
            img_save=True

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
                ax4 = fig.add_subplot(gs[:, :])
                ax4.plot(ts[:i+1], ps[:i+1]/10_000, 'd--',
                         color='black', markersize=4)
                ax4.plot(ts[i:i+1], ps[i:i+1]/10_000, 'd--',
                         color='#7fffd4', markersize=8, markeredgecolor='black')
                # Ax4 figure features
                ax4.set_xlim(200, 750)
                ax4.set_ylim(0, 3)
                ax4.set_ylabel("Pressure [GPa]", fontsize=18)
                ax4.set_xlabel("Temperature [°C]", fontsize=18)
                ax4.tick_params(axis='both', labelsize=18)

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
            ax4 = fig.add_subplot(gs[:, :])
            for i in range(len(phase_data.index)):
                ax4.plot(ts[:i+1], ps[:i+1]/10_000, 'd--',
                         color='black', markersize=4)
                ax4.plot(ts[i:i+1], ps[i:i+1]/10_000, 'd--',
                         color='#7fffd4', markersize=8, markeredgecolor='black')
            # Ax4 figure features
            ax4.set_xlim(200, 750)
            ax4.set_ylim(0, 3)
            ax4.set_ylabel("Pressure [GPa]", fontsize=18)
            ax4.set_xlabel("Temperature [°C]", fontsize=18)
            ax4.tick_params(axis='both', labelsize=18)
            plt.show()

        if gif_save is True:
            # call gif function
            create_gif(phase_data, self.mainfolder, self.filename,
                       group_key, subfolder=subfolder)

    def permeability_plot(self, rock_tag, img_save=False, gif_save=False):

        if gif_save is True:
            img_save=True

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

                ax4 = fig.add_subplot(gs[:, :])
                if False in filtered_permeability.mask[:i+1]:
                    # plot 4
                    ax4.plot(
                        filtered_permeability[:i+1], depth[:i+1]/1000, 'd--', color='black', alpha=0.7)
                    ax4.plot(filtered_permeability[i:i+1], depth[i:i+1]/1000, 'd',
                             color='#7fffd4', markersize=8, markeredgecolor='black')
                else:
                    # plot 4
                    ax4.plot(np.ones(i+1)*1e-30, depth[:i+1]/1000, 'd')

                # add Manning & Ingebritsen to plot 3
                z = np.arange(1, 110, 1)
                c_k = 10**(-14 - 3.2 * np.log10(z))
                ax4.plot(c_k, z, '--', color='#b18e4e', linewidth=3)

                # Ax4 figure features
                ax4.set_title('T:{:.2f} °C P:{:.2f} GPa'.format(
                    ts[i], ps[i]/10000), fontsize=18)
                ax4.set_xlim(1e-26, 1e-14)
                ax4.set_xscale('log')
                ax4.set_ylim(100, 0)
                ax4.set_xlabel("permeability [$m^{2}$]", fontsize=18)
                ax4.set_ylabel("Depth [km]", fontsize=18)
                ax4.fill_between([1e-30, 1e-21], [100, 100],
                                 color="#c5c5c5", edgecolor='black', hatch="/", alpha=0.4)
                ax4.tick_params(axis='both', labelsize=20)
                ax4.annotate("Manning &\nIngebritsen\n(1999)",
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
            ax4 = fig.add_subplot(gs[:, :])
            for i in range(len(phase_data.index)):
                if False in filtered_permeability.mask[:i+1]:
                    # plot 4
                    ax4.plot(
                        filtered_permeability[:i+1], depth[:i+1]/1000, 'd--', color='black', alpha=0.7)
                    ax4.plot(filtered_permeability[i:i+1], depth[i:i+1]/1000, 'd',
                             color='#7fffd4', markersize=8, markeredgecolor='black')
                else:
                    # plot 4
                    ax4.plot(np.ones(i+1)*1e-30, depth[:i+1]/1000, 'd')

            # add Manning & Ingebritsen to plot 3
            z = np.arange(1, 110, 1)
            c_k = 10**(-14 - 3.2 * np.log10(z))
            ax4.plot(c_k, z, '--', color='#b18e4e', linewidth=3)

            # Ax4 figure features
            ax4.set_title('T:{:.2f} °C P:{:.2f} GPa'.format(
                ts[i], ps[i]/10000), fontsize=18)
            ax4.set_xlim(1e-26, 1e-14)
            ax4.set_xscale('log')
            ax4.set_ylim(100, 0)
            ax4.set_xlabel("permeability [$m^{2}$]", fontsize=18)
            ax4.set_ylabel("Depth [km]", fontsize=18)
            ax4.fill_between([1e-30, 1e-21], [100, 100],
                             color="#c5c5c5", edgecolor='black', hatch="/", alpha=0.4)
            ax4.tick_params(axis='both', labelsize=20)
            ax4.annotate("Manning &\nIngebritsen\n(1999)",
                         (4e-17, 15), fontsize=16)
            plt.show()

        # reading images and put it to GIF
        if gif_save is True:
            # call gif function
            create_gif(phase_data, self.mainfolder, self.filename,
                       group_key, subfolder=subfolder)

    def time_int_flux_plot(self, rock_tag, img_save=False, gif_save=False):

        if gif_save is True:
            img_save=True

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

        if len(np.unique(self.rockdic[rock_tag].frac_bool)) == 1 and np.unique(self.rockdic[rock_tag].frac_bool)[-1] == 0:
            filtered_q_ti = np.zeros(len(self.rockdic[rock_tag].frac_bool))
            filtered_int_flux = np.zeros(len(self.rockdic[rock_tag].frac_bool))

        # saving subfolder
        subfolder = 'time_int_flux'

        regional_filtered_flux = filtered_q_ti

        if img_save is True:

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

                ax4 = fig.add_subplot(gs[:, :])
                # ax4.plot(unfiltered_int_flux[:i+1],
                #          depth[:i+1]/1000, 'd--', color='green')
                # ax4.plot(regional_flux[:i+1],
                #          depth[:i+1]/1000, 'd--', color='green', markersize=3, alpha=0.5)
                if False in filtered_int_flux.mask[:i+1]:
                    # plot 4
                    # ax4.plot(filtered_int_flux[:i+1], depth[:i+1]/1000, 'd--', color='black', alpha=0.7)
                    # ax4.plot(filtered_int_flux[i:i+1], depth[i:i+1]/1000, 'd',
                    #          color='#7fffd4', markersize=8, markeredgecolor='black')
                    ax4.plot(
                        regional_filtered_flux[:i+1], depth[:i+1]/1000, 'd--', color='black', alpha=0.7, markersize=10)
                    ax4.plot(regional_filtered_flux[i:i+1], depth[i:i+1]/1000, 'd',
                             color='#7fffd4', markersize=15, markeredgecolor='black')
                else:
                    # plot 4
                    ax4.plot(np.ones(i+1)*1e-30, depth[:i+1]/1000, 'd')

                # Ax4 figure features
                ax4.set_title('T:{:.2f} °C P:{:.2f} GPa'.format(
                    ts[i], ps[i]/10000), fontsize=18)
                ax4.set_xscale('log')
                ax4.set_xlim(10**0, 10**6)
                ax4.set_ylim(100, 0)
                ax4.set_xlabel("Time int.-flux [$m^{3}/m^{2}$]", fontsize=18)
                ax4.set_ylabel("Depth [km]", fontsize=18)

                # pervasive range
                """ax4.fill_between([10**(2), 10**(4)], [100, 100],
                                color="#a9d5b2", edgecolor='black', alpha=0.1)"""
                ax4.fill_between([10**(1), 10**(4)], [100, 100],
                                 cmap="Purples", alpha=0.2, step='mid')
                # channelized range
                ax4.fill_between([10**(4), 10**(6)], [100, 100],
                                 color="#ca638f", alpha=0.1)
                # regional range
                """ax4.fill_between([10**(2.7-0.5), 10**(2.7+0.5)], [100, 100],
                                color="#c5c5c5", edgecolor='black', alpha=0.5)"""

                ax4.tick_params(axis='both', labelsize=20)
                ax4.vlines(10**4, 0, 100, linestyles='--',
                           color='black', linewidth=4)

                props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

                # ax4.annotate("Regional range", (10**(2.7+0.5-0.3), 25), fontsize=16, bbox=props, rotation=90)
                ax4.annotate("Channelized", (10**(4+0.2), 5),
                             fontsize=16, bbox=props)
                ax4.annotate("Pervasive", (10**(1+0.2), 5),
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
            ax4 = fig.add_subplot(gs[:, :])
            for i in range(len(phase_data.index)):
                # ax4.plot(unfiltered_int_flux[:i+1],
                #          depth[:i+1]/1000, 'd--', color='green')
                # ax4.plot(regional_flux[:i+1],
                #          depth[:i+1]/1000, 'd--', color='green', markersize=3, alpha=0.5)
                if False in filtered_int_flux.mask[:i+1]:
                    # plot 4
                    # ax4.plot(filtered_int_flux[:i+1], depth[:i+1]/1000, 'd--', color='black', alpha=0.7)
                    # ax4.plot(filtered_int_flux[i:i+1], depth[i:i+1]/1000, 'd',
                    #          color='#7fffd4', markersize=8, markeredgecolor='black')
                    ax4.plot(
                        regional_filtered_flux[:i+1], depth[:i+1]/1000, 'd--', color='black', alpha=0.7, markersize=10)
                    ax4.plot(regional_filtered_flux[i:i+1], depth[i:i+1]/1000, 'd',
                             color='#7fffd4', markersize=15, markeredgecolor='black')
                else:
                    # plot 4
                    ax4.plot(np.ones(i+1)*1e-30, depth[:i+1]/1000, 'd')

            # Ax4 figure features
            ax4.set_xscale('log')
            ax4.set_xlim(10**0, 10**6)
            ax4.set_ylim(100, 0)
            ax4.set_xlabel("Time int.-flux [$m^{3}/m^{2}$]", fontsize=18)
            ax4.set_ylabel("Depth [km]", fontsize=18)
            # pervasive range
            """ax4.fill_between([10**(2), 10**(4)], [100, 100],
                            color="#a9d5b2", edgecolor='black', alpha=0.1)"""
            ax4.fill_between([10**(1), 10**(4)], [100, 100],
                             cmap="Purples", alpha=0.2, step='mid')
            # channelized range
            ax4.fill_between([10**(4), 10**(6)], [100, 100],
                             color="#ca638f", alpha=0.1)
            # regional range
            """ax4.fill_between([10**(2.7-0.5), 10**(2.7+0.5)], [100, 100],
                            color="#c5c5c5", edgecolor='black', alpha=0.5)"""
            ax4.tick_params(axis='both', labelsize=20)
            ax4.vlines(10**4, 0, 100, linestyles='--',
                       color='black', linewidth=4)
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
            # ax4.annotate("Regional range", (10**(2.7+0.5-0.3), 25), fontsize=16, bbox=props, rotation=90)
            ax4.annotate("Channelized", (10**(4+0.2), 5),
                         fontsize=16, bbox=props)
            ax4.annotate("Pervasive", (10**(1+0.2), 5),
                         fontsize=16, bbox=props)
            ax4.set_title('T:{:.2f} °C P:{:.2f} GPa'.format(
                    ts[i], ps[i]/10000), fontsize=18)
            plt.show()

        # reading images and put it to GIF
        if gif_save is True:
            # call gif function
            create_gif(phase_data, self.mainfolder, self.filename,
                       group_key, subfolder=subfolder)

    def porosity_plot(self, rock_tag, img_save=False, gif_save=False):

        if gif_save is True:
            img_save=True

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
            os.makedirs(f'{self.mainfolder}/img_{self.filename}', exist_ok=True)

            for i in range(len(phase_data.index)):
                os.makedirs(
                    f'{self.mainfolder}/img_{self.filename}/{subfolder}/{group_key}', exist_ok=True)
                # -------------------------------------------------------------------
                fig = plt.figure(constrained_layout=False,
                                facecolor='0.9', figsize=(9, 9))
                gs = fig.add_gridspec(nrows=3, ncols=3, left=0.15,
                                    right=0.95, hspace=0.6, wspace=0.5)

                ax4 = fig.add_subplot(gs[:, :])
                ax4.plot(unfiltered_porosity[:i+1]*100,
                        depth[:i+1]/1000, 'd--', color='green')
                if False in filtered_porosity.mask[:i+1]:
                    # plot 4
                    ax4.plot(
                        filtered_porosity[:i+1]*100, depth[:i+1]/1000, 'd--', color='black', alpha=0.7)
                    ax4.plot(filtered_porosity[i:i+1]*100, depth[i:i+1]/1000, 'd',
                            color='#7fffd4', markersize=8, markeredgecolor='black')
                else:
                    # plot 4
                    ax4.plot(np.ones(i+1)*1e-30, depth[:i+1]/1000, 'd')

                # Ax4 figure features
                ax4.set_title('T:{:.2f} °C P:{:.2f} GPa'.format(
                    ts[i], ps[i]/10000), fontsize=18)
                ax4.set_xlim(0.01, 10)
                ax4.set_xscale('log')
                ax4.set_ylim(100, 0)
                ax4.set_xlabel("Fluid-filled porosity [$vol.%$]", fontsize=18)
                ax4.set_ylabel("Depth [km]", fontsize=18)
                # ax4.set_xlim(10**0, 10**6)
                # ax4.fill_between([10**2.7 - 0.5, 10**2.7 + 0.5], [100, 100],
                #                  color="#c5c5c5", edgecolor='black', hatch="/", alpha=0.4)
                ax4.tick_params(axis='both', labelsize=20)
                # ax4.annotate("Ague 2003 Comp", (10**2.7, 15), fontsize=16)

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
            ax4 = fig.add_subplot(gs[:, :])

            # looping the array
            for i in range(len(phase_data.index)):
                ax4.plot(unfiltered_porosity[:i+1]*100,
                        depth[:i+1]/1000, 'd--', color='green')
                if False in filtered_porosity.mask[:i+1]:
                    # plot 4
                    ax4.plot(
                        filtered_porosity[:i+1]*100, depth[:i+1]/1000, 'd--', color='black', alpha=0.7)
                    ax4.plot(filtered_porosity[i:i+1]*100, depth[i:i+1]/1000, 'd',
                            color='#7fffd4', markersize=8, markeredgecolor='black')
                else:
                    # plot 4
                    ax4.plot(np.ones(i+1)*1e-30, depth[:i+1]/1000, 'd')

            # Ax4 figure features
            ax4.set_title('T:{:.2f} °C P:{:.2f} GPa'.format(
                ts[i], ps[i]/10000), fontsize=18)
            ax4.set_xlim(0.01, 10)
            ax4.set_xscale('log')
            ax4.set_ylim(100, 0)
            ax4.set_xlabel("Fluid-filled porosity [$vol.%$]", fontsize=18)
            ax4.set_ylabel("Depth [km]", fontsize=18)
            # ax4.set_xlim(10**0, 10**6)
            # ax4.fill_between([10**2.7 - 0.5, 10**2.7 + 0.5], [100, 100],
            #                  color="#c5c5c5", edgecolor='black', hatch="/", alpha=0.4)
            ax4.tick_params(axis='both', labelsize=20)
            # ax4.annotate("Ague 2003 Comp", (10**2.7, 15), fontsize=16)
            plt.show()

        # reading images and put it to GIF
        # reading images and put it to GIF
        if gif_save is True:
            # call gif function
            create_gif(phase_data, self.mainfolder, self.filename,
                       group_key, subfolder=subfolder)

    def release_fluid_volume_plot(self, rock_tag, img_save=False, gif_save=False):

        if gif_save is True:
            img_save=True

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

                ax4 = fig.add_subplot(gs[:, :])
                ax4.plot(unfiltered_v_fluid_extr[:i+1],
                         depth[:i+1]/1000, 'd--', color='green')
                if False in filtered_permeability.mask[:i+1]:
                    ax4.plot(filtered_v_fluid_extr[:i+1] /
                             1_000_000, depth[:i+1]/1000, 'd--')
                    ax4.plot(filtered_v_fluid_extr[:i+1]/1_000_000,
                             depth[:i+1]/1000, 'd--', color='black', alpha=0.5)
                    ax4.plot(filtered_v_fluid_extr[i:i+1]/1_000_000, depth[i:i+1]/1000,
                             'd--', color='#7fffd4', markersize=8, markeredgecolor='black')

                # Ax4 figure features
                ax4.set_title("Extracted fluid volume")
                ax4.set_xlim(0.001, 100)
                ax4.set_xscale('log')
                ax4.set_ylim(100, 0)
                # TODO ccm because it is not the geometry scale
                ax4.set_xlabel("Volume extracted [$m^{3}$]", fontsize=12)
                ax4.set_ylabel("Depth [km]", fontsize=12)
                ax4.tick_params(axis='both', labelsize=14)

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
            ax4 = fig.add_subplot(gs[:, :])
            for i in range(len(phase_data.index)):
                ax4.plot(unfiltered_v_fluid_extr[:i+1],
                         depth[:i+1]/1000, 'd--', color='green')
                if False in filtered_permeability.mask[:i+1]:
                    ax4.plot(filtered_v_fluid_extr[:i+1] /
                             1_000_000, depth[:i+1]/1000, 'd--')
                    ax4.plot(filtered_v_fluid_extr[:i+1]/1_000_000,
                             depth[:i+1]/1000, 'd--', color='black', alpha=0.5)
                    ax4.plot(filtered_v_fluid_extr[i:i+1]/1_000_000, depth[i:i+1]/1000,
                             'd--', color='#7fffd4', markersize=8, markeredgecolor='black')

            # Ax4 figure features
            ax4.set_title("Extracted fluid volume")
            ax4.set_xlim(0.001, 100)
            ax4.set_xscale('log')
            ax4.set_ylim(100, 0)
            # TODO ccm because it is not the geometry scale
            ax4.set_xlabel("Volume extracted [$m^{3}$]", fontsize=12)
            ax4.set_ylabel("Depth [km]", fontsize=12)
            ax4.tick_params(axis='both', labelsize=14)
            plt.show()

        # reading images and put it to GIF
        if gif_save is True:
            # call gif function
            create_gif(phase_data, self.mainfolder, self.filename,
                       group_key, subfolder=subfolder)

    def phases_stack_plot(self, rock_tag, img_save=False):

        # XMT naming and coloring
        database = self.rockdic[rock_tag].database
        phases2 = self.rockdic[rock_tag].phases2
        legend_phases, color_set = phases_and_colors_XMT(database, phases2)

        phases = self.rockdic[rock_tag].phases
        phase_data = self.rockdic[rock_tag].phase_data
        frac_bool = self.rockdic[rock_tag].frac_bool

        # Input for the variable of interest
        tag_in = input("Please provide what you want to convert to a stack. ['vol%', 'volume[ccm]', 'wt%', 'wt[g]']")
        tag='df_'+tag_in

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
        label_list = legend_phases

        # y.columns = phases
        # y = Merge_phase_group(y)
        # label_list = list(y.index)
        # pal1 = sns.color_palette("tab20", 20)

        # plot
        plt.rc('axes', labelsize=10)
        plt.rc('xtick', labelsize=6)  # fontsize of the x tick labels
        plt.rc('ytick', labelsize=8)  # fontsize of the y tick labels
        ax1 = host_subplot(111)
        #fig.suptitle('Phase changes for P-T-t')
        # main plot


        # ax1.stackplot(line, y, labels=label_list,
        #               colors=pal1, alpha=.8, edgecolor='black')
        ax1.stackplot(line, y, labels=label_list,
                      colors=color_set, alpha=.8, edgecolor='black')

        # legend definition
        handles, labels = ax1.get_legend_handles_labels()
        ax1.legend(
            handles[::-1], labels[::-1], bbox_to_anchor=(1.43, 0.5), loc='right',
            borderaxespad=0.1, title="Stable phases", fontsize=8
        )

        # tick stepping
        con = round(10/len(line), 2)
        con = round(len(line)*con)
        step = round(len(line)/con)
        if step == 0:
            step = 1

        # define x-axis bottom (temperature) tick labels
        ax1.set_xticks(line[::step])

        if len(line) < len(temperatures):
            selected_temp = temperatures[list(line-1)]
            ax1.set_xticklabels(np.around(selected_temp[::step], 1))
        else:
            ax1.set_xticklabels(np.around(temperatures[::step], 1))
        ax1.xaxis.set_minor_locator(AutoMinorLocator())

        # second y axis
        # free fluid content
        twin1 = ax1.twinx()
        y2 = (st_fluid_before)/system_vol_pre*100
        # extraction steps
        # NOTE extraction marker boolean
        mark_extr = extraction_boolean = self.rockdic[rock_tag].extraction_boolean
        mark_extr = np.array(system_vol_pre-system_vol_post, dtype='bool')
        twin1.plot(line, y2, '1--', c='blue')
        if len(frac_bool) > 0:
            if 1 in frac_bool or 2 in frac_bool:
                extension_bool = np.isin(frac_bool, 1)
                shear_bool = np.isin(frac_bool, 2)
                twin1.plot(line[extension_bool], y2[extension_bool], 'Dr')
                twin1.plot(line[shear_bool], y2[shear_bool], 'Dg')
        else:
            pass
        twin1.set_ylabel("Vol% of fluid-filled porosity")
        twin1.set_ymargin(0)

        # define x-axis top (pressure) tick labels
        ax2 = ax1.twin()
        ax2.set_xticks(line[::step])

        if len(line) < len(pressures):
            selected_pres = pressures[list(line-1)]
            ax1.set_xticklabels(np.around(selected_pres[::step], 1))
        else:
            ax2.set_xticklabels(np.around(np.array(pressures[::step])/1000, 1))

        # labeling and style adjustments
        plt.subplots_adjust(right=0.75)
        pressures = list(pressures)
        temperatures = list(temperatures)
        peakp = pressures.index(max(pressures))
        peakt = temperatures.index(max(temperatures))
        if len(temperatures) > peakp+1:
            ax1.axvline(x=line[peakp+1], linewidth=1.5,
                        linestyle='--', color='black')
            ax1.axvline(x=line[peakp], linewidth=1.5,
                        linestyle='--', color='black')
            ax1.axvline(x=line[peakt], linewidth=1.5,
                        linestyle='--', color='red')
        ax1.set_ylabel(tag[3:])
        ax1.set_xlabel('Temperature [°C]')
        ax2.set_xlabel('Pressure in kBar (1$e^{3}$)')
        ax2.axis["right"].major_ticklabels.set_visible(False)
        ax2.axis["right"].major_ticks.set_visible(False)


        # save image
        if img_save is True:
            plt.savefig(Path(self.mainfolder /
                        f"{self.rock_key}_new_stack_plot.png"), dpi=300)
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

        whitelist_attributes = ['temperature', 'pressure', 'depth', 'systemVolPre', 'systemVolPost', 'permeability', 'extracted_fluid_volume', 'porosity', 'time_int_flux2']
        whitelist_phase_data = ['N', 'vol%', 'volume[ccm]', 'wt%', 'wt[g]']

        if x_ax_label in attributes and x_ax_label in whitelist_attributes:
            if x_ax_label == 'temperature':
                xdata = data.temperature
            if x_ax_label == 'pressure':
                xdata = data.pressure
            if x_ax_label == 'depth':
                xdata = data.depth
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
            xdata = input(f"Your x data query requires a phase input from the following list:\n{sel_d.columns}")
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
            ydata = input(f"Your x data query requires a phase input from the following list:\n{sel_d.columns}")
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
        # XMT naming and coloring
        database = self.rockdic[rock_tag].database
        phases2 = self.rockdic[rock_tag].phases2
        legend_phases, color_set = phases_and_colors_XMT(database, phases2)

        # Plotting routine
        summary = self.rockdic[rock_tag].oxygen_data.T.describe()
        # oxyframe = Merge_phase_group(self.rockdic[rock_tag].oxygen_data)
        oxyframe = self.rockdic[rock_tag].oxygen_data.T
        oxyframe.columns = self.rockdic[rock_tag].temperature
        oxyframe.index = legend_phases
        fig, ax111 = plt.subplots(1, 1, figsize=(8, 5))
        for t, phase in enumerate(list(oxyframe.index)):
            ax111.plot(oxyframe.columns, oxyframe.loc[phase], '--d', color=color_set[t], linewidth=0.7, markeredgecolor='black')


        ax111.plot(self.rockdic[rock_tag].temperature, self.rockdic[rock_tag].bulk_deltao_post, '-', c='black')
        ax111.plot(self.rockdic[rock_tag].temperature, self.rockdic[rock_tag].bulk_deltao_pre, '-.', c='black')
        legend_list = list(oxyframe.index) + ["pre bulk", "post bulk"]
        ax111.legend(legend_list, bbox_to_anchor=(1.28, 0.9))
        min_min = min(summary.loc['min'])
        max_max = max(summary.loc['max'])
        min_val = min_min - (max_max-min_min)*0.05
        max_val = max_max + (max_max-min_min)*0.05
        ax111.set_ylim(min_val, max_val)
        ax111.set_ylabel("$\delta^{18}$O [‰ vs. VSMOW]")
        ax111.set_xlabel("Temperature [°C]")
        plt.subplots_adjust(right=0.8)

        if img_save is True:
            plt.savefig(Path(self.mainfolder /
                            f"{self.rock_key}_oxygen_signatures.png"), dpi=300)
        else:
            plt.show()
        plt.clf()
        plt.close()


if __name__ == '__main__':

    # Read variables from a hdf5 output file from ThorPT
    data = ThorPT_hdf5_reader()
    data.open_ThorPT_hdf5()

    # Activate the plotting module - Predefined plots for evaluation
    compPlot = ThorPT_plots(
        data.filename, data.mainfolder, data.rock, data.compiledrock)

    print("The 'data.rock[key]' keys are:")
    for key in data.rock.keys():
        print(key)

    compPlot.boxplot_to_GIF(rock_tag='rock1', img_save=True, gif_save=True)
    compPlot.phases_stack_plot(rock_tag='rock0', img_save=False)
    compPlot.oxygen_isotopes(rock_tag='rock0')
    compPlot.binary_plot(rock_tag='rock0')
    compPlot.pt_path_plot(rock_tag='rock0', img_save=False, gif_save=False)
    compPlot.permeability_plot(rock_tag='rock0', img_save=False, gif_save=False)
    compPlot.time_int_flux_plot(rock_tag='rock0', img_save=False, gif_save=False)
    compPlot.porosity_plot(rock_tag='rock0', img_save=False, gif_save=False)
    compPlot.release_fluid_volume_plot(rock_tag='rock0', img_save=False, gif_save=False)






