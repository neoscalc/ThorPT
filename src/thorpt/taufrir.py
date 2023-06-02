
"""
Written by
Thorsten Markmann
thorsten.markmann@geo.unibe.ch
status: 09.02.2023
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import h5py
import os
import math

from pathlib import Path
from mpl_toolkits.axes_grid1 import host_subplot
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
from tkinter import filedialog
import imageio
import copy


def set_origin():
    dirname = os.path.dirname(os.path.abspath(__file__))
    os.chdir(dirname)


def Cumulative(lists):
    cu_list = []
    length = len(lists)
    cu_list = [sum(lists[0:x:1]) for x in range(0, length+1)]
    return cu_list[1:]


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


def file_opener():
    filein = filedialog.askopenfilename(
        title="Select h5-file to read",
        filetypes=(
            ("hdf5 file", "*.hdf5"),
            ("All files", "*.*"))
    )

    return filein


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

# ANCHOR mineral translation - get database names
# Name translation Database to XMT-Colour Names


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
    main_folder = Path(__file__).parent.absolute()
    data_folder = main_folder / "DataFiles/"

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


class Simple_plot():

    def __init__(self, temperatures, pressures) -> None:
        self.temperatures = temperatures
        self.pressures = pressures

    def path_plot(self):
        plt.scatter(self.temperatures, self.pressures/1000, s=15)
        plt.xlabel("Temperature in DC")
        plt.ylabel("Pressure in kBar")
        plt.savefig(Path("plotting_results/PTt_path.png"), dpi=400)

# ANCHOR - Plot master


class Plot_master():

    def __init__(self, temperatures, pressures, line, time_list, superior_data, var_dic, td_data, ext_data, physc_data, trackv, norm, hyd_data, rock_key, fluid_before, phases, tensile, live_plot=False):
        """[summary]

        Args:
            temperatures ([type]): [description]
            pressures ([type]): [description]
            time_list ([type]): [description]
            superior_data ([list]): [description]
        """

        self.temperatures = temperatures
        self.pressures = pressures
        # generate line array for x axis
        if line is False:
            self.line = np.arange(1, len(temperatures)+1, 1)
        else:
            self.line = line
        self.time_passed = time_list
        self.svol_data = physc_data
        self.ext_fluid = superior_data[2]
        self.physc_prop_frame = var_dic
        self.tracking_vol = fluid_before
        self.normalization = superior_data[4]
        self.hydrate_data = superior_data[5]
        self.live_plot = live_plot
        self.rock_key = rock_key
        self.phases = phases
        self.tensile = tensile

    def physical_system_and_fluid(self):
        # REVIEW second x axis with pressures has offset.......

        temperatures = self.temperatures
        pressures = self.pressures
        sys_dat = self.td_data
        line = np.arange(1, len(self.temperatures)+1, 1)

        system_density_pre = sys_dat[2][0]
        system_density_post = sys_dat[2][1]
        system_vol_pre = sys_dat[0][0]
        system_vol_post = sys_dat[0][1]
        system_weight_pre = sys_dat[1][0]
        system_weight_post = sys_dat[1][1]

        master_norm = self.normalization

        df_var_dictionary = self.physc_prop_frame

        extracted_fluid_data = self.ext_fluid

        if extracted_fluid_data.empty:
            tags = ['N', 'volume/mol', 'volume[ccm]', 'vol%',
                    'wt/mol', 'wt[g]', 'wt%', 'density[g/ccm]']
            extracted_fluid_data.index = tags

        # st_solid, st_fluid_before, st_fluid_after = self.tracking_vol[
        #     0], self.tracking_vol[1], self.tracking_vol[2]
        st_fluid_before, st_fluid_after = self.tracking_vol[1], self.tracking_vol[2]

        # preparing differentiate and derivatives for system volume, density, weight
        volume_ar = np.array(system_vol_pre)
        volume_deriv_pre = np.diff(volume_ar, prepend=volume_ar[0])
        volume_ar = np.array(system_vol_post)
        volume_deriv_post = np.diff(volume_ar, prepend=volume_ar[0])
        norm_vol_deriv_pre = volume_deriv_pre/master_norm
        norm_vol_deriv_post = volume_deriv_post/master_norm

        dat = np.array(system_weight_pre)
        norm_weight_deriv_pre = np.diff(dat, prepend=dat[0])/master_norm
        dat = np.array(system_weight_post)
        norm_weight_deriv_post = np.diff(dat, prepend=dat[0])/master_norm

        dat = np.array(system_density_pre)
        density_change_pre = np.diff(dat, prepend=dat[1])
        dat = np.array(system_density_post)
        density_change_post = np.diff(dat, prepend=dat[1])
        norm_dens_deriv_pre = density_change_pre/master_norm
        norm_dens_deriv_post = density_change_post/master_norm

        # preparing hydrous phase data
        df_h2o_content_dic = self.hydrate_data
        # TODO - change tag
        density_frame = df_var_dictionary['df_density[g/ccm]']
        if 'df_H2O[g]' in df_h2o_content_dic:
            df_h2o_content_dic['df_H2O[g]'].fillna(0, inplace=True)
            hydrous_phases = list(df_h2o_content_dic['df_H2O[g]'].columns)
            hydrates = True
        else:
            hydrates = False
            print("warning: No hydrous phases to concenate")
        temp_dic = {}
        if hydrates is True:
            for item in hydrous_phases:
                if item == 'water.fluid':
                    pass
                else:
                    temp_dic[item] = density_frame.loc[:, item]
            hydrous_phases_dens_frame = pd.concat(
                [pd.Series(v, name=k) for k, v in temp_dic.items()], axis=1)

        """
        check_phase = hydrous_phases_dens_frame.columns
        if 'PHNG_fcel' in check_phase and 'PHNG_mu' in check_phase:
            # hydrous phases and their weight and volume evolution - restored for plotting
            muscv_vol = np.array(
                hydrous_phases_dens_frame['PHNG_fcel']+hydrous_phases_dens_frame['PHNG_mu'])
            muscv_vol_change = np.diff(
                muscv_vol, prepend=muscv_vol[0])/master_norm
            muscv_h2o_weight = np.array(
                df_h2o_content_dic['df_H2O[g]']['PHNG_fcel']+df_h2o_content_dic['df_H2O[g]']['PHNG_mu'])
            muscv_h2o_change = np.diff(
                muscv_h2o_weight, prepend=muscv_h2o_weight[0])/master_norm
        else:
            muscv_h2o_weight = np.zeros(len(temperatures))

        if 'lawsonite' in check_phase:
            laws_vol = np.array(hydrous_phases_dens_frame['lawsonite'])
            laws_vol_change = np.diff(
                laws_vol, prepend=laws_vol[0])/master_norm
            laws_h2o_weight = np.array(
                df_h2o_content_dic['df_H2O[g]']['lawsonite'])
            laws_h2o_change = np.diff(
                laws_h2o_weight, prepend=laws_h2o_weight[0])/master_norm
        else:
            laws_h2o_weigh = np.zeros(len(temperatures))

        if 'ClAMP_gl2' in check_phase:
            glauc_vol = np.array(hydrous_phases_dens_frame['ClAMP_gl2'])
            glauc_vol_change = np.diff(
                glauc_vol, prepend=glauc_vol[0])/master_norm
            glauc_h2o_weight = np.array(
                df_h2o_content_dic['df_H2O[g]']['ClAMP_gl2'])
            glauc_h2o_change = np.diff(
                glauc_h2o_weight, prepend=glauc_h2o_weight[0])/master_norm
        else:
            glauc_h2o_weight = np.zeros(len(temperatures))

        try:
            struct_water = muscv_h2o_weight + laws_h2o_weight + glauc_h2o_weight
        except UnboundLocalError:
            struct_water = np.zeros(len(temperatures))
        """

        stablething = hydrous_phases_dens_frame.columns
        water_dat = []
        for item in stablething:
            temp_weight_arr = np.array(df_h2o_content_dic['df_H2O[g]'][item])
            water_dat.append(temp_weight_arr)
        struct_water = np.zeros(len(water_dat[0]))
        for item in water_dat:
            struct_water = struct_water + item
        #########################################################################################################
        #########################################################################################################
        #########################################################################################################
        #########################################################################################################
        fig_physic, ((ax1, ax2), (ax3, ax4)) = plt.subplots(
            2, 2, figsize=(15, 10), constrained_layout=True)
        ###
        # d1 - main
        twin1 = ax1.twinx()
        ax1.scatter(line, norm_vol_deriv_post *
                    1000, marker='x', color='#4b0082')
        ax1.scatter(line, norm_vol_deriv_pre *
                    1000, marker='+', color='red')
        twin1.plot(line, system_vol_post, '-', color='black')
        twin1.plot(line, system_vol_pre, '-', color='orange')

        # set botom x axis for temperature
        ax1.set_xticks(line[::6])
        ax1.set_xticklabels(temperatures[::6])

        # define x-axis top (pressure) tick labels
        ax12 = ax1.twiny()
        ax12.set_xticks(line[::6])
        ax12.set_xticklabels(np.around(np.array(pressures[::6])/1000, 1))

        # labels
        ax1.set_ylabel("Deriv. System volume [ccm x1000]")
        #ax1.set_xlabel("Temperature in °C")
        twin1.set_ylabel('Volume [ccm]')
        ax12.set_xlabel('Pressure in kBar(1$e^{3}$)')
        ax1.legend(
            ['Deriv. Vol post-ext', 'Deriv. Vol/mol pre-ext'], loc='lower left')
        twin1.legend(['Post-ext Vol', 'Pre-ext Vol'], loc='upper left')
        ax1.axhline(color='black', linestyle="--", linewidth=0.5)
        tab_list = np.array(list(system_vol_post) + list(system_vol_pre))
        min_val = min(tab_list) - (max(tab_list) - min(tab_list))*0.05
        max_val = max(tab_list) + (max(tab_list) - min(tab_list))*0.05
        twin1.set_ylim(min_val, max_val)
        twin1.set_yticks(np.linspace(min_val, max_val, 20))

        ###
        # d2
        twin1 = ax2.twinx()
        ax2.scatter(line, norm_dens_deriv_post *
                    10000, marker='x', color='#4b0082')
        ax2.scatter(line, norm_dens_deriv_pre *
                    10000, marker='+', color='red')
        twin1.plot(line, system_density_post, '-', color='black')
        twin1.plot(line, system_density_pre, '-', color='orange')

        # define x-axis bottom (temperature) tick labels
        ax2.set_xticks(line[::6])
        ax2.set_xticklabels(temperatures[::6])

        # define x-axis top (pressure) tick labels
        ax22 = ax2.twiny()
        ax22.set_xticks(line[::6])
        ax22.set_xticklabels(np.around(np.array(pressures[::6])/1000, 1))

        # labeling
        ax2.set_ylabel("Deriv. System density [g/ccm x10000]")
        # ax2.set_xlabel("Temperature in °C")
        twin1.set_ylabel('System density [g/ccm]')
        ax22.set_xlabel('Pressure in kBar(1$e^{3}$)')

        # legend
        twin1.legend(['Post-ext density', 'Pre-ext density'],
                     loc='lower right')
        ax2.legend(['Deriv. density post-ext',
                    'Deriv. density pre-ext'], loc='upper left')
        ax2.axhline(color='black', linestyle="--", linewidth=0.5)
        min_val = min(system_density_pre) - \
            (max(system_density_pre) - min(system_density_pre))*0.05
        max_val = max(system_density_post) + \
            (max(system_density_post) - min(system_density_post))*0.05
        twin1.set_ylim(min_val, max_val)

        ###
        # d3
        chache_zeros = system_weight_pre-system_weight_pre
        cum_to_fluid = extracted_fluid_data.loc['wt[g]'] + chache_zeros
        cum_to_fluid = cum_to_fluid.fillna(0)
        cum_extracted_fluid_weight = np.cumsum(np.array(cum_to_fluid))

        # main plot data
        twin1 = ax3.twinx()
        ax3.plot(line, cum_extracted_fluid_weight, '-', color='black')
        ax3.plot(line, np.array(system_weight_post)/1000 +
                 np.array(cum_extracted_fluid_weight)/1000, '-', color='orange')
        # TODO - change tag
        if 'water.fluid' in df_var_dictionary['df_wt[g]'].columns:
            pa1 = df_var_dictionary['df_wt[g]']['water.fluid']
            pa2 = extracted_fluid_data.loc['wt[g]']
            pa1[pa2.index] = pa1[pa2.index] - pa2
            ax3.plot(
                line, np.array(pa1), '-',
                color='red'
            )  # fluid weight post
        else:
            pa1 = []
        ax3.plot(line, struct_water, '--', color='green')
        twin1.plot(line, (st_fluid_before) /
                   system_vol_pre*100, '-.')
        twin1.plot(line, (st_fluid_after)/system_vol_post*100, '-.')
        twin1.invert_yaxis()

        # label x axes
        # define x-axis bottom (temperature) tick labels
        ax3.set_xticks(line[::6])
        ax3.set_xticklabels(temperatures[::6])
        # define x-axis top (pressure) tick labels
        ax32 = ax3.twiny()
        ax32.set_xticks(line[::6])
        ax32.set_xticklabels(np.around(np.array(pressures[::6])/1000, 1))

        # labeling
        ax3.set_ylabel("Weight in [g]")
        ax3.set_xlabel("Temperature in °C")
        ax3.legend(['Extr. cum. fluid', 'Extracted fluid + System extracted / 1kg',
                    'Free fluid present in sys', 'struc. bound water'], loc='center left')
        twin1.legend(['Porosity pre-ext', 'Porosity post-ext'],
                     loc='center left', bbox_to_anchor=(0, 0.65))
        twin1.set_ylabel("Porosity in [Vol%]")

        max_list = list(cum_extracted_fluid_weight) + \
            list(struct_water) + list(pa1)
        max_list = np.array(max_list)
        max_list[np.isnan(max_list)] = 0
        max_val = max(max_list) + (max(max_list) - min(max_list))*0.05
        ax3.set_ylim(-0.5, max_val)

        porosity = (st_fluid_before)/system_vol_pre*100
        max_val = max(porosity) + (max(porosity) - min(porosity))*0.05
        twin1.set_ylim(max_val, -0.05)
        ax3.axhline(color='black', linestyle="--", linewidth=0.5)
        twin1.axhline(color='black', linestyle="--", linewidth=0.5)

        ###
        # d4
        twin1 = ax4.twinx()
        ax4.scatter(line, norm_weight_deriv_post *
                    1000, marker='x', color='#4b0082')
        ax4.scatter(line, norm_weight_deriv_pre *
                    1000, marker='+', color='red')
        twin1.plot(line, system_weight_post, '-', color='black')
        twin1.plot(line, system_weight_pre, '-', color='orange')

        # label x axes
        # define x-axis bottom (temperature) tick labels
        ax4.set_xticks(line[::6])
        ax4.set_xticklabels(temperatures[::6])
        # define x-axis top (pressure) tick labels
        ax42 = ax4.twiny()
        ax42.set_xticks(line[::6])
        ax42.set_xticklabels(np.around(np.array(pressures[::6])/1000, 1))

        # labeling
        ax4.set_ylabel('Deriv. System weight [g x1000]')
        ax4.set_xlabel('Temperature in °C')
        twin1.set_ylabel('System weight [g]')
        ax4.legend(['Deriv. weight post-ext', 'Deriv. weight pre-ext'],
                   loc='lower left')
        ax4.axhline(color='black', linestyle="--", linewidth=0.5)
        twin1.legend(['Post-ext weight', 'Pre-ext weight'],
                     loc='lower center', bbox_to_anchor=(0.6, 0))

        min_val = min(system_weight_post) - \
            (max(system_weight_post) - min(system_weight_post))*0.05
        max_val = max(system_weight_pre) + \
            (max(system_weight_pre) - min(system_weight_pre))*0.05
        twin1.set_ylim(min_val, max_val)

        # saving
        main_folder = Path(__file__).parent.absolute()
        plt.savefig(Path(main_folder / "plotting_results" /
                    "physc_props.png"), dpi=400)

    def porosity_plot(self):
        temperatures = self.temperatures
        sys_dat = self.td_data
        system_vol_pre = sys_dat[0][0]
        system_vol_post = sys_dat[0][1]
        system_weight_pre = sys_dat[1][0]
        extracted_fluid_data = self.ext_fluid
        st_fluid_before = self.tracking_vol[1]
        max_vol = max(system_vol_pre)

        df_var_dictionary = self.physc_prop_frame
        # preparing hydrous phase data
        df_h2o_content_dic = self.hydrate_data
        # TODO - change tag
        density_frame = df_var_dictionary['df_density[g/ccm]']
        if 'df_H2O[g]' in df_h2o_content_dic:
            df_h2o_content_dic['df_H2O[g]'].fillna(0, inplace=True)
            hydrous_phases = list(df_h2o_content_dic['df_H2O[g]'].columns)
            hydrates = True
        else:
            hydrates = False
            print("warning: No hydrous phases to concenate")
        temp_dic = {}
        if hydrates is True:
            for item in hydrous_phases:
                if item == 'water.fluid':
                    pass
                else:
                    temp_dic[item] = density_frame.loc[:, item]
            hydrous_phases_dens_frame = pd.concat(
                [pd.Series(v, name=k) for k, v in temp_dic.items()], axis=1)

        stablething = hydrous_phases_dens_frame.columns
        water_dat = []
        for item in stablething:
            temp_weight_arr = np.array(df_h2o_content_dic['df_H2O[g]'][item])
            water_dat.append(temp_weight_arr)
        struct_water = np.zeros(len(water_dat[0]))
        for item in water_dat:
            struct_water = struct_water + item

        chache_zeros = system_weight_pre-system_weight_pre
        cum_to_fluid = extracted_fluid_data.loc['wt[g]'] + chache_zeros
        cum_to_fluid = cum_to_fluid.fillna(0)
        cum_extracted_fluid_weight = np.cumsum(np.array(cum_to_fluid))

        fig, (ax1, ax2) = plt.subplots(
            2, 1, figsize=(9, 9), constrained_layout=True)

        twin1 = ax1.twinx()
        ax1.plot(temperatures, cum_extracted_fluid_weight, '-', color='black')
        ax1.plot(temperatures, struct_water, '--', color='green')
        ax1.legend(
            ["Cumulative extracted fluid", "Structural bound water"],
            loc='upper center', bbox_to_anchor=(0.24, 1))
        twin1.plot(temperatures, (st_fluid_before) /
                   system_vol_pre*100, '-.', c='blue')
        twin1.legend(["Fluid filled prosity"],
                     loc='upper center', bbox_to_anchor=(0.2, 0.88))

        twin2 = ax2.twinx()
        twin2.plot(temperatures, (st_fluid_before) /
                   system_vol_pre*100, '-.', c='green')
        ax2.plot(temperatures, system_vol_pre/max_vol*100)
        ax2.plot(temperatures, system_vol_post/max_vol*100)
        ax2.legend(["Sys vol pre", "Sys vol post"])

        ax1.set_ylabel("weight [g]")
        ax2.set_xlabel("Temperature in °C")
        ax2.set_ylabel("Vol% normalized to start")
        twin1.set_ylabel("Vol%")

        main_folder = Path(__file__).parent.absolute()
        plt.savefig(
            Path(main_folder / "plotting_results" / "porosity.png"), dpi=400)

    def adv_stack_plot1_2(self, frac_bool, tag='df_vol%'):
        # compile data for plotting
        system_vol_pre = np.array(self.svol_data['system_vol_pre'])
        system_vol_post = np.array(self.svol_data['system_vol_post'])
        st_fluid_before = self.tracking_vol
        # max_vol = max(system_vol_pre)
        temperatures = self.temperatures
        pressures = self.pressures
        line = self.line

        y = self.physc_prop_frame[tag].fillna(value=0)
        y.columns = self.phases
        y = Merge_phase_group(y)
        label_list = list(y.index)
        pal1 = sns.color_palette("tab20", 20)

        # plot
        plt.rc('axes', labelsize=10)
        plt.rc('xtick', labelsize=6)  # fontsize of the x tick labels
        plt.rc('ytick', labelsize=8)  # fontsize of the y tick labels
        ax1 = host_subplot(111)
        #fig.suptitle('Phase changes for P-T-t')
        # main plot
        ax1.stackplot(line, y, labels=label_list,
                      colors=pal1, alpha=.8, edgecolor='black')

        # legend definition
        handles, labels = ax1.get_legend_handles_labels()
        ax1.legend(
            handles[::-1], labels[::-1], bbox_to_anchor=(1.4, 0.5), loc='right',
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

        plt.title(
            f"Tensile strength: {round(self.tensile)}", loc='left', y=0.935)

        # save image
        plt.savefig(Path(main_folder /
                    f"{self.rock_key}_new_stack_plot.png"), dpi=300)
        plt.clf()
        plt.close()

    def stack_plot(self, tag='df_wt[g]'):
        """Creating a stacked plot of every stable phase calculated over the P-T-t path

        Args:
            tag (str, optional): Can be any entry of the physical variable
                available as a dataframe. Defaults to 'df_volume[ccm]'.
        """

        pal1 = sns.color_palette("Paired")
        fig, ax1 = plt.subplots(1, 1, figsize=(8, 5))
        fig.suptitle('Phase changes for P-T-t')
        x = list(self.physc_prop_frame[tag].index)
        y = self.physc_prop_frame[tag].fillna(value=0)
        y = Merge_phase_group(y)

        if tag == "df_vol%":
            norm_cum = []
            for tt in y.columns:
                norm_cum.append(y[tt].sum())
            norm_cum = np.array(norm_cum)

            for i, col in enumerate(y.columns):
                y[col] = y[col] * 100/norm_cum[i]

        label_list = list(y.index)
        ax1.stackplot(x, y, labels=label_list, colors=pal1, alpha=1)

        handles, labels = ax1.get_legend_handles_labels()
        ax1.legend(handles[::-1], labels[::-1], bbox_to_anchor=(1.28, 0.5), loc='right',
                   borderaxespad=0.1, title="Stable phases")
        ax1.set_xlabel('Temperature [°C]')
        ax1.set_ylabel(tag)
        plt.subplots_adjust(right=0.8)
        main_folder = Path(__file__).parent.absolute()
        plt.savefig(Path(main_folder / "plotting_results" /
                    "stack_plot.png"), dpi=300)

    def oxy_iso_plot(self, oxygen_data_frame, save_bulk_oxygen_pre, save_bulk_oxygen_post):

        # isotope plot
        summary = oxygen_data_frame.describe()
        # oxyframe = Merge_phase_group_oxy(oxygen_data_frame.T)
        oxyframe = oxygen_data_frame.T

        # generate line array for x axis
        line = self.line

        # colors
        pal1 = sns.color_palette("tab20", 20)

        # main plot
        ax1 = host_subplot(111)
        i = 0
        for phase in list(oxyframe.index):
            ax1.plot(oxyframe.columns,
                     oxyframe.loc[phase], '-d', alpha=0.6, markersize=1, c=pal1[i])
            i = i+1
        ax1.plot(line, save_bulk_oxygen_pre, '-b')
        ax1.plot(line, save_bulk_oxygen_post, 'xb')

        # legend
        legend_list = list(oxyframe.index) + ["pre bulk", "post bulk"]
        ax1.legend(legend_list, bbox_to_anchor=(1.01, 1), fontsize='x-small')

        con = round(10/len(line), 2)
        con = round(len(line)*con)
        step = round(len(line)/con)
        # define x-axis bottom (temperature) tick labels
        ax1.set_xticks(line[::step])
        ax1.set_xticklabels(np.around(self.temperatures[::step], 1))

        # define x-axis top (pressure) tick labels
        ax2 = ax1.twin()
        ax2.set_xticks(line[::step])
        ax2.set_xticklabels(
            np.around(np.array(self.pressures[::step])/10000, 1))

        # labeling and style
        min_min = min(summary.loc['min'])
        max_max = max(summary.loc['max'])
        min_val = min_min - (max_max-min_min)*0.05
        max_val = max_max + (max_max-min_min)*0.05
        ax1.set_ylim(min_val, max_val)
        ax2.axis["right"].major_ticklabels.set_visible(False)
        ax1.set_ylabel("$\delta^{18}$O [‰ vs. VSMOW]")
        ax1.set_xlabel("Temperature [°C]")
        ax2.set_xlabel('Pressure [$GPa$]')
        plt.rc('xtick', labelsize=6)
        plt.subplots_adjust(right=0.8)

        for i, phase in enumerate(oxyframe.index):
            ax1.annotate(legend_list[i], [
                self.temperatures[-1], oxyframe.iloc[:, -1][phase]])

        plt.savefig(Path(main_folder /
                    f"{self.rock_key}_d18O_plot.png"), dpi=300)
        plt.clf()
        plt.close()

    # Before the change to loop P-T paths
    def oxy_iso_plot_backup220406(self, oxygen_data_frame, save_bulk_oxygen_pre, save_bulk_oxygen_post):

        # isotope plot
        summary = oxygen_data_frame.T.describe()
        # oxyframe = Merge_phase_group_oxy(oxygen_data_frame.T)
        oxyframe = oxygen_data_frame
        line = np.arange(1, len(self.temperatures)+1, 1)
        # main plot
        pal1 = sns.color_palette("tab20", 20)
        fig, ax111 = plt.subplots(1, 1, figsize=(8, 5))
        i = 0
        for phase in list(oxyframe.index):
            ax111.plot(oxyframe.columns,
                       oxyframe.loc[phase], '-d', alpha=0.6, markersize=1, c=pal1[i])
            i = i+1
        ax111.plot(line, save_bulk_oxygen_pre[1:], '-b')
        ax111.plot(line, save_bulk_oxygen_post, 'xb')

        # legend
        legend_list = list(oxyframe.index) + ["pre bulk", "post bulk"]
        ax111.legend(legend_list, bbox_to_anchor=(1.28, 0.9))

        # labeling and style
        min_min = min(summary.loc['min'])
        max_max = max(summary.loc['max'])
        min_val = min_min - (max_max-min_min)*0.05
        max_val = max_max + (max_max-min_min)*0.05
        ax111.set_ylim(min_val, max_val)
        ax111.set_ylabel("$\delta^{18}$O [‰ vs. VSMOW]")
        ax111.set_xlabel("Temperature [°C]")
        plt.subplots_adjust(right=0.8)

        for i, phase in enumerate(oxyframe.index):
            ax111.annotate(legend_list[i], [
                           self.temperatures[-1], oxyframe.iloc[:, -1][phase]])

        main_folder = Path(__file__).parent.absolute()
        plt.savefig(Path(main_folder / "plotting_results" /
                    "oxygen_plot.png"), dpi=300)

    def oxy2_plot(self, oxygen_data_frame, save_bulk_oxygen_pre, save_bulk_oxygen_post):
        # isotope plot
        summary = oxygen_data_frame.T.describe()
        oxyframe = Merge_phase_group(oxygen_data_frame.T)
        fig, ax111 = plt.subplots(1, 1, figsize=(8, 5))
        for phase in list(oxyframe.index):
            ax111.plot(
                oxyframe.columns,
                oxyframe.loc[phase], '--', alpha=0.6
            )
        ax111.plot(self.temperatures, save_bulk_oxygen_post, '-', c='black')
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

        main_folder = Path(__file__).parent.absolute()
        plt.savefig(Path(main_folder / "plotting_results" /
                    "oxygen_plot2.png"), dpi=300)

    def ol_atg_oxygen(self, oxygen_data_frame):
        interest_data = pd.DataFrame()
        names = []
        phases = list(oxygen_data_frame.index)
        # selection of important phases to plot
        for phase in phases:
            if 'OLIVINEi_fo' in phase:
                names.append(phase)
                interest_data = pd.concat(
                    [interest_data, oxygen_data_frame.loc[phase]], axis=1)
            if 'atg' in phase:
                names.append(phase)
                interest_data = pd.concat(
                    [interest_data, oxygen_data_frame.loc[phase]], axis=1)
            if 'BR_br' in phase:
                names.append(phase)
                interest_data = pd.concat(
                    [interest_data, oxygen_data_frame.loc[phase]], axis=1)
            if 'water.fluid' in phase:
                names.append(phase)
                interest_data = pd.concat(
                    [interest_data, oxygen_data_frame.loc[phase]], axis=1)
        interest_data = interest_data.T
        interest_data.columns = oxygen_data_frame.columns

        # merging olivine and using diff of Olivine-Atg
        i = 0
        for name in names:
            if 'OLIVINE' in name:
                arr = interest_data.loc[name]
                arr = arr.fillna('A')
                if i == 0:
                    olivine_arr = arr
                else:
                    olivine_arr = pd.concat(
                        [olivine_arr, arr], axis=1).sum(axis=1)
        olivine_arr = pd.to_numeric(olivine_arr, errors='coerce')
        delta_ol_atg = np.array(olivine_arr) - \
            np.array(interest_data.loc['SERP_atg'])

        fig, ax20 = plt.subplots(1, 1, figsize=(9, 5))
        for phase in list(interest_data.index):
            ax20.plot(interest_data.columns,
                      interest_data.loc[phase], 'd-', alpha=0.6)
        ax20.plot(interest_data.columns, delta_ol_atg, '-.', c='black')
        ax20.legend(list(interest_data.index)+['diff ol-atg'])
        ax20.set_ylabel("$\delta^{18}$O [‰ vs. VSMOW]")
        ax20.set_xlabel("Temperatue [°C]")
        main_folder = Path(__file__).parent.absolute()
        plt.savefig(Path(main_folder / "plotting_results" /
                    "olivine_and_antg.png"), dpi=300)

    def chemp_plot(self, data):
        # chemical potential data
        f_potential = data
        # state numbers of rows if 3 plots per row
        num = math.ceil(len(f_potential.index)/3)
        # subplot start
        fig, axs = plt.subplots(num, 3, figsize=(
            15, 10), constrained_layout=True)

        ticks = self.line

        for i, item in enumerate(f_potential.index):
            if i < 3:
                test_arr = np.array(f_potential.loc[item])
                axs[0, i].plot(self.line, test_arr, 'x--')
                if min(test_arr) < 0.8*np.nanmedian(test_arr) or max(test_arr) < 1.2*np.nanmedian(test_arr):
                    pass
                else:
                    axs[0, i].set_ylim(
                        [0.8*np.nanmedian(test_arr), 1.2*np.nanmedian(test_arr)])
                axs[0, i].set_title(f"{item} potential in the system")
                axs[0, i].set_xlabel("Temperature in DC")
                axs[0, i].set_ylabel("Chem. potential in kJ/mol")
            if i >= 3 and i < 6:
                test_arr = np.array(f_potential.loc[item])
                axs[1, i-3].plot(self.line, test_arr, 'x--')
                if min(test_arr) < 0.8*np.nanmedian(test_arr) or max(test_arr) < 1.2*np.nanmedian(test_arr):
                    pass
                else:
                    axs[1, i-3].set_ylim(
                        [0.8*np.nanmedian(test_arr), 1.2*np.nanmedian(test_arr)])
                axs[1, i-3].set_title(f"{item} potential in the system")
                axs[1, i-3].set_xlabel("Temperature in DC")
                axs[1, i-3].set_ylabel("Chem. potential in kJ/mol")
            if i >= 6 and i < 9:
                test_arr = np.array(f_potential.loc[item])
                axs[2, i-6].plot(self.line, test_arr, 'x--')
                if min(test_arr) < 0.8*np.nanmedian(test_arr) or max(test_arr) < 1.2*np.nanmedian(test_arr):
                    pass
                else:
                    axs[2, i-6].set_ylim(
                        [0.8*np.nanmedian(test_arr), 1.2*np.nanmedian(test_arr)])
                axs[2, i-6].set_title(f"{item} potential in the system")
                axs[2, i-6].set_xlabel("Temperature in DC")
                axs[2, i-6].set_ylabel("Chem. potential in kJ/mol")
            if i >= 9:
                test_arr = np.array(f_potential.loc[item])
                axs[3, i-9].plot(self.line, test_arr, 'x--')
                if min(test_arr) < 0.8*np.nanmedian(test_arr) or max(test_arr) < 1.2*np.nanmedian(test_arr):
                    pass
                else:
                    axs[3, i-9].set_ylim(
                        [0.8*np.nanmedian(test_arr), 1.2*np.nanmedian(test_arr)])
                axs[3, i-9].set_title(f"{item} potential in the system")
                axs[3, i-9].set_xlabel("Temperature in DC")
                axs[3, i-9].set_ylabel("Chem. potential in kJ/mol")

        for ax in axs.flat:
            # set the tick numbers for each subplot
            ax.set_xticks(ticks=ticks)

        plt.setp(axs, xticklabels=np.round(self.temperatures[::5], 1))

        plt.savefig(Path(main_folder /
                         f"{self.rock_key}_potentials.png"), dpi=400)
        plt.clf()
        plt.close()

    def compare_extraction(self, line, all_system_vol_pre, all_system_vol_post, all_fluid_before, all_tensile, all_diffs, all_frac_bool, permea):

        colpal = sns.color_palette("flare", int(max(all_tensile)+1))
        fig, ax3000 = plt.subplots(
            1, 5, figsize=(20, 8), constrained_layout=True)
        tw1 = ax3000[1].twiny()
        tw13 = ax3000[2].twinx()

        for i, item in enumerate(all_system_vol_pre):
            # load system volume arrays
            system_vol_pre = np.array(all_system_vol_pre[i])
            system_vol_post = np.array(all_system_vol_post[i])
            st_fluid_before = all_fluid_before[i]

            # fluid filled porosity array
            y2 = (st_fluid_before)/system_vol_pre*100
            if isinstance(line[i], int):
                if line[i] == False:
                    line[i] = np.arange(1, len(self.temperatures)+1, 1)

            # set ticks
            ticks = np.arange(1, len(self.temperatures)+1, 1)[::3]
            # select extraction steps
            mark_extr = np.array(system_vol_pre-system_vol_post, dtype='bool')
            # porosity line plot
            ax3000[0].plot(line[i], y2, '1--',
                           color=colpal[int(all_tensile[i])])
            ax3000[0].plot(line[i][mark_extr], y2[mark_extr], 'Dr')
            ax3000[0].annotate(str(all_tensile[i]), [line[i][-1], y2[-1]])
            # number of extraction vs tensile strength
            ax3000[1].plot(len(y2[mark_extr]), all_tensile[i],
                           'D',  color=colpal[int(all_tensile[i])])
            if len(y2[mark_extr]) > 0:
                tw1.plot(max(y2[mark_extr]), all_tensile[i], 'Xb')
            else:
                tw1.plot(0, all_tensile[i], 'Db')

            if len(arr_line[i]) < len(self.temperatures):
                selected_temperature = self.temperatures[arr_line[i]-1]
                ax3000[2].plot(np.round(selected_temperature), all_diffs[i],
                               'x--', color=colpal[int(all_tensile[i])])
            else:
                ax3000[2].plot(np.round(self.temperatures), all_diffs[i],
                               'x--', color=colpal[int(all_tensile[i])])
            if np.size(permea[i]) == 1 and np.sum(permea[i]) == 0:
                pass
            else:

                # ANCHOR new routine for masked array technique
                tw13.plot(self.temperatures, permea, 'sb')

                """
                if len(mark_extr) < len(self.temperatures):
                    selected_temperature = self.temperatures[arr_line[i]-1]
                    if len(permea[i]) > len(np.round(selected_temperature)[mark_extr]):
                        tw13.plot(np.round(selected_temperature)[
                            mark_extr], permea[i][permea[i] != 0], 'sb')
                    else:
                        tw13.plot(np.round(selected_temperature)[
                            mark_extr], permea[i], 'sb')
                else:
                    if len(permea[i]) > len(np.round(self.temperatures)[mark_extr]):
                        tw13.plot(np.round(self.temperatures)[
                            mark_extr], permea[i][permea[i] != 0], 'sb')
                    else:
                        tw13.plot(np.round(self.temperatures)[
                            mark_extr], permea[i], 'sb')
                """

            ax3000[3].plot(all_diffs[i], y2, 'd',
                           color=colpal[int(all_tensile[i])])
            ax3000[3].plot(all_diffs[i][mark_extr],
                           y2[mark_extr], 'xr', markersize=15)
            y22 = np.ones(len(all_diffs[i]))*all_tensile[i]
            ax3000[4].plot(all_diffs[i], y22, 'd',
                           color=colpal[int(all_tensile[i])])
            ax3000[4].plot(all_diffs[i][mark_extr],
                           y22[mark_extr], 'xr', markersize=15)
            extension_bool = np.isin(all_frac_bool[i], 1)
            shear_bool = np.isin(all_frac_bool[i], 2)
            ulti_bool = np.isin(all_frac_bool[i], 3)
            ax3000[4].plot(all_diffs[i][extension_bool],
                           y22[extension_bool], 'xr', markersize=15)
            ax3000[4].plot(all_diffs[i][shear_bool],
                           y22[shear_bool], 'xg', markersize=15)
            ax3000[4].plot(all_diffs[i][ulti_bool],
                           y22[ulti_bool], 'xy', markersize=15)

        ax3000[0].set_xticks(ticks=ticks)
        plt.setp(ax3000[0], xticklabels=np.round(self.temperatures[::3], 1))
        ax3000[0].set_xlabel("Temperature [°C]")
        ax3000[0].set_ylabel("Vol% of fluid-filled porosity")

        ax3000[1].set_xlabel("Number of extractions")
        ax3000[1].set_ylabel("Tensile strength rock [MPa]")
        ax3000[1].set_title(f"{len(all_system_vol_pre)} rocks are plotted")
        tw1.set_xlabel("Max. fluid-filled porosity [%]")

        ax3000[2].set_xlabel("Temperature [°C]")
        ax3000[2].set_ylabel("Diff. stress [MPa]")
        ax3000[2].set_title(f"{len(all_system_vol_pre)} rocks are plotted")

        ax3000[3].set_xlabel("Diff. stress [MPa]")
        ax3000[3].set_ylabel("Fluid-filled porosity [Vol%]")
        ax3000[3].set_title(f"{len(all_system_vol_pre)} rocks are plotted")

        ax3000[4].set_xlabel("Diff. stress [MPa]")
        ax3000[4].set_ylabel("Tensile strength rock [MPa]")
        ax3000[4].set_title(f"{len(all_system_vol_pre)} rocks are plotted")

        plt.savefig(Path(main_folder /
                         "porosity_plot.png", dpi=400))
        plt.clf()
        plt.close()

    def compare_extraction_2(self, line, all_system_vol_pre, all_system_vol_post, all_fluid_before, all_tensile, all_diffs, all_frac_bool, permea, t_flux):

        colpal = sns.color_palette("flare", int(max(all_tensile)+1))
        fig, ax3000 = plt.subplots(
            1, 5, figsize=(20, 8), constrained_layout=True)
        tw12 = ax3000[1].twinx()
        tw13 = ax3000[2].twinx()

        for i, item in enumerate(all_system_vol_pre):
            # load system volume arrays
            system_vol_pre = np.array(all_system_vol_pre[i])
            system_vol_post = np.array(all_system_vol_post[i])
            st_fluid_before = all_fluid_before[i]

            # fluid filled porosity array
            y2 = (st_fluid_before)/system_vol_pre*100

            # set ticks
            ticks = np.arange(1, len(self.temperatures)+1, 1)[::3]

            # select extraction steps
            mark_extr = np.array(system_vol_pre-system_vol_post, dtype='bool')

            if np.size(t_flux[i]) == 1 and np.sum(t_flux[i]) == 0:
                pass
            else:
                if len(mark_extr) < len(self.temperatures):
                    selected_temperature = self.temperatures[arr_line[i]-1]
                    if len(permea[i]) > len(np.round(selected_temperature)[mark_extr]):
                        permea[i] = permea[i][permea[i] != 0]
                        t_flux[i] = t_flux[i][t_flux[i] != 0]
                else:
                    if len(permea[i]) > len(np.round(self.temperatures)[mark_extr]):
                        permea[i] = permea[i][permea[i] != 0]
                        t_flux[i] = t_flux[i][t_flux[i] != 0]
                # porsity line plot
                ax3000[0].plot(line[i], y2, '1--',
                               color=colpal[int(all_tensile[i])])
                ax3000[0].plot(line[i][mark_extr], y2[mark_extr], 'Dr')
                ax3000[0].annotate(str(all_tensile[i]), [line[-1], y2[-1]])
                # 2
                # diff stress over prograde with time int fluid flux
                if len(all_diffs[i]) < len(self.temperatures):
                    selected_temperature = self.temperatures[arr_line[i]-1]
                    ax3000[1].plot(np.round(selected_temperature), all_diffs[i],
                                   'x--', color=colpal[int(all_tensile[i])])
                else:
                    ax3000[1].plot(np.round(self.temperatures), all_diffs[i],
                                   'x--', color=colpal[int(all_tensile[i])])

                if len(mark_extr) < len(self.temperatures):
                    tw12.plot(np.round(selected_temperature)[
                        mark_extr], t_flux[i], 'sb')
                else:
                    tw12.plot(np.round(self.temperatures)[
                        mark_extr], t_flux[i], 'sb')
                # 3
                # diff stress over prograde with permeabilities
                if len(all_diffs[i]) < len(self.temperatures):
                    ax3000[2].plot(np.round(selected_temperature), all_diffs[i],
                                   'x--', color=colpal[int(all_tensile[i])])
                    tw13.plot(np.round(selected_temperature)[
                        mark_extr], permea[i], 'sb')
                else:
                    ax3000[2].plot(np.round(self.temperatures), all_diffs[i],
                                   'x--', color=colpal[int(all_tensile[i])])
                    tw13.plot(np.round(self.temperatures)[
                        mark_extr], permea[i], 'sb')
                # 4
                ax3000[3].plot(y2[mark_extr], permea[i], 's')

                # 5
                ax3000[4].plot(y2[mark_extr], t_flux[i], 's')

        ax3000[0].set_xticks(ticks=ticks)
        plt.setp(ax3000[0], xticklabels=np.round(self.temperatures[::3], 1))
        ax3000[0].set_xlabel("Temperature [°C]")
        ax3000[0].set_ylabel("Vol% of fluid-filled porosity")

        ax3000[1].set_xlabel("Temperature [°C]")
        ax3000[1].set_ylabel("Diff. stress [MPa]")
        tw12.set_ylabel("Time-int flux [$m^{3}/m^{2}$]")
        tw12.ticklabel_format(axis='y', style='sci', scilimits=(4, 4))
        tw12.set_yscale('log')
        ax3000[1].set_title(f"{len(all_system_vol_pre)} rocks are plotted")

        ax3000[2].set_xlabel("Temperature [°C]")
        ax3000[2].set_ylabel("Diff. stress [MPa]")
        tw13.set_yscale('log')
        ax3000[2].set_title(f"{len(all_system_vol_pre)} rocks are plotted")

        ax3000[3].set_ylabel("Permeability [$m^{2}$]")
        ax3000[3].set_yscale('log')
        ax3000[3].set_xlabel("Fluid-filled porosity [Vol%]")
        ax3000[3].set_title(f"{len(all_system_vol_pre)} rocks are plotted")

        ax3000[4].set_ylabel("Time-int flux [$m^{3}/m^{2}$]")
        ax3000[4].set_xlabel("Fluid-filled porosity [Vol%]")
        ax3000[4].set_title(f"{len(all_system_vol_pre)} rocks are plotted")

        plt.savefig(Path(main_folder/"quantify_fluid_fluxes.png", dpi=400))
        plt.clf()
        plt.close()

    def eval_permeability(self, all_system_vol_pre, all_system_vol_post, man_ing_permea, man_ing_depth, permea, t_flux, virtual_flux, virtual_permea):

        virtual_permea[-1][np.isnan(virtual_permea[-1])] = 0
        virtual_flux[-1][np.isnan(virtual_flux[-1])] = 0
        system_vol_pre = np.array(all_system_vol_pre)
        system_vol_post = np.array(all_system_vol_post)
        mark_extr = np.array(system_vol_pre-system_vol_post, dtype='bool')

        if len(self.line) < len(vtime[1]):
            selected_vtime = vtime[1][self.line-1]
            ttime = np.diff(selected_vtime[mark_extr], prepend=0)
        else:
            ttime = np.diff(vtime[1][mark_extr], prepend=0)

        test_p = mark_extr
        test_p = test_p.astype(float)
        # select True values for flux array (convert boolean to 1 or 0)
        test_flux = mark_extr
        test_flux = test_flux.astype(float)
        if np.size(permea) == 1 and np.sum(permea) == 0:
            pass
        else:
            if len(t_flux) > len(ttime):
                time_av_flux = t_flux[t_flux != 0.]/ttime
                test_p[test_p == 1] = permea[permea != 0]
                ltag = True
            else:
                time_av_flux = t_flux/ttime
                test_p[test_p == 1] = permea
                ltag = False
            # write values on the flux array
            test_flux[test_flux == 1] = time_av_flux

        z = np.arange(1, 110, 1)
        c_k = 10**(-14 - 3.2 * np.log10(z))
        c_k2 = 10**(-11.5 - 3.2 * np.log10(z))

        fig, axs = plt.subplots(1, 2, constrained_layout=True)

        if 'water.fluid' in self.phases:

            if line is not False and len(line) < len(all_depth[-1]):
                selected_all_depth = all_depth[-1][self.line-1]
                axs[0].plot(virtual_perm[-1], selected_all_depth /
                            1000, '--d', markersize=2)
            else:
                axs[0].plot(virtual_perm[-1], all_depth[-1] /
                            1000, '--d', markersize=2)

        if np.size(permea) == 1 and np.sum(permea) == 0:
            pass
        else:
            if ltag == True:
                up_err = permea[permea != 0]*10 - permea[permea != 0]
                down_err = permea[permea != 0] - permea[permea != 0]/10
                axs[0].errorbar(permea[permea != 0], all_depth[-1]
                                [mark_extr]/1000, xerr=[down_err, up_err], fmt='o')
            else:
                up_err = permea*10 - permea
                down_err = permea - permea/10
                if line is not False and len(line) < len(all_depth[-1]):
                    selected_all_depth = all_depth[-1][self.line-1]
                    axs[0].errorbar(permea, selected_all_depth[mark_extr] /
                                    1000, xerr=[down_err, up_err], fmt='o')
                else:
                    axs[0].errorbar(permea, all_depth[-1][mark_extr] /
                                    1000, xerr=[down_err, up_err], fmt='o')
            axs[0].plot(c_k, z, '-', color='#b18e4e', linewidth=3)
            axs[0].plot(c_k2, z, '-', color='#949494', linewidth=3)
            axs[0].set_ylim(100, 0)
            axs[0].set_xlim(1e-24, 1e-14)
            axs[0].set_xscale('log')
            axs[0].set_ylabel("Depth [km]", fontsize=16)
            axs[0].set_xlabel("Permeability [$m^{2}$]", fontsize=16)

        axs[0].plot(man_ing_permea, man_ing_depth, 'P',
                    color='#956b6b', markersize=3, alpha=0.5)

        axs[0].tick_params(axis='both', labelsize=12)
        tw11 = axs[1].twinx()

        if len(test_p) < len(self.temperatures):
            selected_temperature = self.temperatures[line-1]
            axs[1].plot(selected_temperature, test_p, '-x')
            tw11.plot(selected_temperature, test_flux, '-s', color='red')
        else:
            axs[1].plot(self.temperatures, test_p, '-x')
            tw11.plot(self.temperatures, test_flux, '-s', color='red')

        axs[0].fill_between([1e-24, 1e-21], [100, 100],
                            color="#c5c5c5", edgecolor='black', hatch="/")
        axs[1].set_ylabel("Permeability [$m^{2}$]")
        tw11.set_yscale('log')
        tw11.set_ylim(bottom=1e-4, top=1e1)
        tw11.set_ylabel("Time averaged flux [$m^{3}/m^{2}/year$]")
        axs[1].set_yscale('log')
        axs[1].set_xlabel("Temperature [°C]")
        axs[1].set_title(f"T={tensile}", loc='center')
        axs[0].annotate("log k = -14 - 3.2 log z", (3e-18, 19))
        axs[0].annotate("log k = -11.5 - 3.2 log z", (3e-17, 22))

        plt.savefig(Path(main_folder /
                         f"{self.rock_key}_permeability_of_slab.png", dpi=800))
        plt.clf()
        plt.close()


def calc_permeability(time_frame, physical_props, extraction_data, extraction_time):
    mü_water = 1e-4
    vtime = time_frame
    sys_physc = physical_props
    extrt = extraction_time
    ex_v = np.array(extraction_data.iloc[2, :])  # water volume extracted
    time_index = list(vtime[1])
    indi = []
    for time in extrt:
        indi.append(time_index.index(time))  # indizes for extraction timing
    # system volume at moment of extraction
    ex_vsys = sys_physc['system_vol_pre'][indi]
    dens_f = np.array(extraction_data.iloc[7, :])
    dens_r = sys_physc['system_density_post'][indi]
    dens_cont = (dens_r-dens_f)*1000
    dt_list = []
    prependindex = time_index.index(extrt[0])
    ttime_c = np.diff(extrt, prepend=time_index[prependindex-1])*365*24*60*60
    # time step passed between extraction
    for jj, item in enumerate(extrt):
        indi = time_index.index(item)
        cdt = time_index[indi]-time_index[indi-1]
        dt_list.append(cdt)
    ttime_c = np.array(dt_list)*365*24*60*60
    # for 1 cubic meter of rock - gives 1m2 as surface with column length z=1m => 1_000_000cm3 # no
    rlen = 1000_000_000
    # the time-integrated fluid flux for cm3/cm2
    tint_flux = ex_v * rlen/ex_vsys / 10_000
    vein = 0.1
    tint_flux_max = ex_v * rlen/(ex_vsys*vein) / 10_000
    tint_flux_min = ex_v*0.1 * rlen/ex_vsys / 10_000
    # calculate permeability after Bovay 2021 (after Ingebritsen and Manning 1999)
    # permeability is the time-int flux per duration (dt) and dependent of the viscosity of the fluid (mü)
    # and the driving-force gradient in the deep crust (diff of density rock to fluid multiplied with gravity)
    permea = tint_flux/100 * mü_water / ttime_c / 9.81 / dens_cont  # in m^2
    permea_max = tint_flux_max/100 * mü_water / ttime_c / 9.81 / dens_cont  # in m^2
    permea_min = tint_flux_min/100 * mü_water / ttime_c / 9.81 / dens_cont  # in m^2

    # pore fluid velocity
    # porosity of saturated porous medium as fluid volume per unit volume rock
    perc = ex_v/ex_vsys
    pore_fluid_v = tint_flux/perc/ttime_c*365*24*60*60  # cm/year
    return permea, tint_flux


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


def boxplot_to_GIF(
        file_name, main_folder, phase_data, group_key, legend_phases, color_set
):

    os.makedirs(f'{main_folder}/img_{file_name}', exist_ok=True)

    for i in range(len(phase_data.index)):
        mole_fractions = phase_data.iloc[i, :] / \
            np.sum(phase_data.iloc[i, :])*100
        os.makedirs(
            f'{main_folder}/img_{file_name}/phase_modes/{group_key}', exist_ok=True)
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
        plt.savefig(f'{main_folder}/img_{file_name}/phase_modes/{group_key}/img_{i}.png',
                    transparent=False, facecolor='white')
        # plt.pause(2)
        plt.close()

    # reading images and put it to GIF
    frames = []
    for i in range(len(phase_data.index)):
        image = imageio.v2.imread(
            f'{main_folder}/img_{file_name}/phase_modes/{group_key}/img_{i}.png')
        frames.append(image)
    imageio.mimsave(
        f'{main_folder}/img_{file_name}/phase_modes/{group_key}/phase_box.gif', frames, fps=2)


def permeability_to_GIF(
        file_name, main_folder, phase_data, group_key, filtered_permeability,
        depth, ts, ps, time
):

    os.makedirs(f'{main_folder}/img_{file_name}', exist_ok=True)

    for i in range(len(phase_data.index)):
        os.makedirs(
            f'{main_folder}/img_{file_name}/permeability/{group_key}', exist_ok=True)
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
        plt.savefig(f'{main_folder}/img_{file_name}/permeability/{group_key}/img_{i}.png',
                    transparent=False, facecolor='white')
        # plt.pause(2)
        plt.close()

    # reading images and put it to GIF
    frames = []
    for i in range(len(phase_data.index)):
        image = imageio.v2.imread(
            f'{main_folder}/img_{file_name}/permeability/{group_key}/img_{i}.png')
        frames.append(image)
    imageio.mimsave(
        f'{main_folder}/img_{file_name}/permeability/{group_key}/rock_permeabiltiy.gif', frames, fps=2)


def time_int_flux_to_GIF(
        file_name, main_folder, phase_data, group_key, filtered_int_flux, unfiltered_int_flux, regional_filtered_flux,
    regional_flux,
    depth, ts, ps
):

    os.makedirs(f'{main_folder}/img_{file_name}', exist_ok=True)

    for i in range(len(phase_data.index)):
        os.makedirs(
            f'{main_folder}/img_{file_name}/int_flux/{group_key}', exist_ok=True)
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
        ax4.vlines(10**4, 0, 100, linestyles='--', color='black', linewidth=4)

        props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

        # ax4.annotate("Regional range", (10**(2.7+0.5-0.3), 25), fontsize=16, bbox=props, rotation=90)
        ax4.annotate("Channelized", (10**(4+0.2), 5), fontsize=16, bbox=props)
        ax4.annotate("Pervasive", (10**(1+0.2), 5), fontsize=16, bbox=props)

    # Finishing plot - save
        plt.savefig(f'{main_folder}/img_{file_name}/int_flux/{group_key}/img_{i}.png',
                    transparent=False, facecolor='white')
        # plt.pause(2)
        plt.close()

    # reading images and put it to GIF
    frames = []
    for i in range(len(phase_data.index)):
        image = imageio.v2.imread(
            f'{main_folder}/img_{file_name}/int_flux/{group_key}/img_{i}.png')
        frames.append(image)
    imageio.mimsave(
        f'{main_folder}/img_{file_name}/int_flux/{group_key}/rock_int_flux.gif', frames, fps=2)


def porosity_to_GIF(
        file_name, main_folder, phase_data, group_key, filtered_porosity, unfiltered_porosity,
        depth, ts, ps
):

    os.makedirs(f'{main_folder}/img_{file_name}', exist_ok=True)

    for i in range(len(phase_data.index)):
        os.makedirs(
            f'{main_folder}/img_{file_name}/porosity/{group_key}', exist_ok=True)
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
        plt.savefig(f'{main_folder}/img_{file_name}/porosity/{group_key}/img_{i}.png',
                    transparent=False, facecolor='white')
        # plt.pause(2)
        plt.close()

    # reading images and put it to GIF
    frames = []
    for i in range(len(phase_data.index)):
        image = imageio.v2.imread(
            f'{main_folder}/img_{file_name}/porosity/{group_key}/img_{i}.png')
        frames.append(image)
    imageio.mimsave(
        f'{main_folder}/img_{file_name}/porosity/{group_key}/rock_porosity.gif', frames, fps=2)


def volume_ext_to_GIF(
        file_name, main_folder, phase_data, group_key, filtered_v_fluid_extr, unfiltered_v_fluid_extr,
        depth, ts, ps
):

    os.makedirs(f'{main_folder}/img_{file_name}', exist_ok=True)

    for i in range(len(phase_data.index)):
        os.makedirs(
            f'{main_folder}/img_{file_name}/volume/{group_key}', exist_ok=True)
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
        plt.savefig(f'{main_folder}/img_{file_name}/volume/{group_key}/img_{i}.png',
                    transparent=False, facecolor='white')
        # plt.pause(2)
        plt.close()

    # reading images and put it to GIF
    frames = []
    for i in range(len(phase_data.index)):
        image = imageio.v2.imread(
            f'{main_folder}/img_{file_name}/volume/{group_key}/img_{i}.png')
        frames.append(image)
    imageio.mimsave(
        f'{main_folder}/img_{file_name}/volume/{group_key}/extraction_volume.gif', frames, fps=2)


def pt_path_GIF(
        file_name, main_folder, phase_data, group_key, ts, ps
):

    os.makedirs(f'{main_folder}/img_{file_name}', exist_ok=True)

    for i in range(len(phase_data.index)):
        os.makedirs(
            f'{main_folder}/img_{file_name}/PTpath/{group_key}', exist_ok=True)
        # -------------------------------------------------------------------
        fig = plt.figure(constrained_layout=False,
                         facecolor='0.9', figsize=(9, 9))
        gs = fig.add_gridspec(nrows=3, ncols=3, left=0.15,
                              right=0.95, hspace=0.6, wspace=0.5)

        # plot 6
        ax4 = fig.add_subplot(gs[:, :])
        ax4.plot(ts[:i+1], ps[:i+1]/10_000, 'd--', color='black', markersize=4)
        ax4.plot(ts[i:i+1], ps[i:i+1]/10_000, 'd--',
                 color='#7fffd4', markersize=8, markeredgecolor='black')
        # Ax4 figure features
        ax4.set_xlim(200, 750)
        ax4.set_ylim(0, 3)
        ax4.set_ylabel("Pressure [GPa]", fontsize=18)
        ax4.set_xlabel("Temperature [°C]", fontsize=18)
        ax4.tick_params(axis='both', labelsize=18)

    # Finishing plot - save
        plt.savefig(f'{main_folder}/img_{file_name}/PTpath/{group_key}/img_{i}.png',
                    transparent=False, facecolor='white')
        # plt.pause(2)
        plt.close()

    # reading images and put it to GIF
    frames = []
    for i in range(len(phase_data.index)):
        image = imageio.v2.imread(
            f'{main_folder}/img_{file_name}/PTpath/{group_key}/img_{i}.png')
        frames.append(image)
    imageio.mimsave(
        f'{main_folder}/img_{file_name}/PTpath/{group_key}/PT_path.gif', frames, fps=2)


def gridplot_to_GIF(
        file_name, main_folder, phase_data, group_key, legend_phases, color_set,
        filtered_v_fluid_extr, unfiltered_v_fluid_extr, filtered_flux, unfiltered_flux, filtered_permeability,
        depth, ts, ps, time):

    # cheking the time step in the model
    time_steps = np.unique(np.diff(time[1]))
    if len(time_steps) > 1:
        print("Time alert")
        model_time_step = np.mean(time_steps)
    else:
        model_time_step = time_steps[0]

    # creating diroctory to store images
    os.makedirs(f'{main_folder}/img_{file_name}', exist_ok=True)
    tremor_flag = True
    if tremor_flag is True:
        tremor_data = pd.read_excel(
            r"C:\Users\Markmann\PhD\Data\03_Proj01_PetroModelling\tremor_data_pnsn_cascadia.xlsx")
        test_diff = int(tremor_data['depth'].max()-tremor_data['depth'].min())
        # plt.hist(tremor_data['depth'], bins=int(test_diff/2), orientation='horizontal')

    for i in range(len(phase_data.index)):
        mole_fractions = phase_data.iloc[i, :] / \
            np.sum(phase_data.iloc[i, :])*100
        os.makedirs(
            f'{main_folder}/img_{file_name}/{group_key}', exist_ok=True)
        # -------------------------------------------------------------------
        fig = plt.figure(constrained_layout=False,
                         facecolor='0.9', figsize=(9, 9))
        gs = fig.add_gridspec(nrows=3, ncols=3, left=0.10,
                              right=0.80, hspace=0.6, wspace=0.5)
        # All the single plots
        # plot 1
        y_offset = 0
        ax1 = fig.add_subplot(gs[2, 2])
        for t, val in enumerate(np.array(mole_fractions)):
            if val == np.float64(0):
                ax1.bar(
                    1, val, bottom=0, color=color_set[t], edgecolor='black', linewidth=0.1, label=legend_phases[t])
            else:
                ax1.bar(1, val, bottom=y_offset,
                        color=color_set[t], edgecolor='black', linewidth=0.1, label=legend_phases[t])
            y_offset = y_offset + val

        # plot 2
        ax2 = fig.add_subplot(gs[-1, -2])
        ax3 = fig.add_subplot(gs[-1, 0])
        ax4 = fig.add_subplot(gs[:-1, 0:2])
        ax3.plot(unfiltered_flux[:i+1], depth[:i+1]/1000, 'd--', color='green')
        ax2.plot(unfiltered_v_fluid_extr[:i+1],
                 depth[:i+1]/1000, 'd--', color='green')
        if False in filtered_permeability.mask[:i+1]:
            ax2.plot(filtered_v_fluid_extr[:i+1] /
                     1_000_000, depth[:i+1]/1000, 'd--')
            ax2.plot(filtered_v_fluid_extr[:i+1]/1_000_000,
                     depth[:i+1]/1000, 'd--', color='black', alpha=0.5)
            ax2.plot(filtered_v_fluid_extr[i:i+1]/1_000_000, depth[i:i+1]/1000,
                     'd--', color='#7fffd4', markersize=8, markeredgecolor='black')
            # plot 3
            ax3.plot(filtered_flux[:i+1], depth[:i+1] /
                     1000, 'd--', color='black', alpha=0.5)
            ax3.plot(filtered_flux[i:i+1], depth[i:i+1]/1000, 'd',
                     color='#7fffd4', markersize=8, markeredgecolor='black')
            # plot 4
            ax4.plot(
                filtered_permeability[:i+1], depth[:i+1]/1000, 'd--', color='black', alpha=0.7)
            ax4.plot(filtered_permeability[i:i+1], depth[i:i+1]/1000, 'd',
                     color='#7fffd4', markersize=8, markeredgecolor='black')
        else:
            # plot 4
            ax4.plot(np.ones(i+1)*1e-30, depth[:i+1]/1000, 'd')

        # plot 5
        ax5 = fig.add_subplot(gs[0:2, -1])
        ax5.hist(tremor_data['depth'], bins=int(test_diff/2),
                 orientation='horizontal', color='#bc3f67', ec='black')

        # plot 6
        ax6 = fig.add_subplot(gs[1:2, -1])
        ax6.plot(ts[:i+1], ps[:i+1]/10_000, 'd--', color='orange')

        # add Manning & Ingebritsen to plot 3
        z = np.arange(1, 110, 1)
        c_k = 10**(-14 - 3.2 * np.log10(z))
        ax4.plot(c_k, z, '--', color='#b18e4e', linewidth=3)

        # legend
        handles, labels = ax1.get_legend_handles_labels()
        legend_properties = {'weight': 'bold'}
        ax1.legend(handles[::-1], labels[::-1],
                   title='Phases', bbox_to_anchor=(1.2, 1.4))
        ax1.get_xaxis().set_visible(False)
        ax1.get_yaxis().set_visible(False)
        ax1.set_aspect(0.02)

        # Axs figure features
        ax4.set_title('T:{:.2f} °C P:{:.2f} GPa'.format(ts[i], ps[i]/10000))
        ax4.set_xlim(1e-26, 1e-14)
        ax4.set_xscale('log')
        ax4.set_ylim(100, 0)
        ax4.set_xlabel("permeability [$m^{2}$]", fontsize=18)
        ax4.set_ylabel("Depth [km]", fontsize=18)
        ax4.fill_between([1e-30, 1e-21], [100, 100],
                         color="#c5c5c5", edgecolor='black', hatch="/", alpha=0.4)
        ax4.tick_params(axis='both', labelsize=16)
        ax4.annotate("Manning &\nIngebritsen\n(1999)", (5e-17, 15))

        ax2.set_title("Extracted fluid volume")
        ax2.set_xlim(0, max(v_fluid_extr)/1_000_000)
        ax2.set_ylim(100, 0)
        # TODO ccm because it is not the geometry scale
        ax2.set_xlabel("Volume extracted [$m^{3}$]", fontsize=12)
        ax2.set_ylabel("Depth [km]", fontsize=12)
        ax2.tick_params(axis='both', labelsize=14)

        ax3.set_ylim(100, 0)
        if len(filtered_flux.compressed()) == 0:
            pass
        else:
            ax3.set_xscale('log')
        ax3.set_xlabel(
            "Time-int. fluid flux\n[$m^{3}$ $m^{-2}$ $yr^{-1}$]", fontsize=12)
        ax3.set_ylabel("Depth [km]", fontsize=12)
        # open/closed system line
        ax3.vlines(1e4/100/model_time_step, 0, 100, 'black')
        ax3.annotate("closed system", (1e-4+1e-4*2, 44),
                     (1e-4+1e-4*2, 44), rotation=90, fontsize=7)
        ax3.annotate("open system", (1e-4+1e-4*7, 40),
                     (1e-4+1e-4*7, 40), rotation=90, fontsize=7)
        ax3.tick_params(axis='both', labelsize=14)
        ax3.set_xlim(1e-8, 1e-2)

        ax5.set_ylim(100, 0)
        ax5.get_yaxis().set_visible(False)
        ax5.get_xaxis().set_visible(False)
        ax5.set_title("Abundance of\nseismic activity", fontsize=12)

        ax6.set_xlim(200, 750)
        ax6.set_ylim(0, 3)
        ax6.set_ylabel("Pressure [GPa]")
        ax6.set_xlabel("Temperature [°C]")

        # Finishing plot - save
        plt.savefig(f'{main_folder}/img_{file_name}/{group_key}/img_{i}.png',
                    transparent=False, facecolor='white')
        # plt.pause(2)
        plt.close()

    # reading images and put it to GIF
    frames = []
    for i in range(len(phase_data.index)):
        image = imageio.v2.imread(
            f'{main_folder}/img_{file_name}/{group_key}/img_{i}.png')
        frames.append(image)
    imageio.mimsave(
        f'{main_folder}/img_{file_name}/{group_key}/phase_box.gif', frames, fps=3)


#ANCHOR - Script
if __name__ == '__main__':
    # set the file
    o_file = file_opener()
    # set dictionary path as main folder
    main_folder = Path(o_file).parent.absolute()

    # porosity plot compilation storage
    all_system_vol_pre = []
    all_system_vol_post = []
    all_fluid_before = []
    all_tensile = []
    all_diffs = []
    all_frac_bool = []
    permea = []
    t_flux = []
    all_depth = []
    virtual_perm = []
    virtual_flux = []
    arr_line = []
    geometry = []

    man_ing_permea = np.array([-1.877540777917189629e+01, -1.868506900878293564e+01,
                               -1.800752823086574494e+01, -1.878670012547051371e+01,
                               -1.809786700125470560e+01, -1.819949811794228367e+01, -
                               1.830112923462986174e+01, -1.765746549560853040e+01,
                               -1.850439146800501788e+01, -1.721706398996235876e+01, -
                               1.752195734002509298e+01, -1.665244667503136711e+01,
                               -1.697992471769134326e+01, -1.699121706398996068e+01, -
                               1.658469259723964839e+01, -1.607653701380175448e+01,
                               -1.602007528230865674e+01, -1.598619824341279738e+01, -
                               1.567001254705144220e+01, -1.561355081555834445e+01,
                               -1.492471769134253634e+01, -1.434880803011292372e+01, -
                               1.309535759096612395e+01, -1.300501882057716507e+01,
                               -1.499247176913425506e+01, -1.500376411543287247e+01, -
                               1.420200752823086532e+01, -1.391969887076536949e+01,
                               -1.378419071518193206e+01])
    man_ing_permea = 1*10**man_ing_permea

    man_ing_depth = -np.array([-2.848404255319148604e+01, -2.538120567375886338e+01,
                               -2.358156028368794210e+01, -1.632092198581560183e+01,
                               -1.656914893617021534e+01, -1.464539007092198375e+01, -
                               1.284574468085106602e+01, -1.396276595744681259e+01,
                               -1.241134751773049771e+01, -1.098404255319148959e+01, -
                               7.260638297872340274e+00, -7.508865248226950229e+00,
                               -6.950354609929082272e+00, -4.654255319148937531e+00, -
                               3.102836879432626205e+00, -6.081560283687942103e+00,
                               -3.164893617021277805e+00, -4.964539007092199085e-01, -
                               2.792553191489361097e+00, -5.585106382978725748e+00,
                               -4.902482269503547485e+00, -2.792553191489361097e+00, -
                               2.668439716312057897e+00, -6.826241134751782624e-01,
                               -3.289007092198581006e+00, -2.171985815602837988e+00, -
                               1.054964539007087865e+00, -1.613475177304962926e+00,
                               -8.687943262411366163e-01])

    # open the thorPT output file
    with h5py.File(o_file, 'r') as f:
        rocks = list(f.keys())
        for group_key in rocks:
            print(f"Plotting in progress for {group_key}")
            # print("Keys: %s" % f[group_key].keys())
            # Get the chemical potential data
            pot_tag = list(f[group_key].attrs['pot_tag'])
            p_data = pd.DataFrame(f[group_key]['pot_data'])

            # P-T-t
            ts = np.array(f[group_key]['temperatures'])
            ps = np.array(f[group_key]['pressures'])
            vtime = np.array(f[group_key]['time_frame'])
            extr_time = np.array(f[group_key]['extr_time'])
            frac_bool = np.array(f[group_key]['fracture bool'])
            if 'line' in f[group_key].keys():
                line = np.array(f[group_key]['line'])
            else:
                line = False
            arr_line.append(line)

            tensile = np.float64(f[group_key]['tensile strength'])
            all_tensile.append(tensile)
            all_diffs.append(np.array(f[group_key]['diff. stress']))
            all_frac_bool.append(frac_bool)
            permea.append(np.array(f[group_key]['permeability']))
            t_flux.append(np.array(f[group_key]['time-int flux']))
            all_depth.append(np.array(f[group_key]['depth']))
            virtual_perm.append(np.array(f[group_key]['virtual permeability']))
            virtual_flux.append(
                np.array(f[group_key]['virtual time-int flux']))
            geometry.append(list(f[group_key]['geometry'].asstr()))

            # retrieve all the thermodyamica physical data
            df_var_d = {}
            for item in f[group_key]['df_var_dictionary'].keys():
                df_var_d[item] = pd.DataFrame(
                    f[group_key]['df_var_dictionary'][item])

            # fluid - extracted fluid data
            geometry
            extr_d = pd.DataFrame(f[group_key]['extracted_fluid_data'])
            fluid_before = np.array(f[group_key]['st_fluid_before'])
            # element data
            td_data = pd.DataFrame(f[group_key]['df_element_total'])

            # oxygen isotope data
            d_oxy = pd.DataFrame(f[group_key]['save_oxygen'])
            d_oxy_phases = list(f[group_key].attrs['oxy_data_phases'])
            d_oxy.columns = d_oxy_phases
            save_bulk_oxygen_pre = np.array(
                f[group_key]['save_bulk_oxygen_pre'])
            save_bulk_oxygen_post = np.array(
                f[group_key]['save_bulk_oxygen_post'])

            p_data.index = pot_tag
            vtime = vtime.T

            # phases
            phases = list(f[group_key].attrs['Phases'])
            phases2 = (f[group_key].attrs['Phases'])
            database = str(f[group_key].attrs['database'])
            # XMT naming and coloring
            phase_set, color_set = phases_and_colors_XMT(database, phases2)

            # get some physics arrays
            sys_physc = {}
            for item in f[group_key]['sys_physicals']:
                sys_physc[item] = np.array(f[group_key]['sys_physicals'][item])

            # store system volumes and fluid of all rocks
            all_system_vol_pre.append(sys_physc['system_vol_pre'])
            all_system_vol_post.append(sys_physc['system_vol_post'])
            all_fluid_before.append(fluid_before)

            # permea2, tint_flux = calc_permeability(
            #     vtime, sys_physc, extr_d, extr_time)
            # Call the plotting module and subroutines
            # ANCHOR - Plotting routines start
            plot_base = Plot_master(ts, ps, line, vtime[1], superior_data=[0, 1, 2, 3, 4, 5, 6], var_dic=df_var_d,
                                    td_data=td_data, ext_data=extr_d, physc_data=sys_physc,
                                    trackv=0, norm=0, hyd_data=0, rock_key=group_key, fluid_before=fluid_before, phases=phases, tensile=tensile)
            """# plot_base.chemp_plot(p_data)
            # plot_base.adv_stack_plot1_2(frac_bool, tag='df_vol%')
            plot_base.oxy_iso_plot(
                oxygen_data_frame=d_oxy, save_bulk_oxygen_pre=save_bulk_oxygen_pre, save_bulk_oxygen_post=save_bulk_oxygen_post)
            # plot_base.eval_permeability(all_system_vol_pre[-1], all_system_vol_post[-1], man_ing_permea, man_ing_depth, permea=permea[-1],
            #                             t_flux=t_flux[-1], virtual_flux=virtual_flux, virtual_permea=virtual_perm)
            # plot_base.eval_permeability(all_system_vol_pre[-1], all_system_vol_post[-1], man_ing_permea, man_ing_depth, permea=permea2,
            #                             t_flux=tint_flux, virtual_flux=virtual_flux, virtual_permea=virtual_perm)
            """
            # bar plot of phases
            phase_data = copy.deepcopy(df_var_d['df_N'])
            phase_data.columns = phases
            phase_data[np.isnan(phase_data) == True] = 0
            frames = []

            # create colormap
            colors = plt.cm.plasma(np.linspace(0, 1, len(phases)))

            # geometry factor
            bloc_a = np.float64(geometry[-1][0])
            bloc_b = np.float64(geometry[-1][1])
            bloc_c = np.float64(geometry[-1][2])
            area = bloc_b*bloc_c
            xxx = bloc_a
            size = bloc_a * bloc_b * bloc_c

            # fluid volume
            v_pre = sys_physc['system_vol_pre']
            v_post = sys_physc['system_vol_post']
            v_fluid_extr = v_pre - v_post
            scaling_factor = size*1_000_000/all_system_vol_pre[-1]
            v_fluid_extr = v_fluid_extr*scaling_factor

            # update 01.04.2023
            unfiltered_porosity = all_fluid_before[0]/v_pre
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
            solid_weight = mass_abs_data.T.sum()-mass_abs_data['water.fluid']
            solid_density = np.array(solid_weight)/np.array(solid_volumes)*1000

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

            # depth array
            depth = np.array(f[group_key]['depth'])

            # extraction marker
            # NOTE - new masking method
            extraction_boolean = np.array(all_frac_bool[-1], dtype=bool)
            extraction_boolean = np.invert(extraction_boolean)
            if False in extraction_boolean:
                filtered_flux = np.ma.masked_array(
                    t_flux[-1]*60*60*24*365, extraction_boolean)
                filtered_permeability = np.ma.masked_array(
                    permea[-1], extraction_boolean)
                filtered_v_fluid_extr = np.ma.masked_array(
                    v_fluid_extr, extraction_boolean)
                filtered_porosity = np.ma.masked_array(
                    unfiltered_porosity, extraction_boolean)
                filtered_int_flux = np.ma.masked_array(
                    unfiltered_int_flux, extraction_boolean)
                filtered_q_ti = np.ma.masked_array(
                    q_ti, extraction_boolean)
            else:
                # fake masked array because no fluid was extracted
                filtered_flux = np.ma.masked_array(
                    np.zeros(len(extraction_boolean)), extraction_boolean)
                filtered_permeability = np.ma.masked_array(
                    np.zeros(len(extraction_boolean)), extraction_boolean)
                filtered_v_fluid_extr = np.ma.masked_array(
                    np.zeros(len(v_fluid_extr)), extraction_boolean)
                filtered_v_fluid_extr = np.ma.masked_array(
                    np.zeros(len(extraction_boolean)), extraction_boolean)
                filtered_v_fluid_extr = np.ma.masked_array(
                    np.zeros(len(extraction_boolean)), extraction_boolean)
            unfiltered_flux = virtual_flux[-1]*60*60*24*365

            file_name = o_file.split('/')[-1].split('_')[1]

            if len(np.unique(all_frac_bool[-1])) == 1 and np.unique(all_frac_bool[-1])[-1] == 0:
                filtered_flux = np.zeros(len(all_frac_bool[-1]))
                filtered_permeability = np.zeros(len(all_frac_bool[-1]))
                filtered_v_fluid_extr = np.zeros(len(all_frac_bool[-1]))
                unfiltered_flux = np.zeros(len(all_frac_bool[-1]))
                filtered_q_ti = np.zeros(len(all_frac_bool[-1]))
                filtered_int_flux = np.zeros(len(all_frac_bool[-1]))
            else:
                gif_all = False
                if gif_all is True:
                    # ANCHOR deactivated GIF
                    """gridplot_to_GIF(
                        file_name=file_name,
                        main_folder=main_folder,
                        phase_data=phase_data,
                        group_key=group_key,
                        legend_phases=phase_set,
                        color_set=color_set,
                        filtered_v_fluid_extr=filtered_v_fluid_extr,
                        unfiltered_v_fluid_extr=v_fluid_extr,
                        filtered_flux=filtered_flux,
                        unfiltered_flux=unfiltered_flux,
                        filtered_permeability=filtered_permeability,
                        depth=depth,
                        ts=ts,
                        ps=ps,
                        time=vtime)"""

                    boxplot_to_GIF(file_name, main_folder, phase_data,
                                   group_key, legend_phases=phase_set, color_set=color_set)

                    """permeability_to_GIF(
                        file_name=file_name,
                        main_folder=main_folder,
                        phase_data=phase_data,
                        group_key=group_key,
                        filtered_permeability=filtered_permeability,
                        depth=depth,
                        ts=ts,
                        ps=ps,
                        time=vtime)"""
                    """porosity_to_GIF(
                        file_name=file_name,
                        main_folder=main_folder,
                        phase_data=phase_data,
                        group_key=group_key,
                        filtered_porosity=filtered_porosity,
                        unfiltered_porosity=unfiltered_porosity,
                        depth=depth,
                        ts=ts,
                        ps=ps)"""

                    time_int_flux_to_GIF(
                        file_name=file_name,
                        main_folder=main_folder,
                        phase_data=phase_data,
                        group_key=group_key,
                        filtered_int_flux=filtered_int_flux,
                        unfiltered_int_flux=unfiltered_int_flux,
                        regional_filtered_flux=filtered_q_ti,
                        regional_flux=q_ti,
                        depth=depth,
                        ts=ts,
                        ps=ps)

                    """volume_ext_to_GIF(
                        file_name=file_name,
                        main_folder=main_folder,
                        phase_data=phase_data,
                        group_key=group_key,
                        filtered_v_fluid_extr=filtered_v_fluid_extr,
                        unfiltered_v_fluid_extr=v_fluid_extr,
                        depth=depth,
                        ts=ts,
                        ps=ps)"""
            if group_key == 'rock0':
                pt_path_GIF(file_name, main_folder,
                            phase_data, group_key, ts, ps)

            print(f"Plotting finished for {group_key}")

        print("Closing plotting protocol...")

    line = np.arange(1, len(ts)+1, 1)

    # extraction versus differential stress
    diff_array = []
    extraction_array = []
    for item in all_diffs:
        diff_array.append(np.unique(item)[-1])

    for item in all_frac_bool:
        extraction_array.append(len(item[item != 0]))
    fig, ax = plt.subplots()
    ax.plot(diff_array, extraction_array, 'd')
    ax.set_yticks(np.arange(0, max(extraction_array)+1, 1))
    ax.tick_params(axis='both', labelsize=18)
    ax.set_ylabel("Number of extractions")
    ax.set_xlabel("Differential stress [MPa]")
    plt.show()

    if len(all_diffs[0]) == 0:
        pass
    else:
        plot_base.compare_extraction(arr_line, all_system_vol_pre, all_system_vol_post,
                                     all_fluid_before, all_tensile, all_diffs, all_frac_bool, filtered_permeability)
        # REVIEW scalar error - matplotlib gets infinite - no figure problem
        # plot_base.compare_extraction_2(arr_line, all_system_vol_pre, all_system_vol_post,
        #                                all_fluid_before, all_tensile, all_diffs, all_frac_bool, permea, t_flux)

    bulk_num = 1  # NOTE int(input("How many different bulk rocks?"))
    tensile_len = len(np.unique(all_tensile))

    if bulk_num == 3 and tensile_len == 5:
        """
        fig, ax = plt.subplots(constrained_layout=True)
        z = np.arange(0, 110, 1)
        c_k = 10**(-14 - 3.2 * np.log10(z))
        c_k2 = 10**(-11.5 - 3.2 * np.log10(z))
        ax.plot(c_k, z, '--', color='#b18e4e', linewidth=3)
        ax.plot(c_k2, z, '--', color='#949494', linewidth=3, alpha = 0.6)
        ax.plot(man_ing_permea, man_ing_depth, 'P',
                color='#956b6b', markersize=3, alpha=0.5)
        ax.set_ylim(100, 0)
        ax.set_xlim(1e-24, 1e-14)
        ax.set_xscale('log')
        ax.set_ylabel("Depth [$km$]", fontsize=18)
        ax.set_xlabel("Permeability [$m^{2}$]", fontsize=18)
        ax.tick_params(axis='both', labelsize=18)
        ax.fill_between([1e-24, 1e-21], [100, 100],
                        color="#c5c5c5", edgecolor='black', hatch="/")
        ax.annotate("Manning & Ingebritsen\n(1999)", (3e-18, 20))
        ax.annotate("Ingebritsen & Manning\n(2010)", (3e-17, 45))
        """
        coli = 5
        colpal = sns.color_palette("flare", coli)
        fig, ax = plt.subplots(constrained_layout=True)
        for i, item in enumerate(all_system_vol_pre):
            system_vol_pre = np.array(all_system_vol_pre[i])
            system_vol_post = np.array(all_system_vol_post[i])
            mark_extr = np.array(system_vol_pre-system_vol_post, dtype='bool')
            if np.size(permea[i]) == 1 and np.sum(permea[i]) == 0:
                pass
            else:
                if i < coli:
                    ax.plot(permea[i], all_depth[i][mark_extr] / 1000,
                            'd', color=colpal[i], markersize=9)
                if i > coli and i < coli*2:
                    ax.plot(permea[i], all_depth[i][mark_extr]/1000,
                            's', color=colpal[i-coli], markersize=9)
                if i > coli*2 and i < coli*3:
                    ax.plot(permea[i], all_depth[i][mark_extr]/1000,
                            '^', color=colpal[i-coli*2], markersize=9)

        z = np.arange(0, 110, 1)
        c_k = 10**(-14 - 3.2 * np.log10(z))
        c_k2 = 10**(-11.5 - 3.2 * np.log10(z))
        ax.plot(c_k, z, '--', color='#b18e4e', linewidth=3)
        ax.plot(c_k2, z, '--', color='#949494', linewidth=3, alpha=0.6)
        ax.plot(man_ing_permea, man_ing_depth, 'P',
                color='#956b6b', markersize=3, alpha=0.5)
        ax.set_ylim(100, 0)
        ax.set_xlim(1e-24, 1e-14)
        ax.set_xscale('log')
        ax.set_ylabel("Depth [$km$]", fontsize=18)
        ax.set_xlabel("Permeability [$m^{2}$]", fontsize=18)
        ax.tick_params(axis='both', labelsize=18)
        ax.fill_between([1e-24, 1e-21], [100, 100],
                        color="#c5c5c5", edgecolor='black', hatch="/")
        # legend handles
        [sym, ], [sym1, ], [sym2, ] = [plt.plot(z, "d", markersize=5, color='black'),
                                       plt.plot(z, "s", markersize=5,
                                                color='black'),
                                       plt.plot(z, "^", markersize=5,
                                                color='black'),
                                       ]
        [red_dot, ], [red_dot1, ], [red_dot2, ], [red_dot3, ], [red_dot4, ] = [
            plt.plot(z, "o", markersize=5, color=colpal[0]),
            plt.plot(z, "o", markersize=5, color=colpal[1]),
            plt.plot(z, "o", markersize=5, color=colpal[2]),
            plt.plot(z, "o", markersize=5, color=colpal[3]),
            plt.plot(z, "o", markersize=5, color=colpal[4])]
        # create legends
        leg1 = [red_dot, red_dot1, red_dot2, red_dot3, red_dot4]
        leg2 = [sym, sym1, sym2]
        first_legend = ax.legend(
            leg1, np.unique(all_tensile), loc='upper left', title="tensile strength", fontsize=12)
        plt.gca().add_artist(first_legend)
        ax.legend(leg2, ["bulk1", "bulk2", "bulk3"],
                  loc='lower left', fontsize=12)
        # ax.annotate("log k = -14 - 3.2 log z\n(Manning & Ingebritsen 1999)", (6.5e-19, 29))
        ax.annotate("Manning & Ingebritsen\n(1999)", (3e-18, 20))
        ax.annotate("Ingebritsen & Manning\n(2010)", (3e-17, 45))

        plt.savefig(Path(main_folder / "slab permeability.png", dpi=1200))
        plt.clf()
        plt.close()

        # working for continous extraction
        """
        coli = 5
    colpal = sns.color_palette("flare", coli)
    fig, ax = plt.subplots(constrained_layout=True)
    for i, item in enumerate(all_system_vol_pre):
        system_vol_pre = np.array(all_system_vol_pre[i])
        system_vol_post = np.array(all_system_vol_post[i])
        mark_extr = np.array(system_vol_pre-system_vol_post, dtype='bool')
        if np.size(permea[i]) == 1 and np.sum(permea[i]) == 0:
            pass
        else:
            if i < coli:
                ax.plot(permea[i][1:], all_depth[i][mark_extr] /1000,
                'd', color=colpal[i], markersize=9)
            if i > coli and i < coli*2:
                ax.plot(permea[i][1:], all_depth[i][mark_extr]/1000,
                        's', color=colpal[i-coli], markersize=9)
            if i > coli*2 and i < coli*3:
                ax.plot(permea[i], all_depth[i][mark_extr]/1000,
                        '^', color=colpal[i-coli*2], markersize=9)
    z = np.arange(0, 110, 1)
    c_k = 10**(-14 - 3.2 * np.log10(z))
    c_k2 = 10**(-11.5 - 3.2 * np.log10(z))
    ax.plot(c_k, z, '--', color='#b18e4e', linewidth=3)
    ax.plot(c_k2, z, '--', color='#949494', linewidth=3, alpha = 0.6)
    ax.plot(man_ing_permea, man_ing_depth, 'P',
            color='#956b6b', markersize=3, alpha=0.5)
    ax.set_ylim(100, 0)
    ax.set_xlim(1e-24, 1e-14)
    ax.set_xscale('log')
    ax.set_ylabel("Depth [$km$]", fontsize=18)
    ax.set_xlabel("Permeability [$m^{2}$]", fontsize=18)
    ax.tick_params(axis='both', labelsize=18)
    ax.fill_between([1e-24, 1e-21], [100, 100],
                    color="#c5c5c5", edgecolor='black', hatch="/")
    # legend handles
    [sym, ], [sym1, ], [sym2, ] = [plt.plot(z, "d", markersize=5, color='black'),
                                plt.plot(z, "s", markersize=5,
                                            color='black'),
                                plt.plot(z, "^", markersize=5,
                                            color='black'),
                                ]
    [red_dot, ], [red_dot1, ], [red_dot2, ], [red_dot3, ], [red_dot4, ] = [
        plt.plot(z, "o", markersize=5, color=colpal[0]),
        plt.plot(z, "o", markersize=5, color=colpal[1]),
        plt.plot(z, "o", markersize=5, color=colpal[2]),
        plt.plot(z, "o", markersize=5, color=colpal[3]),
        plt.plot(z, "o", markersize=5, color=colpal[4])]
    # create legends
    leg1 = [red_dot, red_dot1, red_dot2, red_dot3, red_dot4]
    leg2 = [sym, sym1, sym2]
    first_legend = ax.legend(
        leg1, np.unique(all_tensile), loc='upper left', title="tensile strength", fontsize=12)
    plt.gca().add_artist(first_legend)
    ax.legend(leg2, ["bulk1", "bulk2", "bulk3"],
            loc='lower left', fontsize=12)
    # ax.annotate("log k = -14 - 3.2 log z\n(Manning & Ingebritsen 1999)", (6.5e-19, 29))
    ax.annotate("Manning & Ingebritsen\n(1999)", (3e-18, 20))
    ax.annotate("Ingebritsen & Manning\n(2010)", (3e-17, 45))
    plt.savefig(Path(main_folder / "slab permeability.png", dpi=1200))
    plt.clf()
    plt.close()
        """

    if bulk_num == 5 and tensile_len == 6:
        coli = 6
        colpal = sns.color_palette("flare", coli)
        fig, ax = plt.subplots(constrained_layout=True)
        for i, item in enumerate(all_system_vol_pre):
            system_vol_pre = np.array(all_system_vol_pre[i])
            system_vol_post = np.array(all_system_vol_post[i])
            mark_extr = np.array(system_vol_pre-system_vol_post, dtype='bool')
            if np.size(permea[i]) == 1 and np.sum(permea[i]) == 0:
                pass
            else:
                if i < coli:
                    ax.plot(permea[i], all_depth[i][mark_extr] /
                            1000, 'd', color=colpal[i], markersize=9)
                if i > coli and i < coli*2:
                    ax.plot(permea[i], all_depth[i][mark_extr]/1000,
                            's', color=colpal[i-coli], markersize=9)
                if i > coli*2 and i < coli*3:
                    ax.plot(permea[i], all_depth[i][mark_extr]/1000,
                            '^', color=colpal[i-coli*2], markersize=9)
                if i > coli*3 and i < coli*4:
                    ax.plot(permea[i], all_depth[i][mark_extr]/1000,
                            'o', color=colpal[i-coli*3], markersize=9)
                if i > coli*4:
                    ax.plot(permea[i], all_depth[i][mark_extr]/1000,
                            'P', color=colpal[i-coli*4], markersize=9)
        z = np.arange(-10, 110, 1)
        c_k = 10**(-14 - 3.2 * np.log10(z))
        ax.plot(c_k, z, '-', color='#b18e4e', linewidth=3)
        ax.plot(man_ing_permea, man_ing_depth, 'P',
                color='#956b6b', markersize=3, alpha=0.5)
        ax.set_ylim(100, 0)
        ax.set_xlim(1e-24, 1e-14)
        ax.set_xscale('log')
        ax.set_ylabel("Depth $[km]$", fontsize=18)
        ax.set_xlabel("Permeability [$m^{2}$]", fontsize=18)
        ax.tick_params(axis='both', labelsize=18)
        ax.fill_between([1e-24, 1e-21], [100, 100],
                        color="#c5c5c5", edgecolor='black', hatch="/")
        # legend handles
        # plot euqal to number of bulk rocks
        [sym, ], [sym1, ], [sym2, ], [sym3, ], [sym4, ] = [plt.plot(z, "d", markersize=5, color='black'),
                                                           plt.plot(
                                                               z, "s", markersize=5, color='black'),
                                                           plt.plot(
                                                               z, "^", markersize=5, color='black'),
                                                           plt.plot(
                                                               z, "o", markersize=5, color='black'),
                                                           plt.plot(z, "P", markersize=5, color='black')]
        # symbols euqal to number of tensile strength values
        [red_dot, ], [red_dot1, ], [red_dot2, ], [red_dot3, ], [red_dot4, ], [red_dot5, ] = [
            plt.plot(z, "o", markersize=5, color=colpal[0]),
            plt.plot(z, "o", markersize=5, color=colpal[1]),
            plt.plot(z, "o", markersize=5, color=colpal[2]),
            plt.plot(z, "o", markersize=5, color=colpal[3]),
            plt.plot(z, "o", markersize=5, color=colpal[4]),
            plt.plot(z, "o", markersize=5, color=colpal[5])]
        # create legends
        # tensile strength
        leg1 = [red_dot, red_dot1, red_dot2, red_dot3, red_dot4, red_dot5]
        # number of bulk rocks
        leg2 = [sym, sym1, sym2, sym3, sym4]
        first_legend = ax.legend(leg1, np.unique(
            all_tensile), loc='upper left', title="tensile strength", handletextpad=0.1, fontsize=12)
        plt.gca().add_artist(first_legend)
        ax.legend(leg2, ["bulk1", "bulk2", "bulk3", "bulk4", "bulk5"],
                  loc='lower left', handletextpad=0.1, fontsize=12)
        ax.annotate("log k = -14 - 3.2 log z", (3e-18, 19))

        plt.savefig(Path(main_folder / "slab permeability.png", dpi=1200))
        plt.clf()
        plt.close()

    if bulk_num == 4 and tensile_len == 5:
        coli = 5
        colpal = sns.color_palette("flare", coli)
        fig, ax = plt.subplots(constrained_layout=True)
        for i, item in enumerate(all_system_vol_pre):
            system_vol_pre = np.array(all_system_vol_pre[i])
            system_vol_post = np.array(all_system_vol_post[i])
            mark_extr = np.array(system_vol_pre-system_vol_post, dtype='bool')
            if np.size(permea[i]) == 1 and np.sum(permea[i]) == 0:
                pass
            else:
                if i < coli:
                    ax.plot(permea[i], all_depth[i][mark_extr] /
                            1000, 'd', color=colpal[i], markersize=9)
                if i > coli and i < coli*2:
                    ax.plot(permea[i], all_depth[i][mark_extr]/1000,
                            's', color=colpal[i-coli], markersize=9)
                if i > coli*2 and i < coli*3:
                    ax.plot(permea[i], all_depth[i][mark_extr]/1000,
                            '^', color=colpal[i-coli*2], markersize=9)
                if i > coli*3 and i < coli*4:
                    ax.plot(permea[i], all_depth[i][mark_extr]/1000,
                            'o', color=colpal[i-coli*3], markersize=9)
                if i > coli*4:
                    ax.plot(permea[i], all_depth[i][mark_extr]/1000,
                            'P', color=colpal[i-coli*4], markersize=9)
        z = np.arange(-10, 110, 1)
        c_k = 10**(-14 - 3.2 * np.log10(z))
        ax.plot(c_k, z, '-', color='#b18e4e', linewidth=3)
        ax.plot(man_ing_permea, man_ing_depth, 'P',
                color='#956b6b', markersize=3, alpha=0.5)
        ax.set_ylim(100, 0)
        ax.set_xlim(1e-24, 1e-14)
        ax.set_xscale('log')
        ax.set_ylabel("Depth $[km]$", fontsize=18)
        ax.set_xlabel("Permeability [$m^{2}$]", fontsize=18)
        ax.tick_params(axis='both', labelsize=18)
        ax.fill_between([1e-24, 1e-21], [100, 100],
                        color="#c5c5c5", edgecolor='black', hatch="/")
        # legend handles
        # plot euqal to number of bulk rocks
        [sym, ], [sym1, ], [sym2, ], [sym3, ] = [plt.plot(z, "d", markersize=5, color='black'),
                                                 plt.plot(
                                                     z, "s", markersize=5, color='black'),
                                                 plt.plot(
                                                     z, "^", markersize=5, color='black'),
                                                 plt.plot(z, "o", markersize=5, color='black')]
        # symbols euqal to number of tensile strength values
        [red_dot, ], [red_dot1, ], [red_dot2, ], [red_dot3, ], [red_dot4, ] = [
            plt.plot(z, "o", markersize=5, color=colpal[0]),
            plt.plot(z, "o", markersize=5, color=colpal[1]),
            plt.plot(z, "o", markersize=5, color=colpal[2]),
            plt.plot(z, "o", markersize=5, color=colpal[3]),
            plt.plot(z, "o", markersize=5, color=colpal[4])]
        # create legends
        # tensile strength
        leg1 = [red_dot, red_dot1, red_dot2, red_dot3, red_dot4]
        # number of bulk rocks
        leg2 = [sym, sym1, sym2, sym3]
        first_legend = ax.legend(leg1, np.unique(
            all_tensile), loc='upper left', title="tensile strength", handletextpad=0.1, fontsize=12)
        plt.gca().add_artist(first_legend)
        ax.legend(leg2, ["bulk1", "bulk2", "bulk3", "bulk4"],
                  loc='lower left', handletextpad=0.1, fontsize=12)
        ax.annotate("log k = -14 - 3.2 log z", (3e-18, 19))

        plt.savefig(Path(main_folder / "slab permeability.png", dpi=1200))
        plt.clf()
        plt.close()
