
"""
Written by
Thorsten Markmann
thorsten.markmann@unibe.ch
status: 16.07.2024
"""

import shutil
import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import h5py
import copy

from pathlib import Path

# from thorpt_thorPT.valhalla.tunorrad import *
# from thorpt_thorPT.valhalla.Pathfinder import *
# from valhalla.tunorrad import *
# from valhalla.Pathfinder import *
from .tunorrad import *
from .Pathfinder import *
from dataclasses import dataclass
from tqdm import tqdm



def set_origin():
    """
    Sets the current working directory to the directory of the current file.
    """
    dirname = os.path.dirname(os.path.abspath(__file__))
    os.chdir(dirname)


def file_save_path():
    """
    Prompts the user to select a file path to save the data.

    Returns:
        str: The selected file path.
    """
    # write the data to a file
    # first get the filename
    validFile = False
    while not validFile:
        fileout = filedialog.asksaveasfilename(
            title="Select file to save the data",
            filetypes=[("Simple text files (.hdf5)", "*.hdf5")],
            defaultextension='hdf5'
        )
        if len(fileout) == 0:
            # nothing selected, pop up message to retry
            messagebox.showinfo(
                "Filename error", "Please select a filename to save the data.")
        else:
            validFile = True
    return fileout


def whole_rock_convert_3(ready_mol_bulk=0):
    """
    Converts the whole rock composition to a format passable to theriak.

    Args:
        ready_mol_bulk (pd.DataFrame or float, optional): The whole rock composition. If a DataFrame is provided,
            it should have a single column representing the mole fractions of each element. Defaults to 0.

    Returns:
        str: The converted whole rock composition in a format passable to theriak.
    """
    # Rock input, water input and more chemical starting conditions!!!
    # //////////// whole rock composition selection passable to theriak ////////////

    # excess oxygen test
    #         SIO2 TIO2 AL2O3 FEO FE2O3 MNO MGO CAO NA2O K2O H2O CO2
    oxy_frame = {'SI': 2, 'TI': 2, 'AL': 3/2, 'FE': 1, 'MN': 1,
                 'MG': 1, 'CA': 1, 'NA': 1/2, 'K': 1/2, 'H': 1/2, 'C': 2}

    if isinstance(ready_mol_bulk, pd.DataFrame):
        if len(ready_mol_bulk.columns) > 1:
            ready_mol_bulk = ready_mol_bulk.iloc[:, 0]
    test = ready_mol_bulk * \
        pd.DataFrame(oxy_frame, oxy_frame.keys()).iloc[0, :]
    oxy_diff = ready_mol_bulk['O'] - test.sum()
    if oxy_diff < 0.01 and oxy_diff > 0.0:
        simple_iron = test.sum()
    else:
        simple_iron = ready_mol_bulk['O']
    # simple_iron = np.round(simple_iron, 7)
    # plt.plot(temperature, oxy_diff, '^r')
    # print(f"{ready_mol_bulk['O'] } - {test.sum()} = {oxy_diff}")

    bulk = ready_mol_bulk
    new_bulk1 = []
    scan_element = []
    scan_val = []
    for el in bulk.index:
        if el == 'E':
            pass
        elif el == 'O':
            scan_element.append(el)
            scan_val.append(simple_iron)
            new_bulk1.append(el+'('+str(simple_iron)+')')
        else:
            scan_element.append(el)
            val = bulk.loc[el]
            if val < 0.0 or val < 1e-06:
                val = 0
            scan_val.append(np.round(val, 6))
            new_bulk1.append(scan_element[-1]+'('+str(scan_val[-1])+')')
    new_bulk = ''.join(new_bulk1) + "    " + "*theriak-out preprocessed bulk"

    return new_bulk


def whole_rock_to_weight_normalizer(rock_bulk=[32.36, 0.4, 8.78, 2.91, 0.0, 0.0, 1.45, 23.16, 1.96, 1.66],
                                    ready_bulk=False, init_water=0.2, init_carbon=0.41):
    """
    Converts the whole rock composition to weight-normalized bulk composition.

    Args:
        rock_bulk (list): List of the bulk composition of the rock. Default is [32.36, 0.4, 8.78, 2.91, 0.0, 0.0, 1.45, 23.16, 1.96, 1.66].
        ready_bulk (bool): Flag indicating if the bulk composition is ready. Default is False.
        init_water (float): Initial water content. Default is 0.2.
        init_carbon (float): Initial carbon content. Default is 0.41.

    Returns:
        tuple: A tuple containing the weight-normalized bulk composition and the total oxygen content.

    """

    # Rock input, water input and more chemical starting conditions!!!
    # //////////// whole rock composition selection passable to theriak ////////////
    # Oxides and stochimetry
    #                Si Ti Al Fe2 Fe3 Mn Mg Ca Na K + H C O
    num_cation_ext = [1, 1, 2, 1,  2,  1, 1, 1, 2, 2,  2, 1]

    #             SIO2   TIO2 AL2O3 FEO  FE2O3 MNO  MGO   CAO    NA2O  K2O
    num_cation = [1, 1, 2, 1,  2,  1, 1, 1, 2, 2]
    #         SIO2 TIO2 AL2O3 FEO
    #              FE2O3 MNO MGO CAO NA2O K2O H2O CO2
    num_oxy = [2,  2,   3,    1,  3,    1,  1,  1,  1,   1,  1,  2]
    oxide_mass = [60.09, 79.866, 101.96, 71.85,
                  159.68, 70.94, 40.3, 56.08, 61.98, 94.2]
    # SI        TI       AL       FE2      FE3      MN
    #  MG      CA       NA       K        H        C        O
    element_mass = np.array([28.08550, 47.88000, 26.98154, 55.84700, 55.84700,
                            54.93085, 24.30500, 40.07800, 22.98977, 39.09830, 1.00794, 12.01100, 15.99940])
    # Si Al Fe2 Fe3
    nbcat = [1, 2, 1, 1, 1, 1, 2, 1, 2, 2, 1]
    nbO = [2, 3, 1, 1, 1, 1, 1, 2, 1, 1, 2]
    # rock
    # water and carbon input
    init_water = init_water
    init_carbon = init_carbon

    # Testing for Fe speciations and adapting excess oxygen
    nO_FeO = rock_bulk[3]/oxide_mass[3]*num_oxy[3]
    nO_Fe2O3 = rock_bulk[4]/oxide_mass[4]*num_oxy[4]
    # num of oxygen in FeO-tot comes needs wt% calculated from amounts of FeO and Fe2O3
    # amount of Fe2O3 needs weighted conversion of 0.899 (from Molar Mass and cation ratio)
    wt_FeO_tot = rock_bulk[3] + (oxide_mass[3] /
                                 oxide_mass[4]*num_cation_ext[4]) * rock_bulk[4]
    nO_FeO_tot = wt_FeO_tot / oxide_mass[3]*num_oxy[3]
    # excess oxygen = FeO-oxy + Fe2O3-oxy - FeO tot-oxy
    # is zero if no Fe2O3 is in bulk rock
    excess_oxy = nO_FeO + nO_Fe2O3 - nO_FeO_tot

    # Rock bulk to mol
    rock_bulk_mol_orig = np.array(rock_bulk)/np.array(oxide_mass)*num_cation
    mol_frame = np.concatenate(
        (rock_bulk_mol_orig, init_water, init_carbon), axis=None)
    # Recalculate oxygen
    normal_oxy = excess_oxy + np.sum(mol_frame/num_cation_ext*num_oxy)
    mol_frame = np.concatenate((mol_frame, normal_oxy), axis=None)

    if ready_bulk == 'Marco':
        bulk = np.around(mol_frame, 5)
        check = [bulk[0], bulk[2], bulk[3], bulk[5], bulk[6],
                 bulk[7], bulk[8], bulk[1], bulk[9], bulk[10], bulk[11]]
        new_bulk = (
            f"SI({bulk[0]})AL({bulk[2]})FE({bulk[3]})"
            f"MN({bulk[5]})MG({bulk[6]})CA({bulk[7]})"
            f"NA({bulk[8]})TI({bulk[1]})K({bulk[9]})"
            f"H({100})C({bulk[11]})O(?)O({excess_oxy})    * MarcoBulk"
        )
        rockOxy = sum(np.array(check)/np.array(nbcat)*np.array(nbO))

    elif ready_bulk is False:
        # element masses for kg/mol
        element_mass = element_mass/1e3
        # Normed rock bulk mol
        norm_mass = np.sum(mol_frame*element_mass)
        # - normalize rock mass to 1kg, conversion of the moles
        mol_frame_norm = mol_frame*1/norm_mass
        excess_oxy = excess_oxy*1/norm_mass
        # mol_frame_norm = mol_frame*0.01/norm_mass
        # excess_oxy = excess_oxy*0.01/norm_mass

        # Update bulk
        bulk = mol_frame_norm
        # NOTE rounding is enabled
        bulk = np.around(mol_frame_norm, 5)
        check = [bulk[0], bulk[2], bulk[3], bulk[5], bulk[6],
                 bulk[7], bulk[8], bulk[1], bulk[9], bulk[10], bulk[11]]
        # REVIEW "Calculate the proper Fe amount"
        new_bulk = (
            f"SI({bulk[0]})AL({bulk[2]})FE({bulk[3]})"
            f"MN({bulk[5]})MG({bulk[6]})CA({bulk[7]})"
            f"NA({bulk[8]})TI({bulk[1]})K({bulk[9]})"
            f"H({bulk[10]})C({bulk[11]})O(?)O({excess_oxy})    * CalculatedBulk"
        )

        rockOxy = sum(np.array(check)/np.array(nbcat)*np.array(nbO))

    return new_bulk, rockOxy


def oxygen_isotope_recalculation(isotope_data, oxygen_data):
    """
    Recalculates the oxygen isotope composition based on isotope data and oxygen data.

    Args:
        isotope_data (list): List of dictionaries containing isotope data.
        oxygen_data (pandas.DataFrame): DataFrame containing oxygen data.

    Returns:
        float: The recalculated oxygen isotope composition.

    """
    temp_oxygen_data = pd.DataFrame(
        isotope_data[-1]['delta_O'],
        index=isotope_data[-1]['Phases']
    )
    # Oxygen signature from last fractionation
    last_oxygen = isotope_data[-1]
    # read moles of oxygen
    collect_phases = []
    for phase in last_oxygen['Phases']:
        collect_phases.append(oxygen_data.loc['O'][phase])
    phase_oxygen = np.array(collect_phases)
    phase_doxy = np.array(temp_oxygen_data.iloc[:, 0])
    # test for nan in oxygen data - nan is kicked out and phases in phase oxygen are adapted equally
    phase_oxygen = phase_oxygen[np.logical_not(np.isnan(phase_doxy))]
    phase_doxy = phase_doxy[np.logical_not(np.isnan(phase_doxy))]
    new_O_bulk = sum(phase_oxygen*phase_doxy / sum(phase_oxygen))
    return new_O_bulk


def fluid_injection_isotope_recalculation(
        isotope_data, oxygen_data, input_deltaO,
        input_hydrogen, input_oxygen,
        interaction_factor=1, fluid_name_tag='water.fluid'):
    """
    Recalculates the bulk oxygen isotope signature after fluid injection.

    Args:
        isotope_data (list): List of dictionaries containing isotope data for each fractionation step.
        oxygen_data (pandas.DataFrame): DataFrame containing oxygen data for different phases.
        input_deltaO (pandas.DataFrame): DataFrame containing input delta O values for different phases.
        input_hydrogen (float): Input hydrogen value.
        input_oxygen (float): Input oxygen value.
        interaction_factor (float, optional): Interaction factor. Defaults to 1.

    Returns:
        float: The new bulk oxygen isotope signature.
    """

    # last oxygen isotope data saved to temporary variable
    temp_oxygen_data = pd.DataFrame(
        isotope_data[-1]['delta_O'],
        index=isotope_data[-1]['Phases']
            )

    # Oxygen signature from last fractionation
    last_oxygen = isotope_data[-1]

    # read moles of oxygen
    collect_phases = []
    for phase in last_oxygen['Phases']:
        collect_phases.append(oxygen_data.loc['O'][phase])
    phase_oxygen = np.array(collect_phases)
    phase_doxy = np.array(temp_oxygen_data.iloc[:, 0])

    # test for nan in oxygen data - nan is kicked out and phases in phase oxygen are adapted equally
    phase_oxygen = phase_oxygen[np.logical_not(np.isnan(phase_doxy))]
    phase_doxy = phase_doxy[np.logical_not(np.isnan(phase_doxy))]

    # read oxygen isotope data from input
    temp_input = pd.DataFrame(
        input_deltaO['delta_O'],
        index=input_deltaO['Phases']
            )

    # influx data to reclaculate with
    input_deltaO = np.float64(temp_input.loc[fluid_name_tag, 0])
    print(f"Fluid influx is d18O = {input_deltaO}")
    input_oxygen = interaction_factor * input_oxygen

    # recalculation by factors
    bulk_deltaO = sum(phase_oxygen*phase_doxy / sum(phase_oxygen))
    new_O_bulk = (input_oxygen*input_deltaO + sum(phase_oxygen)*bulk_deltaO ) / (input_oxygen+sum(phase_oxygen))
    # print("New bulk oxygen isotope signature is: ", new_O_bulk)

    if fluid_name_tag in oxygen_data.columns:
        #fluid incoming and fluid present in infiltrated rock
        fluid_mix_d18O = (input_oxygen*input_deltaO + oxygen_data.loc['O'][fluid_name_tag]*temp_oxygen_data.loc[fluid_name_tag][0]) / (input_oxygen+oxygen_data.loc['O'][fluid_name_tag])
    else:
        #fluid incoming but no fluid present in infiltrated rock
        fluid_mix_d18O = input_deltaO

    print(f"Fluid mix d18O is: {fluid_mix_d18O}")

    #changing the fluid isotope composition after interaction with the rock
    if fluid_name_tag in oxygen_data.columns:
        #fluid incoming and fluid present in infiltrated rock
        fluid_end_d18O = fluid_mix_d18O - (sum(phase_oxygen) - oxygen_data.loc['O'][fluid_name_tag]) / (input_oxygen+oxygen_data.loc['O'][fluid_name_tag]) * (new_O_bulk-bulk_deltaO)
        free_fluid_oxygen = (input_oxygen+oxygen_data.loc['O'][fluid_name_tag]) # the oxygen of all the free fluid at the moment of interaction
    else:
        #fluid incoming but no fluid present in infiltrated rock
        fluid_end_d18O = fluid_mix_d18O - (sum(phase_oxygen)) / (input_oxygen) * (new_O_bulk-bulk_deltaO)
        free_fluid_oxygen = (input_oxygen) # the oxygen of all the free fluid at the moment of interaction

    print(f"Fluid after interaction d18O is: {fluid_end_d18O}")
    

    # test with append
    phase_oxygen = np.append(phase_oxygen, input_oxygen)
    phase_doxy = np.append(phase_doxy, input_deltaO)
    bulk_deltaO = sum(phase_oxygen*phase_doxy / sum(phase_oxygen))
    # print("New bulk oxygen isotope signature is: ", bulk_deltaO)

    return new_O_bulk, fluid_end_d18O, free_fluid_oxygen

# Progressbar init
def progress(percent=0, width=40):
    """
    Display a progress bar with a given percentage and width.

    Parameters:
    percent (float): The percentage of progress (default is 0).
    width (int): The width of the progress bar (default is 40).

    Returns:
    None
    """
    left = width * percent // 100
    right = width - left
    tags = "#" * int(left)
    spaces = " " * int(right)
    percents = f"{percent:.1f}%"
    print("\r[", tags, spaces, "]", percents, sep="", end="", flush=True)


def read_trace_element_content():
    # get script file path
    # Get the script file location
    script_path = os.path.dirname(os.path.abspath(__file__))
    # Navigate one folder back
    # parent_folder = os.path.dirname(script_path)
    # Enter DataFiles folder
    new_folder_path = os.path.join(script_path, 'DataFiles')

    # read txt file trace_element_budget.txt
    trace_element_bulk = pd.read_csv(
        os.path.join(new_folder_path, 'trace_element_budget.txt'),
        sep=', ', header=0, index_col=0, engine='python'
    )
    return trace_element_bulk

def modify_trace_element_content(trace_element_bulk, trace_element_distirbution, min_name):
    """
    Modifies the trace element content based on the given modification.

    Args:
        trace_element_bulk (pandas.DataFrame): The trace element content.
        trace_element_distirbution (dict): The trace element distribution.
        min_name (str): The mineral name.
    Returns:
        pandas.DataFrame: The modified trace element content.
    """

    # subtract the trace element content from the bulk from the mineral defined by min_name
    # first search the name of min_name in the trace_element_bulk
    if min_name in trace_element_distirbution.index:
        trace_element_bulk = trace_element_bulk - np.array(trace_element_distirbution.loc[min_name])
    
    return trace_element_bulk

@dataclass
class rockactivity:
    """
    A class representing rock activity.

    Attributes:
        function: The function attribute.
        react: The react attribute.
    """
    function: any
    react: any


class ThorPT_Routines():

    """
    Module for different modelling routines.
    Modelling sequence is following:
    1) Initialize rock and execute minimization"
    2) Data frame storage"
    3) Updating storage"
    4) Testing for free fluid & saving system+fluid volume properties"
    5) Oxygen isotopes module"
    6) Mineral fractionation module"
    7) Metastable garnet"
    8) Mechanical failure model"
    """

    def __init__(self,
            temperatures, pressures, master_rock, rock_origin,
            track_time, track_depth, garnet_fractionation, path_methods,
            lowest_permeability, speed, angle, time_step, theriak, debugging_recorder
            ):
        """
        Initialize the ThorPT class.

        Parameters:
        temperatures (list): List of temperatures.
        pressures (list): List of pressures.
        master_rock (dict): Dictionary containing rock information.
        rock_origin (str): Origin of the rock.
        track_time (bool): Flag indicating whether to track time.
        track_depth (bool): Flag indicating whether to track depth.
        garnet_fractionation (bool): Flag indicating whether to perform garnet fractionation.
        path_methods (list): List of path methods.
        lowest_permeability (float): Lowest permeability value.
        speed (float): Speed value.
        angle (float): Angle value.
        time_step (float): Time step value.
        theriak (str): Theriak value.
        """
        # Output variables
        self.temperatures = temperatures
        self.pressures = pressures
        self.rock_dic = master_rock
        self.track_time = track_time
        self.track_depth = track_depth
        self.rock_origin = rock_origin
        self.garnet_fractionation = garnet_fractionation
        # self.mechanical_methods = mechanical_methods
        self.path_methods = path_methods
        self.minimum_permeability = lowest_permeability
        self.speed = speed
        self.angle = angle
        self.time_step = time_step
        self.theriak = theriak
        self.debugging_recorder = debugging_recorder

        # trace element bulk by default
        self.trace_element_bulk = read_trace_element_content()

    def unreactive_multi_rock(self):
        """
        Perform calculations for unreactive multi-rock system.

        This method performs calculations for an unreactive multi-rock system. It initializes thermodynamic conditions,
        calculates thermo data using the theriak wrapper, stores and merges the calculated data, and performs additional
        calculations such as MicaPotassium, SystemFluidTest, oxygen-isotope module, and mineral fractionation.

        Returns:
            None
        """

        # Main variables petrology
        temperatures = self.temperatures
        pressures = self.pressures
        master_rock = self.rock_dic
        rock_origin = self.rock_origin
        track_time = self.track_time
        track_depth = self.track_depth

        # Main variables mechanical model
        lowest_permeability = self.minimum_permeability
        # NOTE shear_stress = self.shear_stress

        """
        # REVIEW - mechanical methods
        # Methods
        factor_method = self.mechanical_methods[0]
        steady_method = self.mechanical_methods[1]
        dynamic_method = self.mechanical_methods[2]
        coulomb = self.mechanical_methods[3]
        coulomb_permea = self.mechanical_methods[4]
        coulomb_permea2 = self.mechanical_methods[5]
        """
        # initialize the trace element composition of the bulk rock
        for item in master_rock.keys():
            # trace element distribution
            master_rock[item]['init_trace_element_bulk'] = self.trace_element_bulk

        # Main variables fractionation
        grt_frac = self.garnet_fractionation

        v_fluid_cubic_track = []
        v_fluid_cubicp_track = []

        # //////////////////////////////////////////////////
        # /////////////////////////////////////////////////
        # ////////////////////////////////////////////////
        count = 0
        """k = 0
        kk = len(temperatures)*len(master_rock)
        progress(int(k/kk)*100)"""
        print("Script: unreactive_multi_rock")
        for num, temperature in enumerate(tqdm(temperatures, desc="Processing modelling steps")):

            # print('\n')
            # print("New calculation iteration")
            # print("===================")
            # print(f"==== 1) time = {track_time[num]} years,\n==== 2) depth = {track_depth[num]}.")

            # //////////////////////////////////////////////////////////////////////////
            # preparing bulk rock for calculation - normalized to 1kg of rock or it is Marco
            for item in master_rock:
                if num < 1:
                    master_rock[item]['new_bulk'], rockOxy = whole_rock_to_weight_normalizer(
                        rock_bulk=master_rock[item]['bulk'][:-2],
                        init_water=float(master_rock[item]['bulk'][-2]),
                        init_carbon=float(master_rock[item]['bulk'][-1])
                        )
                    # master_rock[item]['new_bulk'] = 'AL(0.0325)MG(1.212)FE(0.045)SI(0.6598)H(10)O(?)    * PierreBulk'
                else:
                    master_rock[item]['new_bulk'] = whole_rock_convert_3(
                        ready_mol_bulk=master_rock[item]['df_element_total']['total:']
                    )
                    # print("Bulk rock stop")

                # store the current used bulk rock to the backup dictionary
                rock_origin[item]['bulk'].append(master_rock[item]['new_bulk'])

                # print(f"{item} Bulk rock composition checked.")
            # print(f"All bulk rock compositions were checked. No error found")
            # print("\n")
            all_rock_keys_list = list(master_rock.keys())
            # LINK 1) Initialisation of the rock system
            # Initialize thermodynamic conditions called from P-T-t path - for each rock in the dictionary
            for jjj, item in enumerate(tqdm(master_rock, desc="Model calculation Gibbs energy min., Oxygen fractionation, Trace Elements:")):
                # Start of modelling scheme for rock
                print("\n")
                # print("¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦")
                # print("v v v v v v v v v v v v v v v v v v v v v v v v")
                # print("\N{Runic Letter Mannaz Man M}\N{Runic Letter Isaz Is Iss I}\N{Runic Letter Naudiz Nyd Naud N}\N{Runic Letter Isaz Is Iss I}\N{Runic Letter Mannaz Man M}\N{Runic Letter Isaz Is Iss I}\N{Runic Letter Algiz Eolhx}\N{Runic Letter Ansuz A}\N{Runic Letter Tiwaz Tir Tyr T}\N{Runic Letter Isaz Is Iss I}\N{Runic Letter Othalan Ethel O}\N{Runic Letter Naudiz Nyd Naud N}")
                # print(f"-> Forward modeling step initiated - {item}")
                # display modelling progress
                """ic = k/kk*100
                k += 1
                progress(ic)
                print("\n")"""

                # print(master_rock[item]['new_bulk'])
                """master_rock[item]['new_bulk'] = 'AL(0.0325)MG(1.212)FE(0.045)SI(0.6598)H(10)O(?)    * PierreBulk'"""

                # tracking theriak input before minimization
                # save temperature, pressure and master_rock[item]['new_bulk']
                # test if varibale is dictionary
                if isinstance(master_rock[item]['theriak_input_record'], dict) == False:
                    master_rock[item]['theriak_input_record'] = {}
                    master_rock[item]['theriak_input_record']['temperature'] = [temperature]
                    master_rock[item]['theriak_input_record']['pressure'] = [pressures[num]]
                    master_rock[item]['theriak_input_record']['bulk'] = [master_rock[item]['new_bulk']]
                    # print("Tracking theriak -> Create dictionary -> first entry")

                # test for empty dictionary
                elif isinstance(master_rock[item]['theriak_input_record'], dict) == True:
                    master_rock[item]['theriak_input_record']['temperature'].append(temperature)
                    master_rock[item]['theriak_input_record']['pressure'].append(pressures[num])
                    master_rock[item]['theriak_input_record']['bulk'].append(master_rock[item]['new_bulk'])
                    # print("Tracking theriak -> dictionary exists -> add entry")


                # Send information to the theriak wrapper
                """
                # test for metastable calculation
                if num < 1:
                    master_rock[item]['new_bulk'] = 'O(22.398)AL(1.877)CA(1.486)FE(0.596)H(2.993)K(0.032)MG(1.725)NA(0.66)SI(6.709)TI(0.206)    *metastable test'
                """

                master_rock[item]['minimization'] = Therm_dyn_ther_looper(self.theriak,
                        master_rock[item]['database'], master_rock[item]['new_bulk'],
                        temperature, pressures[num], master_rock[item]['df_var_dictionary'],
                        master_rock[item]['df_h2o_content_dic'], master_rock[item]['df_element_total'],
                        num, fluid_name_tag=master_rock[item]['database_fluid_name'])
                    

            # //////////////////////////////////////////////////////////////////////////

                # calculating difference between new Volumes and previous P-T-step volumes - for derivate, not difference!!!
                if num < 1:
                    master_rock[item]['master_norm'].append(np.sqrt(
                        (temperature-temperatures[num])**2 + (pressures[num]-pressures[num])**2))
                else:
                    master_rock[item]['master_norm'].append(np.sqrt(
                        (temperature-temperatures[num-1])**2 + (pressures[num]-pressures[num-1])**2))

            # //////////////////////////////////////////////////////////////////////////

            # Calculating and passing thermo data by theriak and theriak wrapper
            # ----> self.df_var_dictionary, self.df_all_elements, self.df_hydrous_data_dic, new_fluid_volumes
                # print("bulk rock minimization")
                if item != 'rock000':
                    master_rock[item]['minimization'].thermodynamic_looping_station(
                        theriak_input_rock_before = master_rock[all_rock_keys_list[jjj-1]]['theriak_input_record']
                        )
                else:
                    master_rock[item]['minimization'].thermodynamic_looping_station()

                # Main dictionary save
                # saving G_sys per mol of system
                master_rock[item]['g_sys'].append(master_rock[item]['minimization'].g_sys /
                                                sum(master_rock[item]['minimization'].df_phase_data.iloc[0, :]))
                # g_sys.append(minimization.g_sys)
                master_rock[item]['pot_data'].append(
                    master_rock[item]['minimization'].pot_frame)

                # Backup dictionary save
                rock_origin[item]['g_sys'] = copy.deepcopy(
                    master_rock[item]['g_sys'][-1])
                rock_origin[item]['pot_data'] = copy.deepcopy(
                    master_rock[item]['pot_data'][-1])

                # //////////////////////////////////////////////////////////////////////////

                # Creating DataFrame structure for "df_var_dictionary" in first itteration
                # LINK - 2) Setup dataframes
                if num < 1:
                    # Volume and Density ouput - Dataframes (df_N, df_Vol% etc)
                    for variable in list(master_rock[item]['minimization'].df_phase_data.index):
                        master_rock[item]['df_var_dictionary']['df_' + str(variable)] = pd.DataFrame()
                    water_cont_ind = ["N", "H2O[pfu]", "H2O[mol]",
                                    "H2O[g]", "wt%_phase", "wt%_solids", "wt%_H2O.solid"]
                    for variable in water_cont_ind:
                        master_rock[item]['df_h2o_content_dic']['df_' +
                                                                str(variable)] = pd.DataFrame()
                    # Copy to backup dictionary
                    rock_origin[item]['df_var_dictionary'] = copy.deepcopy(
                        master_rock[item]['df_var_dictionary'])
                    rock_origin[item]['df_h2o_content_dic'] = copy.deepcopy(
                        master_rock[item]['df_h2o_content_dic'])

                # updating dictionary with newly calculated data
                master_rock[item]['minimization'].merge_dataframe_dic()
                # print("\n")

                # print("////// Energy minimization executed //////")
                # print("\n")

                # //////////////////////////////////////////////////////////////////////////
                # multi-rock loop for updating data storage, MicaPotassium, SystemFluidTest, init oxygen-isotope module, mineral fractionation
                # //////////////////////////////////////////////////////////////////////////
                # LINK - 3) Data storage & merge
                # calling dictionaries and dataframe for up-to-date usage
                # print(f"Running data re-storage, MicaPotassium, SystemFluidTest, oxy-module and mineral fractionation")
                # for item in master_rock:
                master_rock[item]['df_var_dictionary'], master_rock[item]['df_h2o_content_dic'], master_rock[item]['df_element_total'] = (
                    master_rock[item]['minimization'].df_var_dictionary,
                    master_rock[item]['minimization'].df_hydrous_data_dic,
                    master_rock[item]['minimization'].df_all_elements)

                #print(master_rock[item]['df_var_dictionary']['df_N'])

                master_rock[item]['df_element_total'] = master_rock[item]['df_element_total'].iloc[:, :-1]
                # hydrogen content of the system before extraction
                if 'H' in master_rock[item]['df_element_total'].index:
                    master_rock[item]['total_hydrogen'] = master_rock[item]['df_element_total']['total:']['H']
                else:
                    master_rock[item]['total_hydrogen'] = 0
                master_rock[item]['st_elements'] = pd.concat(
                    [master_rock[item]['st_elements'], master_rock[item]['df_element_total']['total:']], axis=1)

                # Backup dictionary - merging the data
                for kkey in rock_origin[item]['df_var_dictionary'].keys():
                    cdata = pd.concat(
                        [rock_origin[item]['df_var_dictionary'][kkey], master_rock[item]['df_var_dictionary'][kkey].iloc[:, -1]], axis=1)
                    rock_origin[item]['df_var_dictionary'][kkey] = copy.deepcopy(
                        cdata)
                for kkey in rock_origin[item]['df_h2o_content_dic'].keys():
                    if len(master_rock[item]['df_h2o_content_dic'][kkey].index) > 0:
                        cdata = pd.concat(
                            [rock_origin[item]['df_h2o_content_dic'][kkey], master_rock[item]['df_h2o_content_dic'][kkey].iloc[:, -1]], axis=1)
                    else:
                        cdata = rock_origin[item]['df_h2o_content_dic'][kkey]
                    rock_origin[item]['df_h2o_content_dic'][kkey] = copy.deepcopy(
                        cdata)

                cdata = master_rock[item]['df_element_total']
                rock_origin[item]['df_element_total'].append(copy.deepcopy(cdata))

                # //////////////////////////////////////////////////////////////////////////
                # store Mica potassium if stable
                for phase in master_rock[item]['df_element_total'].columns:
                    if 'PHNG' in phase:
                        master_rock[item]['mica_K'].append(
                            [temperature, master_rock[item]['df_element_total'][phase]['K']])

                # //////////////////////////////////////////////////////////////////////////
                # LINK - 4) System and fluid volumes
                # Checking for fluid/solid volumes at t = 0 and t = -1,
                # calculating difference (used for escape/extraction rule e.g., factor method)
                # print("-> Testing for aq fluid in the system")
                master_rock[item]['minimization'].step_on_water()

                # Checking fluid and solid volumes. Storing t(-1) and t(0) data before calling calculation
                # //////////////
                # crucial step in saving the system volumes of solid and fluid phase to update and keep track
                if num < 1:
                    master_rock[item]['fluid_volume_new'] = master_rock[item]['minimization'].new_fluid_Vol
                    master_rock[item]['fluid_volume_before'] = master_rock[item]['fluid_volume_new']
                    master_rock[item]['solid_volume_new'] = master_rock[item]['minimization'].new_solid_Vol
                    master_rock[item]['solid_volume_before'] = master_rock[item]['solid_volume_new']
                else:
                    master_rock[item]['fluid_volume_before'] = master_rock[item]['minimization'].free_water_before
                    master_rock[item]['fluid_volume_new'] = master_rock[item]['minimization'].new_fluid_Vol

                    """if master_rock[item]['minimization'].solid_vol_before != master_rock[item]['st_solid'][-1]:
                        print("\nWARNING: solid volume mismatch\n")"""
                    master_rock[item]['solid_volume_before'] = master_rock[item]['minimization'].solid_vol_before
                    master_rock[item]['solid_volume_before'] = master_rock[item]['st_solid'][-1]

                    master_rock[item]['solid_volume_new'] = master_rock[item]['minimization'].new_solid_Vol


                # //////////////////////////////////////////////////////////////////////////
                # LINK - 5) Oxygen fractionation module
                # isotope fractionation module
                # print("-> Oxygen isotope module initiated")
                master_rock[item]['model_oxygen'] = Isotope_calc(
                    master_rock[item]['df_var_dictionary']['df_N'], master_rock[item]['minimization'].sol_sol_base,
                    master_rock[item]['df_element_total'], oxygen_signature=master_rock[item]['bulk_oxygen'])
                # print(f"Bulk oxygen is {master_rock[item]['bulk_oxygen']}")

                # Isotope calc - function for oxygen isotope signatures
                master_rock[item]['model_oxygen'].frac_oxygen(temperature)

                # storing isotope fractionation result, dic in list appended
                master_rock[item]['save_oxygen'].append(
                    master_rock[item]['model_oxygen'].oxygen_dic)
                master_rock[item]['save_bulk_oxygen_pre'].append(
                    master_rock[item]['bulk_oxygen'])

                ### Backup dictionary - save oxygen data
                rock_origin[item]['save_oxygen'].append(copy.deepcopy(master_rock[item]['model_oxygen'].oxygen_dic))

                # //////////////////////////////////////////////////////////////////////////
                # LINK - 5-2) Trace element module
                master_rock[item]['model_tracers'] = TraceElementDistribution(
                    master_rock[item]['df_var_dictionary']['df_N'],
                    master_rock[item]['minimization'].sol_sol_base,
                    master_rock[item]['df_element_total'],
                    master_rock[item]['init_trace_element_bulk'],
                    database=master_rock[item]['database'])
                # call the distribution of tracers
                trace_df = master_rock[item]['model_tracers'].distribute_tracers(temperature, pressure=pressures[num], iteration=num)
                # save the tracer data
                master_rock[item]['trace_element_data'][(num, pressures[num], temperature)] = trace_df
                master_rock[item]['trace_element_bulk'][(num, pressures[num], temperature)] = master_rock[item]['init_trace_element_bulk']

                # //////////////////////////////////////////////////////////////////////////
                # LINK - 6) Mineral Fractionation
                # mineral (garnet) fractionation - coupled oxygen bulk modification
                if grt_frac == True:
                    if master_rock[item]['database'] == 'ds62mp.txt' or master_rock[item]['database'] == 'td-ds62-mb50-v07.txt':
                        garnet_name = 'GRT'
                    else:
                        garnet_name = 'GARNET'

                    # print("-> Mineral fractionation initiated")
                    # old frac position this line
                    for phase in master_rock[item]['df_element_total'].columns:
                        if '_' in phase:
                            pos = phase.index('_')
                            name = phase[:pos]
                            if name == garnet_name:
                                # Michelles atigorite fractionation
                                # # if name=='GARNET' or name=='SERP' or name=='BR':
                                # modify the oxygen signature of the bulk rock
                                new_bulk_oxygen = master_rock[item]['minimization'].mineral_fractionation(
                                    master_rock[item]['save_oxygen'][-1], name)
                                

                                """
                                plt.figure()
                                norming = np.array([
                                    0.3670, 0.9570, 0.1370, 0.7110, 0.2310, 0.0870, 0.3060, 
                                    0.0580, 0.3810, 0.0851, 0.2490, 0.0356, 0.2480, 0.0381])
                                for i, value in enumerate(master_rock[item]['trace_element_bulk'].values()):
                                    test = value/norming
                                    plt.plot(test.T, '.-', label=i)
                                # plt.xlabel(test.columns)
                                plt.ylabel('Normalized values')
                                plt.yscale('log')
                                plt.legend()
                                """

                                # modify the trace element content
                                # last entry of the trace element data
                                #data_storage_keys = list(master_rock[item]['data_storage'].keys())
                                # Check if there are any entries
                                if np.size(trace_df) > 0:
                                    #last_entry_key = data_storage_keys[-1]
                                    #last_entry_value = master_rock[item]['data_storage'][last_entry_key]

                                    master_rock[item]['init_trace_element_bulk'] = modify_trace_element_content(
                                        trace_element_bulk=master_rock[item]['init_trace_element_bulk'],
                                        trace_element_distirbution=trace_df,
                                        min_name = phase
                                    )

                                master_rock[item]['garnet'].append(
                                    master_rock[item]['minimization'].separate)

                                # Collect the volume of each formed garnet
                                master_rock[item]['meta_grt_volume'].append(master_rock[item]['minimization'].separate.volume)

                                # print(
                                #     f"Selected phase = {phase} with Vol% = {master_rock[item]['df_var_dictionary']['df_vol%'].loc[phase][-2:]}")
                                # print(
                                #     f"Bulk deltaO changed from {round(master_rock[item]['bulk_oxygen'], 3)} to {round(new_bulk_oxygen, 3)}")
                                # print temperature, pressure and new_bulk_oxygen
                                # print(f"Fractionation of {name} at {temperature} and {pressures[num]}. Bulk oxygen changed from {master_rock[item]['bulk_oxygen']} to {new_bulk_oxygen}")
                                master_rock[item]['bulk_oxygen'] = new_bulk_oxygen
                                print("_______________________")
                                master_rock[item]['garnet_check'].append(1)
                            if len(master_rock[item]['garnet_check']) < num:
                                master_rock[item]['garnet_check'].append(0)

                master_rock[item]['df_element_total'] = master_rock[item]['minimization'].df_all_elements
                # LINK - 7) Metastable garnet
                # Adding the metastable garnet impact
                # calculate the metastable garnet for all bits besides the last
                # add the calculated volume to the solid volume of the current step (this is then saved to the st_solid and used next turn)
            for jjj, item in enumerate(tqdm(master_rock, desc="Processing metastable garnet calculation:")):
                grt_flag = True
                if grt_flag is True and len(master_rock[item]['garnet']) > 0 and len(master_rock[item]['garnet_check']) > 1:
                    if master_rock[item]['garnet_check'][-1] == 0:
                        # take al modelled garnets
                        # LINK - Metastable garnet call
                        # print(f"1MStab-Grt {temperature} --- {pressures[num]}")
                        # print("garnet minimization")
                        metastable_garnet = Garnet_recalc(self.theriak, master_rock[item]['garnet'], temperature, pressures[num])
                        metastable_garnet.recalculation_of_garnets(database=master_rock[item]['database'], garnet_name=garnet_name)
                        print(f"Fluid volume = {master_rock[item]['fluid_volume_new']} ccm")
                        print(f"Solid volume = {master_rock[item]['solid_volume_new']} ccm")
                        volume = metastable_garnet.recalc_volume
                        # Adding the metastable volume of the fractionated garnet to the current minimized solied volume!
                        master_rock[item]['solid_volume_new'] += volume
                        master_rock[item]['meta_grt_weight'].append(metastable_garnet.recalc_weight)
                    if len(master_rock[item]['garnet']) > 1 and master_rock[item]['garnet_check'][-1] == 1:
                        # take all garnets but last one
                        # print(f"2MStab-Grt {temperature} --- {pressures[num]}")
                        # print("garnet minimization")
                        metastable_garnet = Garnet_recalc(self.theriak, master_rock[item]['garnet'][:-1], temperature, pressures[num])
                        metastable_garnet.recalculation_of_garnets(database=master_rock[item]['database'], garnet_name=garnet_name)
                        # print(f"Fluid volume = {master_rock[item]['fluid_volume_new']} ccm")
                        # print(f"Solid volume = {master_rock[item]['solid_volume_new']} ccm")
                        volume = metastable_garnet.recalc_volume
                        # Adding the metastable volume of the fractionated garnet to the current minimized solied volume!
                        master_rock[item]['solid_volume_new'] += volume
                        master_rock[item]['meta_grt_weight'].append(metastable_garnet.recalc_weight)
                elif grt_flag is True and len(master_rock[item]['garnet']) > 0:
                    # nothing to do when garnet is stable in first step
                    pass
                else:
                    # Adding zero volume or weight because no metastable garnet is present
                    # calculation only to make it clear
                    master_rock[item]['solid_volume_new'] += np.array(master_rock[item]['meta_grt_volume']).sum()
                    if np.array(master_rock[item]['meta_grt_volume']).sum() != np.float64(0):
                        print("What happend now? Here is suddenly garnet which is not stable???")
                        # keyboard.wait('esc')
                    metastable_garnet_weight = 0

                # keeping track of stored and removed fluid
                # !!!! necessary for mechanical module
                master_rock[item]['st_fluid_before'].append(master_rock[item]['fluid_volume_new'])
                master_rock[item]['st_fluid_after'].append(master_rock[item]['fluid_volume_new'])
                master_rock[item]['st_solid'].append(master_rock[item]['solid_volume_new'])

            print("All rocks passed the petrochemical model set-up, next step physical fracturing model...")

            # //////////////////////////////////////////////////////////////////////////
            # LINK - MECHANICAL FAILURE MODEL
            # Physical fracture model section - active if free fluid is present
            # Checking if free fluid is present. Stores this value, initializes
            # lowest_permeability = 1e-20
            # print(f"Testing systems with the failure model")

            for jjj, item in enumerate(master_rock):
                fluid_name_tag = master_rock[item]['database_fluid_name']
                if fluid_name_tag in list(master_rock[item]['df_element_total'].columns) and master_rock[item]['fluid_volume_new'] > 0:
                    print("-> Fluid extraction test")
                    master_rock[item]['fluid_hydrogen'] = master_rock[item]['df_element_total'][fluid_name_tag]['H']

                    # Prepare fluid extraction
                    # - last entry of st_fluid_after has last fluid volume (number when not extracted or zero when extracted)
                    if len(master_rock[item]['st_fluid_after']) > 1:
                        fluid_before = master_rock[item]['st_fluid_after'][-2]
                    else:
                        fluid_before = master_rock[item]['st_fluid_after'][-1]


                    # Start Extraction Master Module
                    master_rock[item]['fluid_calculation'] = Ext_method_master(
                                        pressures[num], temperature,
                                        master_rock[item]['df_var_dictionary']['df_volume/mol'].loc[fluid_name_tag].iloc[-1],
                                        fluid_before, master_rock[item]['fluid_volume_new'],
                                        master_rock[item]['solid_volume_before'], master_rock[item]['solid_volume_new'],
                                        master_rock[item]['save_factor'], master_rock[item]['master_norm'][-1],
                                        master_rock[item]['minimization'].df_phase_data,
                                        master_rock[item]['tensile strength'],
                                        differential_stress= master_rock[item]['diff. stress'],
                                        friction= master_rock[item]['friction'],
                                        fluid_pressure_mode= master_rock[item]['fluid_pressure_mode'],
                                        fluid_name_tag=fluid_name_tag ,subduction_angle=self.angle,
                                        extraction_treshold = master_rock[item]['extraction treshold'],
                                        extraction_connectivity = master_rock[item]['fluid connectivity']
                                        )
                    # //////////////////////////////////////////////////////////////////////////
                    # ////// Calculation for new whole rock /////////
                    # !!!!Assumption: fracturing of the system
                    # if condition = open system, free water gets extracted
                    fracturing_flag = False # Trigger by default False - active when coulomb module becomes positive
                    failure_mech = master_rock[item]['Extraction scheme']
                    # SECTION selection for fluid release
                    if failure_mech in ['Factor', 'Dynamic', 'Steady', 'Mohr-Coulomb-Permea2', 'Mohr-Coulomb-Griffith'] and num > 0:
                        # if factor_method is True or dynamic_method is True or steady_method is True or coulomb is True or coulomb_permea is True or coulomb_permea2 is True:

                            # virtual momentary fluid flux and permeabiltiy test
                            mü_water = 1e-4

                            # water data
                            v_water = float(master_rock[item]['df_var_dictionary']['df_volume[ccm]'].loc[fluid_name_tag].iloc[-1])
                            d_water = float(master_rock[item]['df_var_dictionary']['df_density[g/ccm]'].loc[fluid_name_tag].iloc[-1])
                            weigth_water = master_rock[item]['df_var_dictionary']['df_wt[g]'].iloc[:, -1][fluid_name_tag]
                            # system data
                            master_rock[item]['meta_grt_weight']
                            v_system = master_rock[item]['solid_volume_new'] + master_rock[item]['fluid_volume_new'] # modelling solid phases + metastable garnet + fluid

                            if len(master_rock[item]['meta_grt_weight']) > 0:
                                weight_sys = float(
                                    master_rock[item]['df_var_dictionary']['df_wt[g]'].iloc[:, -1].sum()
                                    ) + master_rock[item]['meta_grt_weight'][-1] # weight of solid + fluid + metastable garnet
                            else:
                                weight_sys = float(
                                    master_rock[item]['df_var_dictionary']['df_wt[g]'].iloc[:, -1].sum()
                                    ) # weight of solid + fluid + metastable garnet

                            d_system = float(master_rock[item]['df_var_dictionary']['df_wt[g]'].iloc[:, -1].sum())/v_system
                            # solid rock data
                            v_rock = master_rock[item]['solid_volume_new']
                            weight_rock = weight_sys - weigth_water # weight of solids + metastable garnet
                            d_rock = weight_rock/v_rock
                            # density difference
                            density_cont = (d_rock-d_water)*1000 # kg/m3

                            # time interval
                            if num < 1:
                                tc_time_c = float(np.diff(track_time[num:num+2])*365*24*60*60)
                            else:
                                tc_time_c = float(np.diff(track_time[num-1:num+1])*365*24*60*60)

                            # test 4
                            """
                            cubic_meter = 0.6
                            xxx = np.power(1_000_000 * cubic_meter, 1/3) # length of x of cube (100cm for 1 m3)
                            size = xxx**3 # cubic size
                            area = xxx**2 # surface size
                            # fluid flux = drainiage flux throughout a column
                            volume_flux = v_water * size /v_system/area/tc_time_c # cm3 cm-2 s-1
                            volume_flux = volume_flux/100 # m3 m-2 s-1
                            # integrated permeability
                            int_permea = volume_flux*mü_water/9.81/xxx/density_cont # permeability in m2"""

                            # LINK Fluid flux and permeability
                            # test 05
                            bloc_a = np.float64(master_rock[item]['geometry'][0])
                            bloc_b = np.float64(master_rock[item]['geometry'][1])
                            bloc_c = np.float64(master_rock[item]['geometry'][2])
                            area = bloc_b*bloc_c
                            xxx = bloc_a
                            size = bloc_a * bloc_b * bloc_c
                            v_water1 = v_water/1000000 # cm3 to m3
                            v_system1 = v_system/1000000 # cm3 to m3

                            v_water*1_000_000/master_rock[item]['solid_volume_before']

                            """volume_flux = v_water1 * size/v_system1 /area/tc_time_c # m3 m-2 s-1
                            volume_flux = v_water1/tc_time_c /area * 1/v_system1  # m3 m-2 s-1
                            int_permea = volume_flux*mü_water/9.81/1/density_cont # permeability in m2"""
                            volume_flux = 0
                            volume_flux = 0
                            int_permea = 0

                            # unit bloc
                            """bloc_a = 1
                            bloc_b = 1
                            bloc_c = 1
                            area = bloc_b*bloc_c
                            xxx = bloc_a
                            size = bloc_a * bloc_b * bloc_c
                            v_water1 = v_water/1000000 # cm3 to m3
                            v_system1 = v_system/1000000 # cm3 to m3
                            volume_flux = v_water1 * size/v_system1 /area/tc_time_c # m3 m-2 s-1
                            int_permea = volume_flux*mü_water/9.81/density_cont # permeability in m2"""

                            # udpate to virtual test
                            v_permea = int_permea
                            # translation to rock dictionary
                            master_rock[item]['live_fluid-flux'].append(volume_flux)
                            master_rock[item]['live_permeability'].append(int_permea)

                            # old version
                            # v_permea = v_flux/100 * mü_water / tc_time_c / 9.81 / density_cont  # in m2
                            # print(f"-> Virtual permeability test results: {v_permea}")

                            # ##############################################
                            # LINK Coulomb method tests

                            if num > 0:

                                if failure_mech == 'Mohr-Coulomb-Permea2':
                                    print("\t===== Mohr-Couloumb.Permea2 method active =====")
                                    # Call the coulomb method no 2 - fixed diff stress failure test
                                    # LINK - diff stress input here
                                    master_rock[item]['fluid_calculation'].couloumb_method2(
                                        shear_stress=master_rock[item]['shear'],
                                        friction=master_rock[item]['friction'],
                                        cohesion=master_rock[item]['cohesion']
                                        )
                                elif failure_mech == 'Mohr-Coulomb-Griffith':
                                    print("\t===== Mohr-Coulomb-Griffith method active =====")
                                    #test if differntial stress is in master_rock[item] keys
                                    if 'diff. stress' in master_rock[item].keys():
                                        master_rock[item]['fluid_calculation'].mohr_cloulomb_griffith()
                                    else:
                                        master_rock[item]['fluid_calculation'].mohr_cloulomb_griffith(
                                            shear_stress=master_rock[item]['shear']
                                        )

                                    master_rock[item]['failure module'].append(
                                        master_rock[item]['fluid_calculation'].failure_dictionary)

                                # LINK ii) Steady state fluid extraction
                                elif failure_mech == 'Steady':
                                    print("===== steady method active =====")
                                    master_rock[item]['fluid_calculation'].frac_respo = 5
                                    fracturing_flag = True
                                    master_rock[item]['fluid_calculation'].fracture = True
                                    master_rock[item]['failure module'].append("Steady")

                                else:
                                    fracturing_flag = False
                                    master_rock[item]['fracture bool'].append(
                                            master_rock[item]['fluid_calculation'].frac_respo)

                            else:
                                fracturing_flag = False
                                master_rock[item]['fracture bool'].append(
                                        master_rock[item]['fluid_calculation'].frac_respo)

                            # ##############################################
                            # LINK Coulomb mechanical trigger - use trigger here 
                            # Tracking fracturing from coulomb approach methods
                            # Editing trigger
                            # if coulomb is True or coulomb_permea2 is True or coulomb_permea is True:

                            if failure_mech in ['Factor', 'Dynamic', 'Steady', 'Mohr-Coulomb-Permea2', 'Mohr-Coulomb-Griffith']:

                                """# store the differential stresses
                                master_rock[item]['diff. stress'].append(
                                    master_rock[item]['fluid_calculation'].diff_stress)"""

                                # store a bool index for the type of fracturing
                                master_rock[item]['fracture bool'].append(
                                    master_rock[item]['fluid_calculation'].frac_respo)
                                master_rock[item]['fracture_value'] = 1 + \
                                    master_rock[item]['tensile strength'] / \
                                    (pressures[num]/10)

                                # Fracture flag trigger
                                fracturing_flag = master_rock[item]['fluid_calculation'].fracture
                                # print(f"\nThe calculated extensional fracturing fator is: .... {fracture_value}\n")
                                # print(f"Check factor: {fluid_fac}")

                                # ##############################################
                                # LINK Release criteria
                                # Fluid Extraction when the modules before give true fracturing
                                # checking with the mohr-coloumb model and decision for fracturing or not

                                """if fracturing_flag is True and v_permea > lowest_permeability[jjj]:
                                    print("!!! Below minimum permeability!")"""
                                # # FIXME modified extraction criteria - minimum permeability is never reached 06.03.2023
                                if fracturing_flag is True:
                                    print("Enter fluid extraction")
                                    master_rock[item]['fluid_extraction'] = Fluid_master(
                                        phase_data=master_rock[item]['minimization'].df_phase_data.loc[:, fluid_name_tag],
                                        ext_data=master_rock[item]['extracted_fluid_data'],
                                        temperature=num+1,
                                        new_fluid_V=master_rock[item]['fluid_volume_new'],
                                        sys_H=master_rock[item]['total_hydrogen'],
                                        element_frame=master_rock[item]['df_element_total'],
                                        st_fluid_post=master_rock[item]['st_fluid_after'],
                                        fluid_name_tag=fluid_name_tag
                                        )
                                    # print("Call extraction function")
                                    print(f"*** Extracting fluid from {item}")

                                    if failure_mech == 'Mohr-Coulomb-Griffith' and master_rock[item]['fluid_calculation'].frac_respo == 5:
                                        master_rock[item]['fluid_extraction'].hydrogen_partial_ext(master_rock[item]['extraction treshold'])
                                    else:
                                        master_rock[item]['fluid_extraction'].hydrogen_ext_all(master_rock[item]['extraction percentage'])

                                    # print("Write data of extracted fluid")
                                    master_rock[item]['extracted_fluid_data'] = master_rock[item]['fluid_extraction'].ext_data
                                    # save time and system volume to list at extraction
                                    master_rock[item]['df_element_total'] = master_rock[item]['fluid_extraction'].element_frame
                                    master_rock[item]['extr_time'].append(track_time[num])
                                    # step system total volume (for surface ---> time integrated fluid flux)
                                    master_rock[item]['extr_svol'].append(
                                        np.sum(master_rock[item]['df_var_dictionary']['df_volume[ccm]'].iloc[:, -1]))
                                    master_rock[item]['track_refolidv'] = []


                                    # fractionate trace elements from the bulk
                                    if np.size(trace_df) > 0:
                                        #last_entry_key = data_storage_keys[-1]
                                        #last_entry_value = master_rock[item]['data_storage'][last_entry_key]

                                        master_rock[item]['init_trace_element_bulk'] = modify_trace_element_content(
                                        trace_element_bulk=master_rock[item]['init_trace_element_bulk'],
                                        trace_element_distirbution=trace_df,
                                        min_name = phase
                                        )

                                else:
                                    master_rock[item]['track_refolidv'].append(
                                        master_rock[item]['solid_volume_before'])
                                    master_rock[item]['fracture bool'][-1] = 0
                    # Starts no extraction scheme
                    else:
                        print("////// %s No extraction enabled! %s //////")


                    # //////////////////////////////////////////////////////////////////////////
                    # LINK Static: Recalculate bulk delta O after extraction
                    # Recalculate bulk rock oxygen value after possible extraction
                    new_O_bulk = oxygen_isotope_recalculation(master_rock[item]['save_oxygen'], master_rock[item]['df_element_total'])
                    master_rock[item]['bulk_oxygen'] = new_O_bulk
                    # bulk_oxygen = (rockOxy*bulk_oxygen - oxy_mole_fluid *
                    #                 fluid_oxygen)/(rockOxy - oxy_mole_fluid)
                
                else:
                        # No fluid in the system
                        if len(master_rock[item]['st_fluid_after']) > 1:
                            fluid_before = master_rock[item]['st_fluid_after'][-2]
                        else:
                            fluid_before = master_rock[item]['st_fluid_after'][-1]

                        master_rock[item]['fluid_calculation'] = Ext_method_master(
                            pressures[num], temperature,
                            0,
                            fluid_before, master_rock[item]['fluid_volume_new'],
                            master_rock[item]['solid_volume_before'], master_rock[item]['solid_volume_new'],
                            master_rock[item]['save_factor'], master_rock[item]['master_norm'][-1],
                            master_rock[item]['minimization'].df_phase_data,
                            master_rock[item]['tensile strength'],
                            differential_stress= master_rock[item]['diff. stress'],
                            friction= master_rock[item]['friction'],
                            fluid_pressure_mode= master_rock[item]['fluid_pressure_mode'],
                            fluid_name_tag=fluid_name_tag, subduction_angle=self.angle
                            )

                        if master_rock[item]['Extraction scheme'] == 'Mohr-Coulomb-Griffith':
                            print("\t===== Mohr-Coulomb-Griffith method active =====")
                            #test if differntial stress is in master_rock[item] keys
                            if 'diff. stress' in master_rock[item].keys():
                                master_rock[item]['fluid_calculation'].mohr_cloulomb_griffith()
                            else:
                                master_rock[item]['fluid_calculation'].mohr_cloulomb_griffith(
                                    shear_stress=master_rock[item]['shear']
                                )

                            # TODO - double asignment in output data, one is nan one the dry value - no solution found yet
                            """master_rock[item]['failure module'].append(
                                    master_rock[item]['fluid_calculation'].failure_dictionary)"""


                        master_rock[item]['save_factor'].append(0)
                        # master_rock[item]['diff. stress'].append(0)
                        master_rock[item]['fracture bool'].append(0)
                        master_rock[item]['live_fluid-flux'].append(np.nan)
                        master_rock[item]['live_permeability'].append(np.nan)
                        # master_rock[item]['failure module'].append("None activated, no fluid.")
                        # print(f"No free water in the system for {item} - no fracturing model")
                # save bulk oxygen after extraction
                master_rock[item]['save_bulk_oxygen_post'].append(
                    master_rock[item]['bulk_oxygen'])
                # print("\n")

            count += 1
            #k += 1
            """ic = k/kk*100
            print("=====Progress=====")
            progress(ic)"""
            print("\n")

        ar_flow = np.array(v_fluid_cubic_track)
        ar_perma = np.array(v_fluid_cubicp_track)

        """
        plt.figure()
        norm = np.array([
            0.3670, 0.9570, 0.1370, 0.7110, 0.2310, 0.0870, 0.3060, 
            0.0580, 0.3810, 0.0851, 0.2490, 0.0356, 0.2480, 0.0381])
        for i, value in enumerate(master_rock[item]['trace_element_bulk'].values()):
            test = value/norm
            plt.plot(test.T, '.-', label=i)
        # plt.xlabel(test.columns)
        plt.ylabel('Normalized values')
        plt.yscale('log')
        plt.legend()
        """

    def transmitting_multi_rock(self):
        """
        Perform multi-rock transmission calculations.

        This method calculates the transmission of multiple rocks based on various parameters such as temperature, pressure,
        rock composition, and mechanical methods. It iterates over the list of temperatures and performs calculations for each
        temperature.

        Args:
            self: The current object instance.

        Returns:
            None
        """

        # Main variables petrology
        temperatures = self.temperatures
        pressures = self.pressures
        master_rock = self.rock_dic
        rock_origin = self.rock_origin
        track_time = self.track_time
        track_depth = self.track_depth

        # Main variables mechanical model
        lowest_permeability = self.minimum_permeability

        """
        # REVIEW - mechanical methods
        # Methods
        factor_method = self.mechanical_methods[0]
        steady_method = self.mechanical_methods[1]
        dynamic_method = self.mechanical_methods[2]
        coulomb = self.mechanical_methods[3]
        coulomb_permea = self.mechanical_methods[4]
        coulomb_permea2 = self.mechanical_methods[5]
        """

        # initialize the trace element composition of the bulk rock
        for item in master_rock.keys():
            # trace element distribution
            master_rock[item]['init_trace_element_bulk'] = self.trace_element_bulk


        # Main variables fractionation
        grt_frac = self.garnet_fractionation
        print(f"Garnet fractionation in this model is: {grt_frac}")

        # //////////////////////////////////////////////////
        # /////////////////////////////////////////////////
        # ////////////////////////////////////////////////
        count = 0
        for num, temperature in enumerate(tqdm(temperatures, desc="Processing modelling steps")):


            print('\n')
            print("New calculation")
            print("Script: transmitting_multi_rock")
            print("===================")
            print(f"==== 1) time = {track_time[num]} years,\n==== 2) depth = {track_depth[num]}.")

            # //////////////////////////////////////////////////////////////////////////
            # preparing bulk rock for calculation
            rocks = list(master_rock.keys())
            for tt, item in enumerate(rocks):

                rock_react_item = list(master_rock.keys())[tt-1]

                if tt != 0:
                    # LINK - 0) Rock reactivity
                    # Decision for rock reactivity
                    # Reactivity decision is made on the rock below/before the momentary calculation (base rock vs. stack rock)

                    # taking the bulk rock elements and add the extracted fluid from layer below
                    if master_rock[rock_react_item]['reactivity'].react is True:

                        # Bulk rock calculation plus incoming H and O
                        if num < 1:
                            # FIXME - no fluid influx for rock in first step? Test this
                            master_rock[item]['new_bulk'], rockOxy = whole_rock_to_weight_normalizer(
                                rock_bulk=master_rock[item]['bulk'][:-2],
                                init_water=float(master_rock[item]['bulk'][-2]),
                                init_carbon=float(master_rock[item]['bulk'][-1])
                                )
                        else:
                            # print for testing which rock reacts into the new
                            print(f"Calculating {item}, {rock_react_item} is reactive")
                            if item != list(master_rock.keys())[tt]:
                                print("ERROR - reactivity: Momentary rock is not concordant!!!")
                                # keyboard.wait('esc')
                            # fix double total in column error
                            master_rock[item]['df_element_total']['total:']
                            bulka = master_rock[item]['df_element_total']['total:']

                            # double "total" problem
                            if isinstance(bulka, pd.DataFrame):
                                    if len(bulka.columns) > 1:
                                        bulka = bulka.iloc[:, 0]

                            # add the H and O that is transfered
                            # NOTE - Addition of moles dependend on geometry
                            # - need a factor because thermodynamic modellign uses 1 kg but geometry defines a volume
                            # - different geometries will have different impact on each other
                            external_rock_volume = (master_rock[rock_react_item]['st_solid'][-1] + master_rock[rock_react_item]['st_fluid_before'][-1])/1_000_000
                            external_rock_geometry = master_rock[rock_react_item]['geometry']
                            external_rock_geometry = np.float64(external_rock_geometry[0])*np.float64(external_rock_geometry[1])*np.float64(external_rock_geometry[2])
                            external_volume_change_factor = (
                                master_rock[rock_react_item]['st_solid'][-1] + master_rock[rock_react_item]['st_fluid_before'][-1]
                                                             ) / (
                                (master_rock[rock_react_item]['st_solid'][0] + master_rock[rock_react_item]['st_fluid_before'][0])
                                                             )

                            internal_volume = (master_rock[item]['fluid_volume_new'] + master_rock[item]['solid_volume_new'])/1_000_000
                            internal_geometry = master_rock[item]['geometry']
                            internal_geometry = np.float64(internal_geometry[0])*np.float64(internal_geometry[1])*np.float64(internal_geometry[2])
                            internal_volume_change_factor = (
                                master_rock[item]['fluid_volume_new'] + master_rock[item]['solid_volume_new']
                                                             ) / (
                                (master_rock[item]['st_solid'][0] + master_rock[item]['st_fluid_before'][0])
                                                             )

                            multiplicator_extern = external_rock_geometry / external_rock_volume

                            multiplicator_intern = internal_geometry / internal_volume

                            fluid_influx_factor = external_rock_geometry * internal_volume / external_rock_volume / internal_geometry * external_volume_change_factor / internal_volume_change_factor
                            print(f"Fluid influx factor is {fluid_influx_factor}")

                            # test if 'H' and 'O' are in the bulk rock index
                            if 'H' not in bulka.index:
                                bulka['H'] = 0
                            # Add the H and O to the bulk rock accounting for the geometry
                            bulka['H'] += (master_rock[rock_react_item]['fluid_hydrogen'][-1]*fluid_influx_factor)
                            bulka['O'] += (master_rock[rock_react_item]['fluid_oxygen'][-1]*fluid_influx_factor)
                            print("Fluid influx added to bulk rock")
                            print(f"Fluid influx is H = {master_rock[rock_react_item]['fluid_hydrogen'][-1]*fluid_influx_factor} mol and O = {master_rock[rock_react_item]['fluid_oxygen'][-1]*fluid_influx_factor} mol")
                            print(f"New calc H = {bulka['H']} mol and O = {bulka['O']} mol")
                            # save fluid influx to the rock for later use
                            master_rock[item]["fluid_influx_data"] = 0

                            # calculate the new bulk rock composition
                            master_rock[item]['new_bulk'] = whole_rock_convert_3(ready_mol_bulk=bulka)

                            # Recalculate bulk delta-oxygen after fluid input
                            new_O_bulk, fluid_end_d18O, free_fluid_oxygen = fluid_injection_isotope_recalculation(
                                        master_rock[item]['save_oxygen'],
                                        master_rock[item]['df_element_total'],
                                        master_rock[rock_react_item]['save_oxygen'][-1],
                                        master_rock[rock_react_item]['fluid_hydrogen'][-1]*fluid_influx_factor,
                                        master_rock[rock_react_item]['fluid_oxygen'][-1]*fluid_influx_factor,
                                        fluid_name_tag=master_rock[rock_react_item]['database_fluid_name']
                                        )

                            # Overwrite for new bulk rock oxygen signature
                            master_rock[item]['bulk_oxygen_before_influx'].append(np.copy(master_rock[item]['bulk_oxygen']))
                            master_rock[item]['bulk_oxygen'] = new_O_bulk
                            master_rock[item]['bulk_oxygen_after_influx'].append(np.copy(master_rock[item]['bulk_oxygen']))

                            print(f"Before bulk oxygen is {master_rock[item]['bulk_oxygen_before_influx'][-1]}")
                            print(f"New bulk oxygen is {master_rock[item]['bulk_oxygen_after_influx'][-1]}")

                            # recalculate the trace element bulk after infiltration
                            # ------------------------------------------------------
                            # access the fluid trace element data from the infiltrating fluid - calculate the influx amount multiplied by the geometry factor
                            last_entry_key = list(master_rock[rock_react_item]['trace_element_data'].keys())[-1]
                            tracer_addition = master_rock[rock_react_item]['trace_element_data'][
                                last_entry_key].loc[master_rock[rock_react_item]['database_fluid_name']] * fluid_influx_factor
                            # update the trace element bulk with the new influx
                            master_rock[item]['init_trace_element_bulk'] = master_rock[item]['init_trace_element_bulk'] + tracer_addition

                    else:
                        # Bulk rock calculation - normal
                        if num < 1:
                            master_rock[item]['new_bulk'], rockOxy = whole_rock_to_weight_normalizer(
                                rock_bulk=master_rock[item]['bulk'][:-2],
                                init_water=float(master_rock[item]['bulk'][-2]),
                                init_carbon=float(master_rock[item]['bulk'][-1])
                                )
                        else:
                            master_rock[item]['new_bulk'] = whole_rock_convert_3(
                                ready_mol_bulk=master_rock[item]['df_element_total']['total:']
                                )
                        # Storing delta-oxygen bulk rock in any case before fluid influx
                        master_rock[item]['bulk_oxygen_before_influx'].append(np.copy(master_rock[item]['bulk_oxygen']))
                        master_rock[item]['bulk_oxygen_after_influx'].append(np.copy(master_rock[item]['bulk_oxygen']))
                else:
                    # Bulk rock calculation - normal
                    if num < 1:
                        master_rock[item]['new_bulk'], rockOxy = whole_rock_to_weight_normalizer(
                            rock_bulk=master_rock[item]['bulk'][:-2],
                            init_water=float(master_rock[item]['bulk'][-2]),
                            init_carbon=float(master_rock[item]['bulk'][-1])
                            )
                    else:
                        master_rock[item]['new_bulk'] = whole_rock_convert_3(
                            ready_mol_bulk=master_rock[item]['df_element_total']['total:']
                            )
                    master_rock[item]['bulk_oxygen_before_influx'].append(0)
                    master_rock[item]['bulk_oxygen_after_influx'].append(0)

                print("¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦")
                print("v v v v v v v v v v v v v v v v v v v v v v v v")
                print(f"-> Forward modeling step initiated - {item}")
                    # print("Bulk rock stop")
                # store the current used bulk rock to the backup dictionary
                rock_origin[item]['bulk'].append(master_rock[item]['new_bulk'])
                # print(f"{item} Bulk rock composition checked.")
                print(f"Bulk rock composition checked. No error found")
                print("\n")


                # tracking theriak input before minimization
                if isinstance(master_rock[item]['theriak_input_record'], dict) == False:
                    master_rock[item]['theriak_input_record'] = {}
                    master_rock[item]['theriak_input_record']['temperature'] = [temperature]
                    master_rock[item]['theriak_input_record']['pressure'] = [pressures[num]]
                    master_rock[item]['theriak_input_record']['bulk'] = [master_rock[item]['new_bulk']]
                    # print("Tracking theriak -> Create dictionary -> first entry")

                # test for empty dictionary
                elif isinstance(master_rock[item]['theriak_input_record'], dict) == True:
                    master_rock[item]['theriak_input_record']['temperature'].append(temperature)
                    master_rock[item]['theriak_input_record']['pressure'].append(pressures[num])
                    master_rock[item]['theriak_input_record']['bulk'].append(master_rock[item]['new_bulk'])
                    print("Tracking theriak -> dictionary exists -> add entry")

                # _____________________________________________________________________________
                # 1) Initialize rock
                # LINK 1) Initialisation of the rock system
                master_rock[item]['minimization'] = Therm_dyn_ther_looper(self.theriak,
                    master_rock[item]['database'], master_rock[item]['new_bulk'],
                    temperature, pressures[num], master_rock[item]['df_var_dictionary'],
                    master_rock[item]['df_h2o_content_dic'], master_rock[item]['df_element_total'],
                    num, fluid_name_tag=master_rock[item]['database_fluid_name'])
                # //////////////////////////////////////////////////////////////////////////
                # Master norm values
                # calculating difference between new Volumes and previous P-T-step volumes - for derivate, not difference!!!
                if num < 1:
                    master_rock[item]['master_norm'].append(np.sqrt(
                        (temperature-temperatures[num])**2 + (pressures[num]-pressures[num])**2))
                else:
                    master_rock[item]['master_norm'].append(np.sqrt(
                        (temperature-temperatures[num-1])**2 + (pressures[num]-pressures[num-1])**2))

                # //////////////////////////////////////////////////////////////////////////
                # 1) Minimization
                # Calculating and passing thermo data by theriak and theriak wrapper
                # ----> self.df_var_dictionary, self.df_all_elements, self.df_hydrous_data_dic, new_fluid_volumes
                master_rock[item]['minimization'].thermodynamic_looping_station()
                # Main dictionary save
                # saving G_sys per mol of system
                master_rock[item]['g_sys'].append(master_rock[item]['minimization'].g_sys /
                                                sum(master_rock[item]['minimization'].df_phase_data.iloc[0, :]))
                # g_sys.append(minimization.g_sys)
                master_rock[item]['pot_data'].append(
                    master_rock[item]['minimization'].pot_frame)
                # Backup dictionary save
                rock_origin[item]['g_sys'] = copy.deepcopy(
                    master_rock[item]['g_sys'][-1])
                rock_origin[item]['pot_data'] = copy.deepcopy(
                    master_rock[item]['pot_data'][-1])

                # //////////////////////////////////////////////////////////////////////////
                # 2) Creating DataFrame structure for "df_var_dictionary" in first itteration
                # LINK - 2) Setup dataframes
                if num < 1:
                    # Volume and Density ouput - Dataframes (df_N, df_Vol% etc)
                    for variable in list(master_rock[item]['minimization'].df_phase_data.index):
                        master_rock[item]['df_var_dictionary']['df_' +
                                                            str(variable)] = pd.DataFrame()
                    water_cont_ind = ["N", "H2O[pfu]", "H2O[mol]",
                                    "H2O[g]", "wt%_phase", "wt%_solids", "wt%_H2O.solid"]
                    for variable in water_cont_ind:
                        master_rock[item]['df_h2o_content_dic']['df_' +
                                                                str(variable)] = pd.DataFrame()
                    # Copy to backup dictionary
                    rock_origin[item]['df_var_dictionary'] = copy.deepcopy(
                        master_rock[item]['df_var_dictionary'])
                    rock_origin[item]['df_h2o_content_dic'] = copy.deepcopy(
                        master_rock[item]['df_h2o_content_dic'])
                # updating dictionary with newly calculated data
                master_rock[item]['minimization'].merge_dataframe_dic()
                print("\n")
                print("////// Energy minimization executed //////")
                print("\n")

                # //////////////////////////////////////////////////////////////////////////
                # multi-rock loop for updating data storage, MicaPotassium, SystemFluidTest, init oxygen-isotope module, mineral fractionation
                # //////////////////////////////////////////////////////////////////////////
                # 3)
                # LINK - 3) Data storage & merge
                # calling dictionaries and dataframe for up-to-date usage
                # print(f"Running data re-storage, MicaPotassium, SystemFluidTest, oxy-module and mineral fractionation")
                # for item in master_rock:
                master_rock[item]['df_var_dictionary'], master_rock[item]['df_h2o_content_dic'], master_rock[item]['df_element_total'] = (
                    master_rock[item]['minimization'].df_var_dictionary,
                    master_rock[item]['minimization'].df_hydrous_data_dic,
                    master_rock[item]['minimization'].df_all_elements)
                master_rock[item]['df_element_total'] = master_rock[item]['df_element_total'].iloc[:, :-1]
                # hydrogen content of the system before extraction
                if 'H' in master_rock[item]['df_element_total'].index:
                    master_rock[item]['total_hydrogen'] = master_rock[item]['df_element_total']['total:']['H']
                else:
                    master_rock[item]['total_hydrogen'] = 0
                master_rock[item]['st_elements'] = pd.concat(
                    [master_rock[item]['st_elements'], master_rock[item]['df_element_total']['total:']], axis=1)
                # Backup dictionary - merging the data
                for kkey in rock_origin[item]['df_var_dictionary'].keys():
                    cdata = pd.concat(
                        [rock_origin[item]['df_var_dictionary'][kkey], master_rock[item]['df_var_dictionary'][kkey].iloc[:, -1]], axis=1)
                    rock_origin[item]['df_var_dictionary'][kkey] = copy.deepcopy(
                        cdata)
                for kkey in rock_origin[item]['df_h2o_content_dic'].keys():
                    # Testing if dataframe is empty
                    if rock_origin[item]['df_h2o_content_dic'][kkey].empty is True:
                        cdata = rock_origin[item]['df_h2o_content_dic'][kkey]
                    # If not empty, merge data
                    else:
                        cdata = pd.concat(
                            [rock_origin[item]['df_h2o_content_dic'][kkey], master_rock[item]['df_h2o_content_dic'][kkey].iloc[:, -1]], axis=1)
                    rock_origin[item]['df_h2o_content_dic'][kkey] = copy.deepcopy(
                        cdata)
                cdata = master_rock[item]['df_element_total']
                rock_origin[item]['df_element_total'].append(copy.deepcopy(cdata))
                # //////////////////////////////////////////////////////////////////////////
                # store Mica potassium if stable
                for phase in master_rock[item]['df_element_total'].columns:
                    if 'PHNG' in phase:
                        master_rock[item]['mica_K'].append(
                            [temperature, master_rock[item]['df_element_total'][phase]['K']])
                # //////////////////////////////////////////////////////////////////////////
                # 4)
                # LINK - 4) System and fluid volumes
                # Checking for fluid/solid volumes at t = 0 and t = -1,
                # calculating difference (used for escape/extraction rule e.g., factor method)
                # print("-> Testing for aq fluid in the system")
                master_rock[item]['minimization'].step_on_water()
                # Checking fluid and solid volumes. Storing t(-1) and t(0) data before calling calculation
                if num < 1:
                    master_rock[item]['fluid_volume_new'] = master_rock[item]['minimization'].new_fluid_Vol
                    master_rock[item]['fluid_volume_before'] = master_rock[item]['fluid_volume_new']
                    master_rock[item]['solid_volume_new'] = master_rock[item]['minimization'].new_solid_Vol
                    master_rock[item]['solid_volume_before'] = master_rock[item]['solid_volume_new']
                else:
                    master_rock[item]['fluid_volume_before'] = master_rock[item]['minimization'].free_water_before
                    master_rock[item]['fluid_volume_new'] = master_rock[item]['minimization'].new_fluid_Vol
                    # if master_rock[item]['minimization'].solid_vol_before != master_rock[item]['st_solid'][-1]:
                    #     print("\nWARNING: solid volume mismatch\n")
                    master_rock[item]['solid_volume_before'] = master_rock[item]['minimization'].solid_vol_before
                    master_rock[item]['solid_volume_before'] = master_rock[item]['st_solid'][-1]
                    master_rock[item]['solid_volume_new'] = master_rock[item]['minimization'].new_solid_Vol
                # //////////////////////////////////////////////////////////////////////////
                # 5)
                # LINK - 5) Oxygen fractionation module
                # isotope fractionation module
                # print("-> Oxygen isotope module initiated")

                # print("d18O before oxygen isotope module")
                # print(f"Value is {master_rock[item]['bulk_oxygen']}")
                master_rock[item]['model_oxygen'] = Isotope_calc(
                    master_rock[item]['df_var_dictionary']['df_N'], master_rock[item]['minimization'].sol_sol_base,
                    master_rock[item]['df_element_total'], oxygen_signature=master_rock[item]['bulk_oxygen'])
                # print(f"Bulk oxygen is {master_rock[item]['bulk_oxygen']}")
                # Isotope calc - function for oxygen isotope signatures
                master_rock[item]['model_oxygen'].frac_oxygen(temperature)
                # storing isotope fractionation result, dic in list appended
                master_rock[item]['save_oxygen'].append(
                    master_rock[item]['model_oxygen'].oxygen_dic)
                master_rock[item]['save_bulk_oxygen_pre'].append(
                    master_rock[item]['bulk_oxygen'])
                ### Backup dictionary - save oxygen data
                rock_origin[item]['save_oxygen'].append(copy.deepcopy(master_rock[item]['model_oxygen'].oxygen_dic))

                # LINK - 5-2) Trace element module
                master_rock[item]['model_tracers'] = TraceElementDistribution(
                    master_rock[item]['df_var_dictionary']['df_N'],
                    master_rock[item]['minimization'].sol_sol_base,
                    master_rock[item]['df_element_total'],
                    master_rock[item]['init_trace_element_bulk'],
                    database=master_rock[item]['database'])
                # call the distribution of tracers
                trace_df = master_rock[item]['model_tracers'].distribute_tracers(temperature, pressure=pressures[num], iteration=num)
                # save the tracer data
                master_rock[item]['trace_element_data'][(num, pressures[num], temperature)] = trace_df
                master_rock[item]['trace_element_bulk'][(num, pressures[num], temperature)] = master_rock[item]['init_trace_element_bulk']

                # !!! When fluid is consumed after interaction the d18O of the fluid in equilibirum with the system is defined by the equilibration calculation
                # fluid volume new < fluid volume extern + fluid volume before
                # taking the bulk rock elements and add the extracted fluid from layer below
                if master_rock[rock_react_item]['reactivity'].react is True and num > 0 and tt != 0:
                    if master_rock[item]['fluid_volume_new'
                                        ] <= master_rock[item]['fluid_volume_before'
                                                        ] + master_rock[rock_react_item]['extracted_fluid_data'
                                                                ].loc['volume[ccm]'].iloc[-1] * fluid_influx_factor:
                        pass

                    else:
                        if master_rock[item]['fluid_volume_new'
                                            ] > master_rock[item]['fluid_volume_before'
                                                        ] + master_rock[rock_react_item]['extracted_fluid_data'
                                                                ].loc['volume[ccm]'].iloc[-1] * fluid_influx_factor:

                            # internal oxygen moles and oxygen isotope signature
                            oxygen_dic = master_rock[item]['model_oxygen'].oxygen_dic
                            df_o_temp = pd.DataFrame(
                                    oxygen_dic['delta_O'],
                                    index=oxygen_dic['Phases']
                                        )
                            internal_fluid_d18o = df_o_temp.loc[fluid_name_tag][0]
                            internal_fluid_oxy = master_rock[item]['df_element_total'][fluid_name_tag].loc['O']

                            # external oxygen moles         = free_fluid_oxygen
                            # interacted isotope signature  = fluid_end_d18O

                            """oxygen_dic = master_rock[rock_react_item]['save_oxygen'][-1]
                            df_o_temp = pd.DataFrame(
                                    oxygen_dic['delta_O'],
                                    index=oxygen_dic['Phases']
                                    )
                            input_fluid_d18o = df_o_temp.loc[fluid_name_tag][0]
                            input_oxygen = master_rock[rock_react_item]['fluid_oxygen'][-1]*fluid_influx_factor"""


                            fluid_mix_d18O = (free_fluid_oxygen*fluid_end_d18O + (internal_fluid_oxy-free_fluid_oxygen)*internal_fluid_d18o) / (free_fluid_oxygen+(internal_fluid_oxy-free_fluid_oxygen))
                            print(f"Fluid mix update d18O is: {fluid_mix_d18O}")

                            # index number of fluid_name_tag in Phases of oxygen_dic
                            index = oxygen_dic['Phases'].index(fluid_name_tag)

                            # overwriting the fluid oxygen isotope composition to the one calculated by the mixing
                            rock_origin[item]['save_oxygen'][-1]['delta_O'][index] = fluid_mix_d18O



                print("d18O after oxygen isotope module")
                print(f"Value is {master_rock[item]['bulk_oxygen']}")
                # //////////////////////////////////////////////////////////////////////////
                # 6)
                # LINK - 6) Mineral Fractionation
                # mineral (garnet) fractionation - coupled oxygen bulk modification
                if grt_frac == True:
                    if master_rock[item]['database'] == 'ds62mp.txt' or master_rock[item]['database'] == 'td-ds62-mb50-v07.txt':
                        garnet_name = 'GRT'
                    else:
                        garnet_name = 'GARNET'
                    # print("-> Mineral fractionation initiated")
                    # old frac position this line
                    for phase in master_rock[item]['df_element_total'].columns:
                        if '_' in phase:
                            pos = phase.index('_')
                            name = phase[:pos]
                            if name == garnet_name:
                                # Michelles atigorite fractionation
                                # # if name=='GARNET' or name=='SERP' or name=='BR':
                                new_bulk_oxygen = master_rock[item]['minimization'].mineral_fractionation(
                                    master_rock[item]['save_oxygen'][-1], name)

                                # modify the trace element content
                                if np.size(trace_df) > 0:
                                    #last_entry_key = data_storage_keys[-1]
                                    #last_entry_value = master_rock[item]['data_storage'][last_entry_key]

                                    master_rock[item]['init_trace_element_bulk'] = modify_trace_element_content(
                                        trace_element_bulk=master_rock[item]['init_trace_element_bulk'],
                                        trace_element_distirbution=trace_df,
                                        min_name = phase
                                    )

                                master_rock[item]['garnet'].append(
                                    master_rock[item]['minimization'].separate)
                                # Collect the volume of each formed garnet
                                master_rock[item]['meta_grt_volume'].append(master_rock[item]['minimization'].separate.volume)
                                # print(
                                #     f"Selected phase = {phase} with Vol% = {master_rock[item]['df_var_dictionary']['df_vol%'].loc[phase][-2:]}")
                                # print(
                                #     f"Bulk deltaO changed from {round(master_rock[item]['bulk_oxygen'], 3)} to {round(new_bulk_oxygen, 3)}")
                                master_rock[item]['bulk_oxygen'] = new_bulk_oxygen
                                print("_______________________")
                                master_rock[item]['garnet_check'].append(1)
                            if len(master_rock[item]['garnet_check']) < num:
                                master_rock[item]['garnet_check'].append(0)

                print("d18O after mineral fractionation module")
                print(f"Value is {master_rock[item]['bulk_oxygen']}")

                # Do not delete - necessary step -
                master_rock[item]['df_element_total'] = master_rock[item]['minimization'].df_all_elements
                # LINK - 7) Metastable garnet
                # Adding the metastable garnet impact
                # calculate the metastable garnet for all bits besides the last
                # add the calculated volume to the solid volume of the current step (this is then saved to the st_solid and used next turn)
                grt_flag = True
                if grt_flag is True and len(master_rock[item]['garnet']) > 0 and len(master_rock[item]['garnet_check']) > 1:
                    if master_rock[item]['garnet_check'][-1] == 0:
                        print("Garnet protocol")
                        # take al modelled garnets
                        # LINK - Metastable garnet call
                        metastable_garnet = Garnet_recalc(self.theriak, master_rock[item]['garnet'], temperature, pressures[num])
                        metastable_garnet.recalculation_of_garnets(database=master_rock[item]['database'], garnet_name=garnet_name)
                        print(f"Fluid volume = {master_rock[item]['fluid_volume_new']} ccm")
                        print(f"Solid volume = {master_rock[item]['solid_volume_new']} ccm")
                        volume = metastable_garnet.recalc_volume
                        # Adding the metastable volume of the fractionated garnet to the current minimized solied volume!
                        master_rock[item]['solid_volume_new'] += volume
                    if len(master_rock[item]['garnet']) > 1 and master_rock[item]['garnet_check'][-1] == 1:
                        print("protocol")
                        # take all garnets but last one
                        metastable_garnet = Garnet_recalc(self.theriak, master_rock[item]['garnet'][:-1], temperature, pressures[num])
                        metastable_garnet.recalculation_of_garnets(database=master_rock[item]['database'], garnet_name=garnet_name)
                        print(f"Fluid volume = {master_rock[item]['fluid_volume_new']} ccm")
                        print(f"Solid volume = {master_rock[item]['solid_volume_new']} ccm")
                        volume = metastable_garnet.recalc_volume
                        # Adding the metastable volume of the fractionated garnet to the current minimized solied volume!
                        master_rock[item]['solid_volume_new'] += volume
                        master_rock[item]['meta_grt_weight'].append(metastable_garnet.recalc_weight)
                else:
                    master_rock[item]['solid_volume_new'] += np.array(master_rock[item]['meta_grt_volume']).sum()
                    if np.array(master_rock[item]['meta_grt_volume']).sum() != np.float64(0):
                        print("What happend now? Here is suddenly garnet which is not stable???")
                        # keyboard.wait('esc')
                    metastable_garnet_weight = 0

                # keeping track of stored and removed fluid
                # !!!! necessary for mechanical module
                master_rock[item]['st_fluid_before'].append(
                    master_rock[item]['fluid_volume_new'])
                master_rock[item]['st_fluid_after'].append(
                    master_rock[item]['fluid_volume_new'])
                master_rock[item]['st_solid'].append(
                    master_rock[item]['solid_volume_new'])

                # //////////////////////////////////////////////////////////////////////////
                # LINK - MECHANICAL FAILURE MODEL
                # 1) Failure and extraction mode selection
                # 2) Oxygen isotope recalculation
                # Physical fracture model section - active if free fluid is present
                # Checking if free fluid is present. Stores this value, initializes
                print(f"Mechanical failure model activated.")
                fluid_name_tag = master_rock[item]['database_fluid_name']
                if fluid_name_tag in list(master_rock[item]['df_element_total'].columns):
                    print("-> Fluid extraction test")
                    # Prepare fluid extraction
                    if master_rock[rock_react_item]['reactivity'].react is True:
                            print("reactivity break")

                    # Prepare fluid extraction
                    # - last entry of st_fluid_after has last fluid volume (number when not extracted or zero when extracted)
                    if len(master_rock[item]['st_fluid_after']) > 1:
                        fluid_before = master_rock[item]['st_fluid_after'][-2]
                    else:
                        fluid_before = master_rock[item]['st_fluid_after'][-1]

                    # Start Extraction Master Module
                    master_rock[item]['fluid_calculation'] = Ext_method_master(
                        pressures[num], temperature,
                        master_rock[item]['df_var_dictionary']['df_volume/mol'].loc[fluid_name_tag].iloc[-1],
                        fluid_before, master_rock[item]['fluid_volume_new'],
                        master_rock[item]['solid_volume_before'], master_rock[item]['solid_volume_new'],
                        master_rock[item]['save_factor'], master_rock[item]['master_norm'][-1],
                        master_rock[item]['minimization'].df_phase_data,
                        master_rock[item]['tensile strength'],
                        differential_stress= master_rock[item]['diff. stress'],
                        friction= master_rock[item]['friction'],
                        fluid_pressure_mode= master_rock[item]['fluid_pressure_mode'],
                        fluid_name_tag=fluid_name_tag, subduction_angle=self.angle,
                        rock_item_tag=item,
                        extraction_treshold = master_rock[item]['extraction treshold'],
                                        extraction_connectivity = master_rock[item]['fluid connectivity']
                        )
                    # //////////////////////////////////////////////////////////////////////////
                    # LINK 1) selection of the failure and fluid extraction
                    # ////// Calculation for new whole rock /////////
                    # !!!!Assumption: fracturing of the system
                    # if condition = open system, free water gets extracted
                    fracturing_flag = False # Trigger by default False - active when coulomb module becomes positive
                    failure_mech = master_rock[item]['Extraction scheme']
                    if failure_mech in ['Factor', 'Dynamic', 'Steady', 'Mohr-Coulomb-Permea2', 'Mohr-Coulomb-Griffith']:
                        # ANCHOR work in progres GRIFFITH
                        # if factor_method is True or dynamic_method is True or steady_method is True or coulomb is True or coulomb_permea is True or coulomb_permea2 is True:
                        # LINK Fluid flux & permeability
                        # Fluid flux check - Virtual calculation
                        # Hypothetically momentary fluid flux and permeabiltiy test
                        mü_water = 1e-4
                        # water data
                        v_water = float(master_rock[item]['df_var_dictionary']['df_volume[ccm]'].loc[fluid_name_tag].iloc[-1])
                        d_water = float(master_rock[item]['df_var_dictionary']['df_density[g/ccm]'].loc[fluid_name_tag].iloc[-1])
                        weigth_water = master_rock[item]['df_var_dictionary']['df_wt[g]'].iloc[:, -1][fluid_name_tag]
                        # system data
                        master_rock[item]['meta_grt_weight']
                        v_system = master_rock[item]['solid_volume_new'] + master_rock[item]['fluid_volume_new'] # modelling solid phases + metastable garnet + fluid
                        # add metastable garnet weight - volume is already up to date
                        if len(master_rock[item]['meta_grt_weight']) > 0:
                            weight_sys = float(
                                master_rock[item]['df_var_dictionary']['df_wt[g]'].iloc[:, -1].sum()
                                ) + master_rock[item]['meta_grt_weight'][-1] # weight of solid + fluid + metastable garnet
                        else:
                            weight_sys = float(
                                master_rock[item]['df_var_dictionary']['df_wt[g]'].iloc[:, -1].sum()
                                ) # weight of solid + fluid + metastable garnet

                        d_system = float(master_rock[item]['df_var_dictionary']['df_wt[g]'].iloc[:, -1].sum())/v_system
                        # solid rock data
                        v_rock = master_rock[item]['solid_volume_new']
                        weight_rock = weight_sys - weigth_water # weight of solids + metastable garnet
                        d_rock = weight_rock/v_rock
                        # density difference
                        density_cont = (d_rock-d_water)*1000 # kg/m3

                        # Reading the time interval
                        if num < 1:
                            tc_time_c = float(np.diff(track_time[num:num+2])*365*24*60*60)
                        else:
                            tc_time_c = float(np.diff(track_time[num-1:num+1])*365*24*60*60)


                        """# test 4
                        cubic_meter = 2 # TODO put as input here - geometry or size of the block!
                        xxx = np.power(1_000_000 * cubic_meter, 1/3) # length of x of cube (100cm for 1 m3)
                        size = xxx**3 # cubic size
                        area = xxx**2 # surface size
                        # fluid flux = drainiage flux throughout a column
                        volume_flux = v_water * size /v_system/area/tc_time_c # cm3 cm-2 s-1
                        volume_flux = volume_flux/100 # m3 m-2 s-1
                        # integrated permeability
                        int_permea = volume_flux*mü_water/9.81*xxx/density_cont # permeability in m2"""

                        # test 05
                        bloc_a = np.float64(master_rock[item]['geometry'][0])
                        bloc_b = np.float64(master_rock[item]['geometry'][1])
                        bloc_c = np.float64(master_rock[item]['geometry'][2])
                        area = bloc_b*bloc_c
                        xxx = bloc_a
                        size = bloc_a * bloc_b * bloc_c
                        v_water1 = v_water/1000000 # cm3 to m3
                        v_system1 = v_system/1000000 # cm3 to m3
                        volume_flux = v_water1 * size/v_system1/area/tc_time_c # m3 m-2 s-1
                        int_permea = volume_flux*mü_water/9.81/xxx/density_cont # permeability in m2

                        # udpate to virtual test
                        v_permea = int_permea
                        master_rock[item]['live_fluid-flux'].append(volume_flux)
                        master_rock[item]['live_permeability'].append(int_permea)
                        # print(f"-> Virtual permeability test results: {v_permea}")

                        # Latest failure criterion 07.02.2023
                        # LINK i) Mohr-Coulomb failure w. diff. stress input and min. permeabiltiy
                        if failure_mech == 'Mohr-Coulomb-Permea2':
                            print("\t===== Mohr-Couloumb.Permea2 method active =====")
                            master_rock[item]['fluid_calculation'].couloumb_method2(
                                shear_stress=master_rock[item]['shear'],
                                friction=master_rock[item]['friction'],
                                cohesion=master_rock[item]['cohesion']
                                )

                        elif failure_mech == 'Mohr-Coulomb-Griffith':
                            print("\t===== Mohr-Coulomb-Griffith method active =====")
                            if 'diff. stress' in master_rock[item].keys():
                                master_rock[item]['fluid_calculation'].mohr_cloulomb_griffith()
                            else:
                                master_rock[item]['fluid_calculation'].mohr_cloulomb_griffith(
                                        shear_stress=master_rock[item]['shear']
                                        )
                            master_rock[item]['failure module'].append(
                                master_rock[item]['fluid_calculation'].failure_dictionary)

                        # LINK ii) Steady state fluid extraction
                        elif failure_mech == 'Steady':
                            print("===== steady method active =====")
                            master_rock[item]['fluid_calculation'].frac_respo = 5
                            fracturing_flag = True
                            master_rock[item]['fluid_calculation'].fracture = True
                            master_rock[item]['failure module'].append("Steady")

                        else:
                            fracturing_flag = False
                            master_rock[item]['fracture bool'].append(
                                    master_rock[item]['fluid_calculation'].frac_respo)

                        # ##############################################
                        # LINK Coulomb mechanical trigger
                        # Tracking fracturing from coulomb approach methods
                        # Editing trigger
                        # if coulomb is True or coulomb_permea2 is True or coulomb_permea is True:
                        if failure_mech in ['Factor', 'Dynamic', 'Steady', 'Mohr-Coulomb-Permea2', 'Mohr-Coulomb-Griffith']:

                            # store a bool index for the type of fracturing
                            master_rock[item]['fracture bool'].append(
                                master_rock[item]['fluid_calculation'].frac_respo)
                            master_rock[item]['fracture_value'] = 1 + \
                                master_rock[item]['tensile strength'] / \
                                (pressures[num]/10)

                            # Fracture flag trigger
                            # activate fracture flag trigger when vol% of fluid is above or equal 10%
                            if v_water/v_system >= 0.1:
                                master_rock[item]['fluid_calculation'].fracture = True
                                fracturing_flag = True
                                master_rock[item]['fracture bool'][-1] = 10
                            else:
                                fracturing_flag = master_rock[item]['fluid_calculation'].fracture
                            # print(f"\nThe calculated extensional fracturing fator is: .... {fracture_value}\n")
                            # print(f"Check factor: {fluid_fac}")

                            # ##############################################
                            # LINK Release criteria
                            # Fluid Extraction when the modules before give true fracturing
                            # checking with the mohr-coloumb model and decision for fracturing or not
                            """if fracturing_flag is True and v_permea > lowest_permeability[tt]:
                                print("!!! Below minimum permeability!")"""
                            # FIXME modified extraction criteria - minimum permeability is never reached 06.03.2023
                            if fracturing_flag is True:
                                print("Enter fluid extraction")
                                # keyboard activation and exite by esc
                                # if master_rock[rock_react_item]['reactivity'].react is True:
                                #     print("Keyboard wait exception")
                                #     keyboard.wait('esc')
                                master_rock[item]['fluid_extraction'] = Fluid_master(
                                    phase_data=master_rock[item]['minimization'].df_phase_data.loc[:, fluid_name_tag],
                                    ext_data=master_rock[item]['extracted_fluid_data'],
                                    temperature=num+1,
                                    new_fluid_V=master_rock[item]['fluid_volume_new'],
                                    sys_H=master_rock[item]['total_hydrogen'],
                                    element_frame=master_rock[item]['df_element_total'],
                                    st_fluid_post=master_rock[item]['st_fluid_after'],
                                    fluid_name_tag=fluid_name_tag
                                    )
                                # backup before extraction
                                master_rock[item]['fluid_hydrogen'].append(master_rock[item]['df_element_total'][fluid_name_tag]['H'].copy())
                                master_rock[item]['fluid_oxygen'].append(master_rock[item]['df_element_total'][fluid_name_tag]['O'].copy())
                                # Execute the extraction
                                if failure_mech == 'Mohr-Coulomb-Griffith' and master_rock[item]['fluid_calculation'].frac_respo == 5:
                                        master_rock[item]['fluid_extraction'].hydrogen_partial_ext(master_rock[item]['extraction treshold'])
                                else:
                                    master_rock[item]['fluid_extraction'].hydrogen_ext_all(master_rock[item]['extraction percentage'])

                                # Read the data of the fluid extracted
                                master_rock[item]['extracted_fluid_data'] = master_rock[item]['fluid_extraction'].ext_data
                                # Read the element frame when extraction
                                master_rock[item]['df_element_total'] = master_rock[item]['fluid_extraction'].element_frame

                                # Control the rock reactivity
                                master_rock[item]['reactivity'].react=True
                                # Save the time of extraction step
                                master_rock[item]['extr_time'].append(track_time[num])
                                # step system total volume (for surface ---> time integrated fluid flux)
                                master_rock[item]['extr_svol'].append(
                                    np.sum(master_rock[item]['df_var_dictionary']['df_volume[ccm]'].iloc[:, -1]))
                                master_rock[item]['track_refolidv'] = []

                                # fractionate trace elements from the bulk
                                if np.size(trace_df) > 0:
                                    #last_entry_key = data_storage_keys[-1]
                                    #last_entry_value = master_rock[item]['data_storage'][last_entry_key]

                                    master_rock[item]['init_trace_element_bulk'] = modify_trace_element_content(
                                        trace_element_bulk=master_rock[item]['init_trace_element_bulk'],
                                        trace_element_distirbution=trace_df,
                                        min_name = phase
                                        )

                            else:
                                print("!!! No release!")
                                master_rock[item]['reactivity'].react=False
                                master_rock[item]['track_refolidv'].append(
                                    master_rock[item]['solid_volume_before'])
                                master_rock[item]['fracture bool'][-1] = 0

                    # OPTION no extraction
                    # Starts no extraction scheme
                    else:
                        master_rock[item]['reactivity'].react=False
                        print("////// %s No extraction enabled! %s //////")

                    # //////////////////////////////////////////////////////////////////////////
                    # LINK Recalculate the oxygen isotope signature
                    # Recalculate bulk rock oxygen value after possible extraction
                    print("Oxygen isotope signature recalculation fluid extraction - before")
                    print(f"Value is {master_rock[item]['bulk_oxygen']}")
                    new_O_bulk = oxygen_isotope_recalculation(
                        master_rock[item]['save_oxygen'],
                        master_rock[item]['df_element_total'])
                    # Overwrite for new bulk rock oxygen signature
                    master_rock[item]['bulk_oxygen'] = new_O_bulk
                    # bulk_oxygen = (rockOxy*bulk_oxygen - oxy_mole_fluid *
                    #                 fluid_oxygen)/(rockOxy - oxy_mole_fluid)
                    print("Oxygen isotope signature recalculation fluid extraction - after")
                    print(f"Value is {master_rock[item]['bulk_oxygen']}")
                else:
                    master_rock[item]['save_factor'].append(0)
                    master_rock[item]['fracture bool'].append(0)
                    master_rock[item]['live_fluid-flux'].append(np.nan)
                    master_rock[item]['live_permeability'].append(np.nan)
                    master_rock[item]['failure module'].append("None activated, no fluid.")
                    master_rock[item]['reactivity'].react=False
                    print(f"No free water in the system for {item} - no fracturing model")
                # save bulk oxygen after extraction
                master_rock[item]['save_bulk_oxygen_post'].append(
                    master_rock[item]['bulk_oxygen'])
                print("\n")

            count += 1

        # LINK ROUTINE END
        # //////////////////////////////////////////////////////////////////////////
        # //////////////////////////// END OF SCRIPT ///////////////////////////////
        # //////////////////////////////////////////////////////////////////////////

        self.rock_dic = master_rock

    def transmitting_multi_rock_altPT(self):
        """
        Perform multi-rock transmission calculations.
        Each rock model has its own PT-path.

        This method calculates the transmission of multiple rocks based on various parameters such as temperature, pressure,
        rock composition, and mechanical methods. It iterates over the list of temperatures and performs calculations for each
        temperature.

        Args:
            self: The current object instance.

        Returns:
            None
        """
        # save the list self.debugging_recorder to txt file
        with open('debugging_recorder.txt', 'w') as f:
            for item in self.debugging_recorder:
                f.write("%s\n" % item)


        # Main variables petrology
        temperatures = self.temperatures
        pressures = self.pressures
        master_rock = self.rock_dic
        rock_origin = self.rock_origin
        track_time = self.track_time
        track_depth = self.track_depth

        # Main variables mechanical model
        lowest_permeability = self.minimum_permeability

        """
        # REVIEW - mechanical methods
        # Methods
        factor_method = self.mechanical_methods[0]
        steady_method = self.mechanical_methods[1]
        dynamic_method = self.mechanical_methods[2]
        coulomb = self.mechanical_methods[3]
        coulomb_permea = self.mechanical_methods[4]
        coulomb_permea2 = self.mechanical_methods[5]
        """


        # Main variables fractionation
        grt_frac = self.garnet_fractionation

        # //////////////////////////////////////////////////
        # /////////////////////////////////////////////////
        # ////////////////////////////////////////////////
        count = 0
        k = 0
        kk = len(temperatures[0])*len(master_rock)
        progress(int(k/kk)*100)
        for num, temperature in enumerate(temperatures):

            print('\n')
            print("New calculation")
            print("Script: transmitting_multi_rock")
            print("===================")
            print(f"==== 1) time = {track_time[num]} years,\n==== 2) depth = {track_depth[num]}.")

            # //////////////////////////////////////////////////////////////////////////
            # preparing bulk rock for calculation
            rocks = list(master_rock.keys())
            for tt, item in enumerate(rocks):

                rock_react_item = list(master_rock.keys())[tt-1]

                if tt != 0:
                    # LINK - 0) Rock reactivity
                    # Decision for rock reactivity
                    # Reactivity decision is made on the rock below/before the momentary calculation (base rock vs. stack rock)

                    # taking the bulk rock elements and add the extracted fluid from layer below
                    if master_rock[rock_react_item]['reactivity'].react is True:

                        # Bulk rock calculation plus incoming H and O
                        if num < 1:
                            # FIXME - no fluid influx for rock in first step? Test this
                            master_rock[item]['new_bulk'], rockOxy = whole_rock_to_weight_normalizer(
                                rock_bulk=master_rock[item]['bulk'][:-2],
                                init_water=float(master_rock[item]['bulk'][-2]),
                                init_carbon=float(master_rock[item]['bulk'][-1])
                                )
                        else:
                            # print for testing which rock reacts into the new
                            print(f"Calculating {item}, {rock_react_item} is reactive")
                            if item != list(master_rock.keys())[tt]:
                                print("ERROR - reactivity: Momentary rock is not concordant!!!")
                                # keyboard.wait('esc')
                            # fix double total in column error
                            master_rock[item]['df_element_total']['total:']
                            bulka = master_rock[item]['df_element_total']['total:']

                            # double "total" problem
                            if isinstance(bulka, pd.DataFrame):
                                    if len(bulka.columns) > 1:
                                        bulka = bulka.iloc[:, 0]

                            # add the H and O that is transfered
                            # FIXME - Addition of moles dependend on geometry
                            # - need a factor because thermodynamic modellign uses 1 kg but geometry defines a volume
                            # - different geometries will have different impact on each other
                            external_rock_volume = (master_rock[rock_react_item]['st_solid'][-1] + master_rock[rock_react_item]['st_fluid_before'][-1])/1_000_000
                            external_rock_geometry = master_rock[rock_react_item]['geometry']
                            external_rock_geometry = np.float64(external_rock_geometry[0])*np.float64(external_rock_geometry[1])*np.float64(external_rock_geometry[2])

                            internal_volume = (master_rock[item]['fluid_volume_new'] + master_rock[item]['solid_volume_new'])/1_000_000
                            internal_geometry = master_rock[item]['geometry']
                            internal_geometry = np.float64(internal_geometry[0])*np.float64(internal_geometry[1])*np.float64(internal_geometry[2])

                            fluid_influx_factor = external_rock_geometry * internal_volume / external_rock_volume / internal_geometry
                            print(f"Fluid influx factor is {fluid_influx_factor}")

                            # test if 'H' and 'O' are in the bulk rock index
                            if 'H' not in bulka.index:
                                bulka['H'] = 0
                            # Add the H and O to the bulk rock accounting for the geometry
                            bulka['H'] += (master_rock[rock_react_item]['fluid_hydrogen'][-1]*fluid_influx_factor)
                            bulka['O'] += (master_rock[rock_react_item]['fluid_oxygen'][-1]*fluid_influx_factor)
                            print("Fluid influx added to bulk rock")
                            print(f"Fluid influx is {master_rock[rock_react_item]['fluid_hydrogen'][-1]*fluid_influx_factor} mol H and {master_rock[rock_react_item]['fluid_oxygen'][-1]*fluid_influx_factor} mol O")

                            # save fluid influx to the rock for later use
                            master_rock[item]["fluid_influx_data"] = 0

                            # calculate the new bulk rock composition
                            master_rock[item]['new_bulk'] = whole_rock_convert_3(ready_mol_bulk=bulka)

                            # Recalculate bulk delta-oxygen after fluid input
                            new_O_bulk = fluid_injection_isotope_recalculation(
                                        master_rock[item]['save_oxygen'],
                                        master_rock[item]['df_element_total'],
                                        master_rock[rock_react_item]['save_oxygen'][-1],
                                        master_rock[rock_react_item]['fluid_hydrogen'][-1]*fluid_influx_factor,
                                        master_rock[rock_react_item]['fluid_oxygen'][-1]*fluid_influx_factor,
                                        fluid_name_tag=master_rock[rock_react_item]['database_fluid_name']
                                        )

                            # Overwrite for new bulk rock oxygen signature
                            master_rock[item]['bulk_oxygen_before_influx'].append(np.copy(master_rock[item]['bulk_oxygen']))
                            master_rock[item]['bulk_oxygen'] = new_O_bulk
                            master_rock[item]['bulk_oxygen_after_influx'].append(np.copy(master_rock[item]['bulk_oxygen']))

                            print(f"New bulk oxygen is {master_rock[item]['bulk_oxygen_before_influx'][-1]}")
                            print(f"New bulk oxygen is {master_rock[item]['bulk_oxygen_after_influx'][-1]}")

                    else:
                        # Bulk rock calculation - normal
                        if num < 1:
                            master_rock[item]['new_bulk'], rockOxy = whole_rock_to_weight_normalizer(
                                rock_bulk=master_rock[item]['bulk'][:-2],
                                init_water=float(master_rock[item]['bulk'][-2]),
                                init_carbon=float(master_rock[item]['bulk'][-1])
                                )
                        else:
                            master_rock[item]['new_bulk'] = whole_rock_convert_3(
                                ready_mol_bulk=master_rock[item]['df_element_total']['total:']
                                )
                        # Storing delta-oxygen bulk rock in any case before fluid influx
                        master_rock[item]['bulk_oxygen_before_influx'].append(np.copy(master_rock[item]['bulk_oxygen']))
                        master_rock[item]['bulk_oxygen_after_influx'].append(np.copy(master_rock[item]['bulk_oxygen']))
                else:
                    # Bulk rock calculation - normal
                    if num < 1:
                        master_rock[item]['new_bulk'], rockOxy = whole_rock_to_weight_normalizer(
                            rock_bulk=master_rock[item]['bulk'][:-2],
                            init_water=float(master_rock[item]['bulk'][-2]),
                            init_carbon=float(master_rock[item]['bulk'][-1])
                            )
                    else:
                        master_rock[item]['new_bulk'] = whole_rock_convert_3(
                            ready_mol_bulk=master_rock[item]['df_element_total']['total:']
                            )
                    master_rock[item]['bulk_oxygen_before_influx'].append(0)
                    master_rock[item]['bulk_oxygen_after_influx'].append(0)

                print("¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦ ¦")
                print("v v v v v v v v v v v v v v v v v v v v v v v v")
                print(f"-> Forward modeling step initiated - {item}")
                    # print("Bulk rock stop")
                # store the current used bulk rock to the backup dictionary
                rock_origin[item]['bulk'].append(master_rock[item]['new_bulk'])
                # print(f"{item} Bulk rock composition checked.")
                print(f"Bulk rock composition checked. No error found")
                print("\n")
                # display modelling progress
                ic = k/kk*100
                progress(ic)
                k += 1
                print("\n")

                # tracking theriak input before minimization
                if isinstance(master_rock[item]['theriak_input_record'], dict) == False:
                    master_rock[item]['theriak_input_record'] = {}
                    master_rock[item]['theriak_input_record']['temperature'] = [temperature[tt]]
                    master_rock[item]['theriak_input_record']['pressure'] = [pressures[num][tt]]
                    master_rock[item]['theriak_input_record']['bulk'] = [master_rock[item]['new_bulk']]
                    print("Tracking theriak -> Create dictionary -> first entry")

                # test for empty dictionary
                elif isinstance(master_rock[item]['theriak_input_record'], dict) == True:
                    master_rock[item]['theriak_input_record']['temperature'].append(temperature[tt])
                    master_rock[item]['theriak_input_record']['pressure'].append(pressures[num][tt])
                    master_rock[item]['theriak_input_record']['bulk'].append(master_rock[item]['new_bulk'])
                    print("Tracking theriak -> dictionary exists -> add entry")

                # _____________________________________________________________________________
                # 1) Initialize rock
                # LINK 1) Initialisation of the rock system
                master_rock[item]['minimization'] = Therm_dyn_ther_looper(self.theriak,
                    master_rock[item]['database'], master_rock[item]['new_bulk'],
                    temperature[tt], pressures[num][tt], master_rock[item]['df_var_dictionary'],
                    master_rock[item]['df_h2o_content_dic'], master_rock[item]['df_element_total'],
                    num, fluid_name_tag=master_rock[item]['database_fluid_name'])
                # //////////////////////////////////////////////////////////////////////////
                # Master norm values
                # calculating difference between new Volumes and previous P-T-step volumes - for derivate, not difference!!!
                if num < 1:
                    master_rock[item]['master_norm'].append(np.sqrt(
                        (temperature[tt]-temperatures[num][tt])**2 + (pressures[num][tt]-pressures[num][tt])**2))
                else:
                    master_rock[item]['master_norm'].append(np.sqrt(
                        (temperature-temperatures[num-1][tt])**2 + (pressures[num][tt]-pressures[num-1][tt])**2))

                # //////////////////////////////////////////////////////////////////////////
                # 1) Minimization
                # Calculating and passing thermo data by theriak and theriak wrapper
                # ----> self.df_var_dictionary, self.df_all_elements, self.df_hydrous_data_dic, new_fluid_volumes
                master_rock[item]['minimization'].thermodynamic_looping_station()
                # Main dictionary save
                # saving G_sys per mol of system
                master_rock[item]['g_sys'].append(master_rock[item]['minimization'].g_sys /
                                                sum(master_rock[item]['minimization'].df_phase_data.iloc[0, :]))
                # g_sys.append(minimization.g_sys)
                master_rock[item]['pot_data'].append(
                    master_rock[item]['minimization'].pot_frame)
                # Backup dictionary save
                rock_origin[item]['g_sys'] = copy.deepcopy(
                    master_rock[item]['g_sys'][-1])
                rock_origin[item]['pot_data'] = copy.deepcopy(
                    master_rock[item]['pot_data'][-1])

                # //////////////////////////////////////////////////////////////////////////
                # 2) Creating DataFrame structure for "df_var_dictionary" in first itteration
                # LINK - 2) Setup dataframes
                if num < 1:
                    # Volume and Density ouput - Dataframes (df_N, df_Vol% etc)
                    for variable in list(master_rock[item]['minimization'].df_phase_data.index):
                        master_rock[item]['df_var_dictionary']['df_' +
                                                            str(variable)] = pd.DataFrame()
                    water_cont_ind = ["N", "H2O[pfu]", "H2O[mol]",
                                    "H2O[g]", "wt%_phase", "wt%_solids", "wt%_H2O.solid"]
                    for variable in water_cont_ind:
                        master_rock[item]['df_h2o_content_dic']['df_' +
                                                                str(variable)] = pd.DataFrame()
                    # Copy to backup dictionary
                    rock_origin[item]['df_var_dictionary'] = copy.deepcopy(
                        master_rock[item]['df_var_dictionary'])
                    rock_origin[item]['df_h2o_content_dic'] = copy.deepcopy(
                        master_rock[item]['df_h2o_content_dic'])
                # updating dictionary with newly calculated data
                master_rock[item]['minimization'].merge_dataframe_dic()
                print("\n")
                print("////// Energy minimization executed //////")
                print("\n")

                # //////////////////////////////////////////////////////////////////////////
                # multi-rock loop for updating data storage, MicaPotassium, SystemFluidTest, init oxygen-isotope module, mineral fractionation
                # //////////////////////////////////////////////////////////////////////////
                # 3)
                # LINK - 3) Data storage & merge
                # calling dictionaries and dataframe for up-to-date usage
                print(f"Running data re-storage, MicaPotassium, SystemFluidTest, oxy-module and mineral fractionation")
                # for item in master_rock:
                master_rock[item]['df_var_dictionary'], master_rock[item]['df_h2o_content_dic'], master_rock[item]['df_element_total'] = (
                    master_rock[item]['minimization'].df_var_dictionary,
                    master_rock[item]['minimization'].df_hydrous_data_dic,
                    master_rock[item]['minimization'].df_all_elements)
                master_rock[item]['df_element_total'] = master_rock[item]['df_element_total'].iloc[:, :-1]
                # hydrogen content of the system before extraction
                if 'H' in master_rock[item]['df_element_total'].index:
                    master_rock[item]['total_hydrogen'] = master_rock[item]['df_element_total']['total:']['H']
                else:
                    master_rock[item]['total_hydrogen'] = 0
                master_rock[item]['st_elements'] = pd.concat(
                    [master_rock[item]['st_elements'], master_rock[item]['df_element_total']['total:']], axis=1)
                # Backup dictionary - merging the data
                for kkey in rock_origin[item]['df_var_dictionary'].keys():
                    cdata = pd.concat(
                        [rock_origin[item]['df_var_dictionary'][kkey], master_rock[item]['df_var_dictionary'][kkey].iloc[:, -1]], axis=1)
                    rock_origin[item]['df_var_dictionary'][kkey] = copy.deepcopy(
                        cdata)
                for kkey in rock_origin[item]['df_h2o_content_dic'].keys():
                    # Testing if dataframe is empty
                    if rock_origin[item]['df_h2o_content_dic'][kkey].empty is True:
                        cdata = rock_origin[item]['df_h2o_content_dic'][kkey]
                    # If not empty, merge data
                    else:
                        cdata = pd.concat(
                            [rock_origin[item]['df_h2o_content_dic'][kkey], master_rock[item]['df_h2o_content_dic'][kkey].iloc[:, -1]], axis=1)
                    rock_origin[item]['df_h2o_content_dic'][kkey] = copy.deepcopy(
                        cdata)
                cdata = master_rock[item]['df_element_total']
                rock_origin[item]['df_element_total'].append(copy.deepcopy(cdata))
                # //////////////////////////////////////////////////////////////////////////
                # store Mica potassium if stable
                for phase in master_rock[item]['df_element_total'].columns:
                    if 'PHNG' in phase:
                        master_rock[item]['mica_K'].append(
                            [temperature[tt], master_rock[item]['df_element_total'][phase]['K']])
                # //////////////////////////////////////////////////////////////////////////
                # 4)
                # LINK - 4) System and fluid volumes
                # Checking for fluid/solid volumes at t = 0 and t = -1,
                # calculating difference (used for escape/extraction rule e.g., factor method)
                # print("-> Testing for aq fluid in the system")
                master_rock[item]['minimization'].step_on_water()
                # Checking fluid and solid volumes. Storing t(-1) and t(0) data before calling calculation
                if num < 1:
                    master_rock[item]['fluid_volume_new'] = master_rock[item]['minimization'].new_fluid_Vol
                    master_rock[item]['fluid_volume_before'] = master_rock[item]['fluid_volume_new']
                    master_rock[item]['solid_volume_new'] = master_rock[item]['minimization'].new_solid_Vol
                    master_rock[item]['solid_volume_before'] = master_rock[item]['solid_volume_new']
                else:
                    master_rock[item]['fluid_volume_before'] = master_rock[item]['minimization'].free_water_before
                    master_rock[item]['fluid_volume_new'] = master_rock[item]['minimization'].new_fluid_Vol
                    if master_rock[item]['minimization'].solid_vol_before != master_rock[item]['st_solid'][-1]:
                        print("\nWARNING: solid volume mismatch\n")
                    master_rock[item]['solid_volume_before'] = master_rock[item]['minimization'].solid_vol_before
                    master_rock[item]['solid_volume_before'] = master_rock[item]['st_solid'][-1]
                    master_rock[item]['solid_volume_new'] = master_rock[item]['minimization'].new_solid_Vol
                # //////////////////////////////////////////////////////////////////////////
                # 5)
                # LINK - 5) Oxygen fractionation module
                # isotope fractionation module
                # print("-> Oxygen isotope module initiated")
                print("d18O before oxygen isotope module")
                print(f"Value is {master_rock[item]['bulk_oxygen']}")
                master_rock[item]['model_oxygen'] = Isotope_calc(
                    master_rock[item]['df_var_dictionary']['df_N'], master_rock[item]['minimization'].sol_sol_base,
                    master_rock[item]['df_element_total'], oxygen_signature=master_rock[item]['bulk_oxygen'])
                # print(f"Bulk oxygen is {master_rock[item]['bulk_oxygen']}")
                # Isotope calc - function for oxygen isotope signatures
                master_rock[item]['model_oxygen'].frac_oxygen(temperature[tt])
                # storing isotope fractionation result, dic in list appended
                master_rock[item]['save_oxygen'].append(
                    master_rock[item]['model_oxygen'].oxygen_dic)
                master_rock[item]['save_bulk_oxygen_pre'].append(
                    master_rock[item]['bulk_oxygen'])
                ### Backup dictionary - save oxygen data
                rock_origin[item]['save_oxygen'].append(copy.deepcopy(master_rock[item]['model_oxygen'].oxygen_dic))
                print("d18O after oxygen isotope module")
                print(f"Value is {master_rock[item]['bulk_oxygen']}")
                # //////////////////////////////////////////////////////////////////////////
                # 6)
                # LINK - 6) Mineral Fractionation
                # mineral (garnet) fractionation - coupled oxygen bulk modification
                if grt_frac == True:
                    if master_rock[item]['database'] == 'ds62mp.txt' or master_rock[item]['database'] == 'td-ds62-mb50-v07.txt':
                        garnet_name = 'GRT'
                    else:
                        garnet_name = 'GARNET'
                    # print("-> Mineral fractionation initiated")
                    # old frac position this line
                    for phase in master_rock[item]['df_element_total'].columns:
                        if '_' in phase:
                            pos = phase.index('_')
                            name = phase[:pos]
                            if name == garnet_name:
                                # Michelles atigorite fractionation
                                # # if name=='GARNET' or name=='SERP' or name=='BR':
                                new_bulk_oxygen = master_rock[item]['minimization'].mineral_fractionation(
                                    master_rock[item]['save_oxygen'][-1], name)
                                master_rock[item]['garnet'].append(
                                    master_rock[item]['minimization'].separate)
                                # Collect the volume of each formed garnet
                                master_rock[item]['meta_grt_volume'].append(master_rock[item]['minimization'].separate.volume)
                                # print(
                                #     f"Selected phase = {phase} with Vol% = {master_rock[item]['df_var_dictionary']['df_vol%'].loc[phase][-2:]}")
                                # print(
                                #     f"Bulk deltaO changed from {round(master_rock[item]['bulk_oxygen'], 3)} to {round(new_bulk_oxygen, 3)}")
                                master_rock[item]['bulk_oxygen'] = new_bulk_oxygen
                                print("_______________________")
                                master_rock[item]['garnet_check'].append(1)
                            if len(master_rock[item]['garnet_check']) < num:
                                master_rock[item]['garnet_check'].append(0)
                # Do not delete - necessary step -
                master_rock[item]['df_element_total'] = master_rock[item]['minimization'].df_all_elements
                # LINK - 7) Metastable garnet
                # Adding the metastable garnet impact
                # calculate the metastable garnet for all bits besides the last
                # add the calculated volume to the solid volume of the current step (this is then saved to the st_solid and used next turn)
                grt_flag = True
                if grt_flag is True and len(master_rock[item]['garnet']) > 0 and len(master_rock[item]['garnet_check']) > 1:
                    if master_rock[item]['garnet_check'][-1] == 0:
                        print("Garnet protocol")
                        # take al modelled garnets
                        # LINK - Metastable garnet call
                        metastable_garnet = Garnet_recalc(self.theriak, master_rock[item]['garnet'], temperature[tt], pressures[num][tt])
                        metastable_garnet.recalculation_of_garnets(database=master_rock[item]['database'], garnet_name=garnet_name)
                        print(f"Fluid volume = {master_rock[item]['fluid_volume_new']} ccm")
                        print(f"Solid volume = {master_rock[item]['solid_volume_new']} ccm")
                        volume = metastable_garnet.recalc_volume
                        # Adding the metastable volume of the fractionated garnet to the current minimized solied volume!
                        master_rock[item]['solid_volume_new'] += volume
                    if len(master_rock[item]['garnet']) > 1 and master_rock[item]['garnet_check'][-1] == 1:
                        print("protocol")
                        # take all garnets but last one
                        metastable_garnet = Garnet_recalc(self.theriak, master_rock[item]['garnet'][:-1], temperature[tt], pressures[num][tt])
                        metastable_garnet.recalculation_of_garnets(database=master_rock[item]['database'], garnet_name=garnet_name)
                        print(f"Fluid volume = {master_rock[item]['fluid_volume_new']} ccm")
                        print(f"Solid volume = {master_rock[item]['solid_volume_new']} ccm")
                        volume = metastable_garnet.recalc_volume
                        # Adding the metastable volume of the fractionated garnet to the current minimized solied volume!
                        master_rock[item]['solid_volume_new'] += volume
                        master_rock[item]['meta_grt_weight'].append(metastable_garnet.recalc_weight)
                else:
                    master_rock[item]['solid_volume_new'] += np.array(master_rock[item]['meta_grt_volume']).sum()
                    if np.array(master_rock[item]['meta_grt_volume']).sum() != np.float64(0):
                        print("What happend now? Here is suddenly garnet which is not stable???")
                        # keyboard.wait('esc')
                    metastable_garnet_weight = 0

                # keeping track of stored and removed fluid
                # !!!! necessary for mechanical module
                master_rock[item]['st_fluid_before'].append(
                    master_rock[item]['fluid_volume_new'])
                master_rock[item]['st_fluid_after'].append(
                    master_rock[item]['fluid_volume_new'])
                master_rock[item]['st_solid'].append(
                    master_rock[item]['solid_volume_new'])

                print("All rocks passed the petrochemical model set-up")
                # //////////////////////////////////////////////////////////////////////////
                # LINK - MECHANICAL FAILURE MODEL
                # 1) Failure and extraction mode selection
                # 2) Oxygen isotope recalculation
                # Physical fracture model section - active if free fluid is present
                # Checking if free fluid is present. Stores this value, initializes
                print(f"Mechanical failure model activated.")
                fluid_name_tag = master_rock[item]['database_fluid_name']
                if fluid_name_tag in list(master_rock[item]['df_element_total'].columns):
                    print("-> Fluid extraction test")
                    # Prepare fluid extraction
                    if master_rock[rock_react_item]['reactivity'].react is True:
                            print("reactivity break")

                    # Prepare fluid extraction
                    # - last entry of st_fluid_after has last fluid volume (number when not extracted or zero when extracted)
                    if len(master_rock[item]['st_fluid_after']) > 1:
                        fluid_before = master_rock[item]['st_fluid_after'][-2]
                    else:
                        fluid_before = master_rock[item]['st_fluid_after'][-1]

                    # Start Extraction Master Module
                    master_rock[item]['fluid_calculation'] = Ext_method_master(
                        pressures[num][tt], temperature,
                        master_rock[item]['df_var_dictionary']['df_volume/mol'].loc[fluid_name_tag][-1],
                        fluid_before, master_rock[item]['fluid_volume_new'],
                        master_rock[item]['solid_volume_before'], master_rock[item]['solid_volume_new'],
                        master_rock[item]['save_factor'], master_rock[item]['master_norm'][-1],
                        master_rock[item]['minimization'].df_phase_data,
                        master_rock[item]['tensile strength'],
                        differential_stress= master_rock[item]['diff. stress'],
                        friction= master_rock[item]['friction'],
                        fluid_pressure_mode= master_rock[item]['fluid_pressure_mode'],
                        fluid_name_tag=fluid_name_tag, subduction_angle=self.angle,
                        rock_item_tag=item
                        )
                    # //////////////////////////////////////////////////////////////////////////
                    # LINK 1) selection of the failure and fluid extraction
                    # ////// Calculation for new whole rock /////////
                    # !!!!Assumption: fracturing of the system
                    # if condition = open system, free water gets extracted
                    fracturing_flag = False # Trigger by default False - active when coulomb module becomes positive
                    failure_mech = master_rock[item]['Extraction scheme']
                    if failure_mech in ['Factor', 'Dynamic', 'Steady', 'Mohr-Coulomb-Permea2', 'Mohr-Coulomb-Griffith']:
                        # ANCHOR work in progres GRIFFITH
                        # if factor_method is True or dynamic_method is True or steady_method is True or coulomb is True or coulomb_permea is True or coulomb_permea2 is True:
                        # LINK Fluid flux & permeability
                        # Fluid flux check - Virtual calculation
                        # Hypothetically momentary fluid flux and permeabiltiy test
                        mü_water = 1e-4
                        # water data
                        v_water = float(master_rock[item]['df_var_dictionary']['df_volume[ccm]'].loc[fluid_name_tag].iloc[-1])
                        d_water = float(master_rock[item]['df_var_dictionary']['df_density[g/ccm]'].loc[fluid_name_tag].iloc[-1])
                        weigth_water = master_rock[item]['df_var_dictionary']['df_wt[g]'].iloc[:, -1][fluid_name_tag]
                        # system data
                        master_rock[item]['meta_grt_weight']
                        v_system = master_rock[item]['solid_volume_new'] + master_rock[item]['fluid_volume_new'] # modelling solid phases + metastable garnet + fluid
                        # add metastable garnet weight - volume is already up to date
                        if len(master_rock[item]['meta_grt_weight']) > 0:
                            weight_sys = float(
                                master_rock[item]['df_var_dictionary']['df_wt[g]'].iloc[:, -1].sum()
                                ) + master_rock[item]['meta_grt_weight'][-1] # weight of solid + fluid + metastable garnet
                        else:
                            weight_sys = float(
                                master_rock[item]['df_var_dictionary']['df_wt[g]'].iloc[:, -1].sum()
                                ) # weight of solid + fluid + metastable garnet

                        d_system = float(master_rock[item]['df_var_dictionary']['df_wt[g]'].iloc[:, -1].sum())/v_system
                        # solid rock data
                        v_rock = master_rock[item]['solid_volume_new']
                        weight_rock = weight_sys - weigth_water # weight of solids + metastable garnet
                        d_rock = weight_rock/v_rock
                        # density difference
                        density_cont = (d_rock-d_water)*1000 # kg/m3

                        # Reading the time interval
                        if num < 1:
                            tc_time_c = float(np.diff(track_time[num:num+2])*365*24*60*60)
                        else:
                            tc_time_c = float(np.diff(track_time[num-1:num+1])*365*24*60*60)


                        """# test 4
                        cubic_meter = 2 # TODO put as input here - geometry or size of the block!
                        xxx = np.power(1_000_000 * cubic_meter, 1/3) # length of x of cube (100cm for 1 m3)
                        size = xxx**3 # cubic size
                        area = xxx**2 # surface size
                        # fluid flux = drainiage flux throughout a column
                        volume_flux = v_water * size /v_system/area/tc_time_c # cm3 cm-2 s-1
                        volume_flux = volume_flux/100 # m3 m-2 s-1
                        # integrated permeability
                        int_permea = volume_flux*mü_water/9.81*xxx/density_cont # permeability in m2"""

                        # test 05
                        bloc_a = np.float64(master_rock[item]['geometry'][0])
                        bloc_b = np.float64(master_rock[item]['geometry'][1])
                        bloc_c = np.float64(master_rock[item]['geometry'][2])
                        area = bloc_b*bloc_c
                        xxx = bloc_a
                        size = bloc_a * bloc_b * bloc_c
                        v_water1 = v_water/1000000 # cm3 to m3
                        v_system1 = v_system/1000000 # cm3 to m3
                        volume_flux = v_water1 * size/v_system1/area/tc_time_c # m3 m-2 s-1
                        int_permea = volume_flux*mü_water/9.81/xxx/density_cont # permeability in m2

                        # udpate to virtual test
                        v_permea = int_permea
                        master_rock[item]['live_fluid-flux'].append(volume_flux)
                        master_rock[item]['live_permeability'].append(int_permea)
                        print(f"-> Virtual permeability test results: {v_permea}")

                        # Latest failure criterion 07.02.2023
                        # LINK i) Mohr-Coulomb failure w. diff. stress input and min. permeabiltiy
                        if failure_mech == 'Mohr-Coulomb-Permea2':
                            print("\t===== Mohr-Couloumb.Permea2 method active =====")
                            master_rock[item]['fluid_calculation'].couloumb_method2(
                                shear_stress=master_rock[item]['shear'],
                                friction=master_rock[item]['friction'],
                                cohesion=master_rock[item]['cohesion']
                                )

                        elif failure_mech == 'Mohr-Coulomb-Griffith':
                            print("\t===== Mohr-Coulomb-Griffith method active =====")
                            if 'diff. stress' in master_rock[item].keys():
                                master_rock[item]['fluid_calculation'].mohr_cloulomb_griffith()
                            else:
                                master_rock[item]['fluid_calculation'].mohr_cloulomb_griffith(
                                        shear_stress=master_rock[item]['shear']
                                        )
                            master_rock[item]['failure module'].append(
                                master_rock[item]['fluid_calculation'].failure_dictionary)

                        # LINK ii) Steady state fluid extraction
                        elif failure_mech == 'Steady':
                            print("===== steady method active =====")
                            master_rock[item]['fluid_calculation'].frac_respo = 5
                            fracturing_flag = True
                            master_rock[item]['fluid_calculation'].fracture = True
                            master_rock[item]['failure module'].append("Steady")

                        else:
                            fracturing_flag = False
                            master_rock[item]['fracture bool'].append(
                                    master_rock[item]['fluid_calculation'].frac_respo)

                        # ##############################################
                        # LINK Coulomb mechanical trigger
                        # Tracking fracturing from coulomb approach methods
                        # Editing trigger
                        # if coulomb is True or coulomb_permea2 is True or coulomb_permea is True:
                        if failure_mech in ['Factor', 'Dynamic', 'Steady', 'Mohr-Coulomb-Permea2', 'Mohr-Coulomb-Griffith']:

                            # store a bool index for the type of fracturing
                            master_rock[item]['fracture bool'].append(
                                master_rock[item]['fluid_calculation'].frac_respo)
                            master_rock[item]['fracture_value'] = 1 + \
                                master_rock[item]['tensile strength'] / \
                                (pressures[num][tt]/10)

                            # Fracture flag trigger
                            # activate fracture flag trigger when vol% of fluid is above or equal 10%
                            if v_water/v_system >= 0.1:
                                master_rock[item]['fluid_calculation'].fracture = True
                                fracturing_flag = True
                                master_rock[item]['fracture bool'][-1] = 10
                            else:
                                fracturing_flag = master_rock[item]['fluid_calculation'].fracture
                            # print(f"\nThe calculated extensional fracturing fator is: .... {fracture_value}\n")
                            # print(f"Check factor: {fluid_fac}")

                            # ##############################################
                            # LINK Release criteria
                            # Fluid Extraction when the modules before give true fracturing
                            # checking with the mohr-coloumb model and decision for fracturing or not
                            """if fracturing_flag is True and v_permea > lowest_permeability[tt]:
                                print("!!! Below minimum permeability!")"""
                            # FIXME modified extraction criteria - minimum permeability is never reached 06.03.2023
                            if fracturing_flag is True:
                                print("Enter fluid extraction")
                                # keyboard activation and exite by esc
                                # if master_rock[rock_react_item]['reactivity'].react is True:
                                #     print("Keyboard wait exception")
                                #     keyboard.wait('esc')
                                master_rock[item]['fluid_extraction'] = Fluid_master(
                                    phase_data=master_rock[item]['minimization'].df_phase_data.loc[:, fluid_name_tag],
                                    ext_data=master_rock[item]['extracted_fluid_data'],
                                    temperature=num+1,
                                    new_fluid_V=master_rock[item]['fluid_volume_new'],
                                    sys_H=master_rock[item]['total_hydrogen'],
                                    element_frame=master_rock[item]['df_element_total'],
                                    st_fluid_post=master_rock[item]['st_fluid_after'],
                                    fluid_name_tag=fluid_name_tag
                                    )
                                # backup before extraction
                                master_rock[item]['fluid_hydrogen'].append(master_rock[item]['df_element_total'][fluid_name_tag]['H'].copy())
                                master_rock[item]['fluid_oxygen'].append(master_rock[item]['df_element_total'][fluid_name_tag]['O'].copy())
                                # Execute the extraction
                                master_rock[item]['fluid_extraction'].hydrogen_ext_all()
                                # Read the data of the fluid extracted
                                master_rock[item]['extracted_fluid_data'] = master_rock[item]['fluid_extraction'].ext_data
                                # Read the element frame when extraction
                                master_rock[item]['df_element_total'] = master_rock[item]['fluid_extraction'].element_frame

                                # Control the rock reactivity
                                master_rock[item]['reactivity'].react=True
                                # Save the time of extraction step
                                master_rock[item]['extr_time'].append(track_time[num])
                                # step system total volume (for surface ---> time integrated fluid flux)
                                master_rock[item]['extr_svol'].append(
                                    np.sum(master_rock[item]['df_var_dictionary']['df_volume[ccm]'].iloc[:, -1]))
                                master_rock[item]['track_refolidv'] = []
                            else:
                                print("!!! No release!")
                                master_rock[item]['reactivity'].react=False
                                master_rock[item]['track_refolidv'].append(
                                    master_rock[item]['solid_volume_before'])
                                master_rock[item]['fracture bool'][-1] = 0

                    # OPTION no extraction
                    # Starts no extraction scheme
                    else:
                        master_rock[item]['reactivity'].react=False
                        print("////// %s No extraction enabled! %s //////")

                    # //////////////////////////////////////////////////////////////////////////
                    # LINK Recalculate the oxygen isotope signature
                    # Recalculate bulk rock oxygen value after possible extraction
                    print("Oxygen isotope signature recalculation fluid extraction - before")
                    print(f"Value is {master_rock[item]['bulk_oxygen']}")
                    if item == 'rock003':
                        print("break")
                    new_O_bulk = oxygen_isotope_recalculation(
                        master_rock[item]['save_oxygen'],
                        master_rock[item]['df_element_total'])
                    # Overwrite for new bulk rock oxygen signature
                    master_rock[item]['bulk_oxygen'] = new_O_bulk
                    # bulk_oxygen = (rockOxy*bulk_oxygen - oxy_mole_fluid *
                    #                 fluid_oxygen)/(rockOxy - oxy_mole_fluid)
                    print("Oxygen isotope signature recalculation fluid extraction - after")
                    print(f"Value is {master_rock[item]['bulk_oxygen']}")
                else:
                    master_rock[item]['save_factor'].append(0)
                    master_rock[item]['fracture bool'].append(0)
                    master_rock[item]['live_fluid-flux'].append(np.nan)
                    master_rock[item]['live_permeability'].append(np.nan)
                    master_rock[item]['failure module'].append("None activated, no fluid.")
                    master_rock[item]['reactivity'].react=False
                    print(f"No free water in the system for {item} - no fracturing model")
                # save bulk oxygen after extraction
                master_rock[item]['save_bulk_oxygen_post'].append(
                    master_rock[item]['bulk_oxygen'])
                print("\n")

            count += 1
            # k += 1
            ic = k/kk*100
            print("=====Progress=====")
            progress(ic)
            print("\n")

        # LINK ROUTINE END
        # //////////////////////////////////////////////////////////////////////////
        # //////////////////////////// END OF SCRIPT ///////////////////////////////
        # //////////////////////////////////////////////////////////////////////////

        self.rock_dic = master_rock

    def data_reduction(self, defined_path=False):
        """
        Perform data reduction on the given parameters. Runs after each of the main routines. Performs the saving of the data to a hdf5 file.

        Args:
            defined_path (bool, optional): Flag indicating whether a defined path is used for saving the data. Defaults to False.

        Returns:
            None
        """

        # Main variables petrology
        master_rock = self.rock_dic
        temperatures = self.temperatures
        pressures = self.pressures
        track_time = self.track_time
        time_steps = 0
        track_depth = self.track_depth
        conv_speed = 0
        angle = 0
        path_methods = self.path_methods
        pathfinder = self.path_methods[-1]

        # Main variables mechanical model
        lowest_permeability = self.minimum_permeability

        # Gridded bulk option
        bulk_2d = False


        # ------------------- Data reduction ----------------------
        # ---> Summarize time information
        for tt, item in enumerate(master_rock):
            # read the fluid name from the database
            fluid_name_tag = master_rock[item]['database_fluid_name']

            if pathfinder is True:
                cum_time = np.array(track_time)
                time_frame = pd.DataFrame({"Time steps": np.ones(
                    len(cum_time))*time_steps, "Cum Time": cum_time})
            else:
                cum_time = np.cumsum(np.array(track_time))
                time_frame = pd.DataFrame(
                    {"Time steps": track_time, "Cum Time": cum_time})

            # shear_stress = master_rock[item]['shear']
            # ---> Initializing module to format data
            # Chronological dataframe (IMPORTANT for prograde + retrograde)
            line = np.arange(1, len(temperatures)+1, 1)
            ref_len = len(master_rock[item]['df_var_dictionary']['df_N'].columns)

            if ref_len == len(line):
                pass
            else:
                line = [0]
                rocks = list(master_rock.keys())
                for val in master_rock[rocks[tt-1]]['extr_time']:
                    line.append(np.where(master_rock['rock0']['time_frame']['Cum Time'] == val)[0][0])
                line = list(line + np.ones(len(line), int))

            # Setting trackable indice for DataFrame
            master_rock[item]['line'] = line

            format_data = System_status(
                master_rock[item]['df_var_dictionary'],
                master_rock[item]['df_h2o_content_dic'],
                master_rock[item]['df_element_total'],
                master_rock[item]['st_elements'])
            format_data.formatting_data(
                temperatures,
                master_rock[item]['st_solid'], master_rock[item]['st_fluid_before'],
                master_rock[item]['st_fluid_after'], master_rock[item]['extracted_fluid_data'],
                line = line)
            master_rock[item]['sys_physicals'] = format_data.sys_dat

            master_rock[item]['time_frame'] = time_frame

            # ---> Creating DataFrame for oxygen fractionation data
            master_rock[item]['all_oxy'] = pd.DataFrame()
            for entry in master_rock[item]['save_oxygen']:
                df = pd.DataFrame(entry['delta_O'], index=entry['Phases'])
                master_rock[item]['all_oxy'] = pd.concat(
                    [master_rock[item]['all_oxy'], df], axis=1)

            master_rock[item]['all_oxy'].columns = line

            # ---> Test for multiple water assignments - origin in Database
            if 'H2O.liq' in master_rock[item]['all_oxy'].index and fluid_name_tag in master_rock[item]['all_oxy'].index:
                water_1 = master_rock[item]['all_oxy'].loc['H2O.liq']
                water_2 = master_rock[item]['all_oxy'].loc[fluid_name_tag]
                water_2[water_2.isna()] = water_1[water_2.isna()]
                master_rock[item]['all_oxy'].loc[fluid_name_tag] = water_2
                master_rock[item]['all_oxy'] = master_rock[item]['all_oxy'].drop(
                    'H2O.liq', axis=0)

            master_rock[item]['pot_data'] = pd.concat(
                master_rock[item]['pot_data'], axis=1)
            master_rock[item]['pot_data'].columns = line

            # testing for empty dataframe - prevent error in later plotting (int-fluid plot)
            if master_rock[item]['extracted_fluid_data'].empty:
                master_rock[item]['extracted_fluid_data'] = np.zeros(
                    len(temperatures))

            # recompile oxygen data
            o_chache = pd.DataFrame()
            for vala in master_rock[item]['save_oxygen']:
                cc = pd.DataFrame(vala['delta_O']).T
                cc.columns = vala['Phases']
                o_chache = pd.concat([o_chache, cc], axis=0)
            o_chache.index = line
            master_rock[item]['save_oxygen'] = o_chache

            # reorganise the trace element data
            # -----------------------------------------------------
            # Initialize an empty list to store the rows
            rows = []
            for (num, pressure, temperature), df in master_rock[item]['trace_element_bulk'].items():
                # Combine the tuple and the DataFrame into a single row
                for index, row in df.iterrows():
                    row_data = [num, pressure, temperature] + row.tolist()
                    rows.append(row_data)
            # Define the column names
            columns = ['num', 'pressure', 'temperature'] + df.columns.tolist()
            # Convert the list to a DataFrame
            master_rock[item]['trace_element_bulk'] = pd.DataFrame(rows, columns=columns)

            # Initialize an empty dictionary to store the DataFrames for each mineral phase
            mineral_phase_dict = {}

            for (num, pressure, temperature), df in master_rock[item]['trace_element_data'].items():
                # print(num)
                # print(df)
                for phase in df.index:
                    if phase not in mineral_phase_dict:
                        mineral_phase_dict[phase] = pd.DataFrame()

                    row_data = [num, pressure, temperature] + list(df.loc[phase])
                    row_df = pd.DataFrame([row_data], columns=['num', 'pressure', 'temperature'] + df.columns.tolist())
                    mineral_phase_dict[phase] = pd.concat([mineral_phase_dict[phase],row_df], axis=0)

            master_rock[item]['trace_element_data'] = mineral_phase_dict

            """ plot to test the data
            norm = np.array([
            0.3670, 0.9570, 0.1370, 0.7110, 0.2310, 0.0870, 0.3060,
            0.0580, 0.3810, 0.0851, 0.2490, 0.0356, 0.2480, 0.0381])

            # get colormap warmcool based on temperatures
            cmap = plt.get_cmap('plasma')
            colors = cmap(np.linspace(0, 1, len(temperatures)))

            for phase in mineral_phase_dict.keys():
                test = mineral_phase_dict[phase].iloc[:,3:]
                fig, ax = plt.subplots()
                for i, row in enumerate(test.iterrows()):
                    index = np.where(temperatures == mineral_phase_dict[phase].iloc[:,2].iloc[i])[0]
                    mineral_phase_dict[phase].iloc[:,2].iloc[0]
                    dat = np.array(row[1])/norm
                    plt.plot(dat, '.-', color=colors[index[0]])
                plt.title(phase)
                plt.xticks(np.arange(len(test.columns)), test.columns, rotation=90)
                plt.yscale('log')
                # add colorbar
                sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=min(temperatures), vmax=max(temperatures)))
                sm.set_array([])  # You can set an e

                # Add the colorbar to the figure
                cbar = fig.colorbar(sm, ax=ax)
                cbar.set_label('Temperature [°C]')
            """

            # -----------------------------------------------------
            # Store some information of the dataframes
            master_rock[item]['phase_order'] = list(master_rock[item][
                'df_var_dictionary']['df_N'].columns)
            master_rock[item]['el_index'] = np.array(
                master_rock[item]['df_element_total'].index)
            master_rock[item]['el_col'] = np.array(
                master_rock[item]['df_element_total'].columns)
            master_rock[item]['pot_tag'] = np.array(
                master_rock[item]['pot_data'].index)
            master_rock[item]['oxy_data_phases'] = list(
                master_rock[item]['all_oxy'].index)

        # LINK Permeability from extractions
        # rock system or rock system?
        """dat = master_rock['rock0']['live_permeability']
        boolean = np.array(master_rock['rock0']['fracture bool'], dtype=bool)
        boolean = np.invert(boolean)
        test = np.ma.masked_array(dat, boolean)
        plt.plot(master_rock['rock0']['depth']/1000,master_rock['rock0']['live_permeability'], 'd--')
        plt.plot(np.ma.masked_array(master_rock['rock0']['depth']/1000, boolean), test, 'd')
        plt.yscale('log')"""


        for rock in master_rock:
            master_rock[rock]['permeability'] = 0
            master_rock[rock]['time-int flux'] = 0

            # translation from live tracking to storing the values
            master_rock[rock]['virtual time-int flux'] = master_rock[rock]['live_fluid-flux']
            master_rock[rock]['virtual permeability'] = master_rock[rock]['live_permeability']

            if isinstance(master_rock[rock]['extracted_fluid_data'], pd.DataFrame) == True:

                # get boolean mask
                boolean = np.array(master_rock[rock]['fracture bool'], dtype=bool)
                boolean = np.invert(boolean)

                # real modelled fluid-fluxes and permeability
                if len(boolean) == len(master_rock[rock]['live_fluid-flux']):
                    master_rock[rock]['time-int flux'] = np.ma.masked_array(master_rock[rock]['live_fluid-flux'], boolean)
                    master_rock[rock]['permeability'] = np.ma.masked_array(master_rock[rock]['live_permeability'], boolean)
                else:
                    master_rock[rock]['time-int flux'] = np.ma.masked_array(master_rock[rock]['virtual time-int flux'], boolean[1:])
                    master_rock[rock]['permeability'] = np.ma.masked_array(master_rock[rock]['virtual permeability'], boolean[1:])

            else:
                pass


        print('\n=====================\
            ===============================\nCalculations fully passed\n====================\
                ================================')

        # universal name for the fluid phase in the databases
        for rock in master_rock:
            # replace "STEAM" in master_rock[rock]["phase_order"] with "fluid"
            master_rock[rock]["phase_order"] = ["fluid" if x==master_rock[rock]["database_fluid_name"] else x for x in master_rock[rock]["phase_order"]]
            master_rock[rock]["el_col"] = ["fluid" if x==master_rock[rock]["database_fluid_name"] else x for x in master_rock[rock]["el_col"]]
            master_rock[rock]["oxy_data_phases"] = ["fluid" if x==master_rock[rock]["database_fluid_name"] else x for x in master_rock[rock]["oxy_data_phases"]]

        for rock in master_rock:
            master_rock[rock]['st_elements_index'] = list(master_rock[rock]['st_elements'].index)

        # //////////////////////////////////////////////////////////////////////////
        # ------------------- Data storing in hdf5----------------------
        no_go = ['minimization', 'model_oxygen',
                'fluid_calculation', 'fluid_extraction', 'reactivity', 'model_tracers']
        meta_h5 = ['bulk', 'new_bulk', 'database',
                'phase_order', 'el_index', 'el_col', 'pot_tag', 'oxy_data_phases', 'database_fluid_name', 'st_elements_index']
        h5_parameters = ['time_frame', 'cohesion', 'temperatures', 'pressures', 'convergence_speed', 'subuction_angle', 'geometry',
                         'shear', 'friction', 'tensile strength', 'Extraction scheme', 'depth', 'diff. stress', 'line', 'theriak_input_record']
        h5_oxygenisotope = ['all_oxy', 'save_oxygen', 'save_bulk_oxygen_pre', 'save_bulk_oxygen_post', 'bulk_oxygen',
                            'bulk_oxygen_after_influx', 'bulk_oxygen_before_influx']
        h5_trace_element = ['trace_element_data', 'trace_element_bulk']
        h5_fluiddata = ['extr_svol', 'extracted_fluid_data', 'fluid_hydrogen', 'fluid_influx_data', 'fluid_oxygen', 'extr_time']
        h5_mechanicaldata = ['fluid_volume_before', 'fluid_volume_new', 'solid_volume_before', 'solid_volume_new', 'save_factor',
                             'st_fluid_before', 'st_fluid_after', 'st_solid', 'fracture bool', 'fracture_value', 'track_refolidv', 'fluid_pressure_mode']
        h5_system_data = [
            'g_sys',
            'df_element_total',
            'st_elements',
            'pot_data',
            ]
        h5_garnet_data = [
            'garnet_check',
            'meta_grt_volume',
            'meta_grt_weight'
        ]

        # ANCHOR - static path
        if defined_path is False:
            # static
            destination = r"C:\Users\Markmann\PhD\Data\03_Proj01_PetroModelling\Angiboust_OszGarnet"
            # transimitting
            # destination = r"C:\Users\Markmann\PhD\Data\03_Proj01_PetroModelling\230210_Condit-test_channeling_stati_transmitting_transmitting-large\transmitting-large\real test"
            # Syros garnet
            # destination = r"C:\Users\Markmann\PhD\Data\05_Grt-metastable_modelling\application to Syros"

            # Selecting path + file-name + file-type
            f_path = destination + f"\{defined_path[:-4]}" + ".hdf5"
        else:
            # when data_reduction comes with init input
            f_path = file_save_path()
            f_path = Path(f_path)



        with h5py.File(f_path, 'w') as hf:

            if bulk_2d is True:
                hf.attrs.create("bulk_grid", bulk_2d)
            else:
                hf.attrs.create("bulk_grid", False)

            for tt, rock in enumerate(master_rock):
                entries = list(master_rock[rock].keys())
                if np.ndim(temperatures) > 1:
                    hf.create_dataset(f"{rock}/Parameters/temperatures", data=np.array(temperatures).T[tt])
                elif np.ndim(temperatures) == 1:
                    hf.create_dataset(f"{rock}/Parameters/temperatures", data=temperatures)
                else:
                    print(f"Temperature array warning. Number of dimensions = {np.ndim(temperatures)} - is not a valid array.".format(temperatures))
                if np.ndim(pressures) > 1:
                    hf.create_dataset(f"{rock}/Parameters/pressures", data=np.array(pressures).T[tt])
                elif np.ndim(pressures) == 1:
                    hf.create_dataset(f"{rock}/Parameters/pressures", data=pressures)
                else:
                    print(f"Pressure array warning. Number of dimensions = {np.ndim(pressures)} - is not a valid array.".format(pressures))
                hf.create_dataset(f"{rock}/Parameters/convergence_speed", data=conv_speed)
                hf.create_dataset(f"{rock}/Parameters/subuction_angle", data=angle)

                for i, item in enumerate(master_rock[rock]):
                    if item in no_go:
                        pass
                    # save master_rock[item]['theriak_input_record'] to hdf5 in rock/item
                    elif item in meta_h5:
                        if item == 'bulk':
                            hf[rock].attrs.create(item, master_rock[rock][item])
                        if item == 'new_bulk':
                            hf[rock].attrs.create(
                                'last_bulk', master_rock[rock][item])
                        if item == 'database':
                            hf[rock].attrs.create(item, master_rock[rock][item])
                        if item == 'phase_order':
                            hf[rock].attrs.create(
                                'Phases', master_rock[rock][item])
                        if item == 'el_index':
                            hf[rock].attrs.create(item, master_rock[rock][item])
                        if item == 'el_col':
                            hf[rock].attrs.create(item, master_rock[rock][item])
                        if item == 'pot_tag':
                            hf[rock].attrs.create(item, master_rock[rock][item])
                        if item == 'oxy_data_phases':
                            hf[rock].attrs.create(item, master_rock[rock][item])
                        if item == 'line':
                            hf[rock].attrs.create(item, master_rock[rock][item])
                        if item == 'st_elements_index':
                            hf[rock].attrs.create(item, master_rock[rock][item])

                    elif item == 'failure module':
                        for subject in master_rock[rock][item]:
                            if isinstance(subject, dict) is True:
                                item_list = list(subject.keys())
                                dataframe_fracture_module = pd.DataFrame(
                                    np.zeros((len(master_rock[rock][item]), len(item_list))))
                                dataframe_fracture_module.columns = item_list
                                for kkey in item_list:
                                    stored = []
                                    for element in master_rock[rock][item]:
                                        if element == 'None activated, no fluid.':
                                            stored.append(np.nan)
                                        else:
                                            stored.append(element[kkey])
                                    dataframe_fracture_module[kkey] = stored

                                # save dataframe_fracture_module to hdf5 in rock/item
                                # dataframe_fracture_module.to_hdf(f_path, f"{rock}/{item}/failure module", mode='a', append=True)
                                hf.create_dataset(f"{rock}/{item}", data=dataframe_fracture_module)
                                hf[f"{rock}/{item}"].attrs.create('header', item_list)

                                break

                    elif item == 'garnet':
                        tts = []
                        pps = []
                        density = []
                        mass = []
                        massp = []
                        moles = []
                        volp = []
                        volume = []
                        name = []
                        elements = pd.DataFrame()
                        volPmole = []

                        for k in range(len(master_rock[rock]['garnet'])):
                            tts.append(master_rock[rock]['garnet'][k].temperature)
                            pps.append(master_rock[rock]['garnet'][k].pressure)
                            density.append(master_rock[rock]['garnet'][k].density)
                            mass.append(master_rock[rock]['garnet'][k].mass)
                            massp.append(master_rock[rock]['garnet'][k].massp)
                            moles.append(master_rock[rock]['garnet'][k].moles)
                            name.append(master_rock[rock]['garnet'][k].name)
                            volp.append(master_rock[rock]['garnet'][k].volp)
                            volume.append(master_rock[rock]['garnet'][k].volume)
                            volPmole.append(master_rock[rock]['garnet'][k].volPmole)
                            frame = pd.DataFrame(master_rock[rock]['garnet'][k].elements[0])
                            frame.index = master_rock[rock]['garnet'][k].elements[1]
                            elements = pd.concat([elements, frame], axis=1)

                        values1 = ['density', 'elements', 'mass', 'massp', 'moles',
                                'name', 'pressure', 'temperature', 'volp', 'volume', 'VolumePMole']
                        values2 = [density, elements, mass, massp,
                                moles, name, pps, tts, volp, volume, volPmole]
                        for j in range(len(values1)):
                            dataset = hf.create_dataset(
                                f"{rock}/GarnetData/{item}/{values1[j]}", data=values2[j])

                        hf[rock].attrs.create(item, list(elements.index))

                    elif isinstance(master_rock[rock][item], dict) is True:

                        # NOTE - New naming convention for the dictionary keys
                        # New naming convention prepared 2023.12.07
                        if str(master_rock[rock][item].keys()) == "dict_keys(['df_N', 'df_volume/mol', 'df_volume[ccm]', 'df_vol%', 'df_wt/mol', 'df_wt[g]', 'df_wt%', 'df_density[g/ccm]'])":
                            new_dic_names = ['df_N', 'df_vol_mol', 'df_volume', 'df_vol%', 'df_wt_mol', 'df_wt', 'df_wt%', 'df_density']
                            old_dic_names = list(master_rock[rock][item].keys())
                            for jj, entry in enumerate(new_dic_names):
                                dataset = hf.create_dataset(f"{rock}/SystemData/{item}/{entry}",
                                                        data=master_rock[rock][item][old_dic_names[jj]])
                        else:
                            for entry in master_rock[rock][item]:
                                dataset = hf.create_dataset(f"{rock}/SystemData/{item}/{entry}",
                                                            data=master_rock[rock][item][entry])

                    # Structuring the single variables in groups for the hdf5 file (paramters, isotope data, fluid data, mechanical data, other recordings)
                    elif item in h5_parameters:
                        if item == 'theriak_input_record':
                            data_in = pd.DataFrame(master_rock[rock][item])
                            data_in_bulk = np.array(data_in['bulk'])
                            hf.create_dataset(f"{rock}/{item}", data=data_in_bulk)
                            # hf[f"{rock}/{item}"].attrs.create('header', list(master_rock[rock][item].keys()))
                        else:
                            hf.create_dataset(f"{rock}/Parameters/{entries[i]}", data=master_rock[rock][item])
                    elif item in h5_oxygenisotope:
                        hf.create_dataset(f"{rock}/IsotopeData/{entries[i]}", data=master_rock[rock][item])
                    elif item in h5_trace_element:
                        if item == 'trace_element_data':
                            for phase in master_rock[rock][item]:
                                hf.create_dataset(f"{rock}/TraceElementData/{phase}", data=master_rock[rock][item][phase])
                                # add row and column names to the dataset
                                hf[f"{rock}/TraceElementData/{phase}"].attrs.create('header', list(master_rock[rock][item][phase].columns))

                        else:
                            hf.create_dataset(f"{rock}/TraceElementBulk/{entries[i]}", data=master_rock[rock][item])
                            # add row and column names to the dataset
                            hf[f"{rock}/TraceElementBulk/{entries[i]}"].attrs.create('header', list(master_rock[rock][item].columns))
                    elif item in h5_fluiddata:
                        hf.create_dataset(f"{rock}/FluidData/{entries[i]}", data=master_rock[rock][item])
                    elif item in h5_mechanicaldata:
                        hf.create_dataset(f"{rock}/MechanicsData/{entries[i]}", data=master_rock[rock][item])
                    elif item in h5_garnet_data:
                        hf.create_dataset(f"{rock}/GarnetData/{entries[i]}", data=master_rock[rock][item])
                    elif item in h5_system_data:
                        hf.create_dataset(f"{rock}/SystemData/{entries[i]}", data=master_rock[rock][item])
                    elif item == 'master_norm':
                        pass
                    else:
                        hf.create_dataset(f"{rock}/Other/{entries[i]}", data=master_rock[rock][item])

        print('\n=====================\
            ===============================\nHDF5 data saved\n====================\
                ================================')
        # copy the init file to data destination folder
        source = Path(defined_path)
        f_path = str(f_path)
        pos = f_path.rfind("\\")
        destination = Path(f_path[:pos+1])
        # test if file already exists
        if os.path.isfile(destination / source.name):
            pass
        else:
            shutil.copy(source, destination)


