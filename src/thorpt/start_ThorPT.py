
"""
Written by
Thorsten Markmann
thorsten.markmann@unibe.ch
status: 16.07.2024
"""

# External modules
import numpy as np
import pandas as pd
import os
from pathlib import Path
import copy
from dataclasses import dataclass
from tkinter import filedialog
import time
import scipy.special

# ThorPT modules
# from thorpt_thorPT.valhalla.Pathfinder import *
# from thorpt_thorPT.valhalla.routines_ThorPT import *
# from valhalla.Pathfinder import *
# from valhalla.routines_ThorPT import *
# from valhalla.tunorrad import run_theriak as test_theriak
# from valhalla import Pathfinder
# from valhalla import routines_ThorPT
# from valhalla.tunorrad import run_theriak as test_theriak


# NOTE: this is the main file to run the ThorPT routine
# import thorpt_thorPT.valhalla.Pathfinder as nasa
# from thorpt_thorPT.valhalla.routines_ThorPT import *
# from thorpt_thorPT.valhalla.tunorrad import run_theriak as test_theriak

def file_opener():
    """
    Opens a file dialog to select an initialization file to read.

    Returns:
        str: The path of the selected file.
    """
    filein = filedialog.askopenfilename(
        title="Select init file to read",
        filetypes=(
            ("txt file", "*.txt"),
            ("All files", "*.*"))
    )

    return filein


def file_opener_multi():
    """
    Opens a file dialog to select multiple files.

    Returns:
        filein (tuple): A tuple containing the paths of the selected files.
    """
    filein = filedialog.askopenfilenames(
        title="Select init file to read",
        filetypes=(
            ("txt file", "*.txt"),
            ("All files", "*.*"))
    )

    return filein

def set_origin():
    """
    Sets the current working directory to the directory of the script.

    This function retrieves the absolute path of the current script file and changes the working directory to that path.
    """
    dirname = os.path.dirname(os.path.abspath(__file__))
    os.chdir(dirname)

# Function to process bulk composition
def process_bulk(entry, rock_init):
    pos = entry.index(":")
    bulkr = entry[pos+1:].split('\t')[-1]
    bulkr = bulkr[1:-1].split(',')
    bulkr = [float(item) for item in bulkr]
    hydrogen_mole = next(float(inner_entry.split(":")[-1].split('\t')[-1]) for inner_entry in rock_init if 'Mole of H' in inner_entry)
    carbon_mole = next(float(inner_entry.split(":")[-1].split('\t')[-1]) for inner_entry in rock_init if 'Mole of C' in inner_entry)
    bulkr.extend([hydrogen_mole, carbon_mole])
    return bulkr


@dataclass
class rockactivity:
    """
    This class represents a rock activity.

    Attributes:
        function: The function associated with the rock activity.
        react: The reaction associated with the rock activity.
    """
    function: any
    react: any



def run_main_routine():

    # Test if function is not called by a script
    if __name__ != '__main__':
        print("Import the modules of ThorPT during package import")
        from thorpt.valhalla import Pathfinder
        from thorpt.valhalla import routines_ThorPT
        from thorpt.valhalla.tunorrad import run_theriak as test_theriak
        print("Import finished")

    set_origin()

    debugging_recorder = []

    # /////////////////////////////////////////////////////
    # Starting up and select the init file to model
    init_interface = True
    init_interface2 = False
    init_files = []

    # Single file selector
    if init_interface is True:
        in_file = file_opener()
        init_files.append(in_file)

    # Multiple file selector
    elif init_interface2 is True:
        in_file = file_opener_multi()
        init_files = list(in_file)

    # Manual selector
    else:
        init_files = ["_initmulti73_.txt", "_initmulti74_.txt", "_initmulti75_.txt", "_initmulti76_.txt", "_initmulti77_.txt", "_initmulti79_.txt"]



    for init_file in init_files:
        print("Script is starting...\u03BA\u03B1\u03BB\u03B7\u03BC\u03B5\u03C1\u03B1!")

        # /////////////////////////////////////////////////////
        # LINK init file reading
        # Starting and reading the init-file
        main_folder = Path(__file__).parent.absolute()
        file_to_open = main_folder / "DataFiles" / init_file

        with open(file_to_open, 'r') as file:
            init = file.read()

        init = init.splitlines()
        init_data = {}

        # take path, fracturing and further settings from init file applied to all rocks modelled
        path_arguments = False
        path_incerement = False
        conditions = {
            'Theriak:': 'theriak',
            'Path:': 'path',
            'Path-increments:': 'path_increment',
            'Extraction scheme:': 'extraction',
            'Marco:': 'Marco',
            'grt_frac_off:': 'grt_frac_off',
            'Input-arguments:': 'path_arguments',
            'Slab-Temperature-Difference:': 'slab_temp_diff'
        }

        for entry in init:
            if '*' in entry:
                continue
            for condition, key in conditions.items():
                if condition in entry:
                    print(f"{condition}:\n{entry}")
                    pos = entry.index(":")
                    value = entry[pos+1:].split('\t')[-1]
                    if key == 'path_arguments':
                        value = value.split('\t')[-1].split(', ')
                        if 'new' in value:
                            pressure_unit = input("You want to digitize a new P-T path. Please provide the physical unit for the pressure value you intend to use.[kbar/GPa]")
                            value[2] = pressure_unit
                    elif key == 'path_increment':
                        value = value.split('\t')[-1].split(', ')
                    elif key == 'min_permeability' or key == 'shearstress':
                        value = float(value)
                    elif key == 'Marco' or key == 'grt_frac_off':
                        value = True
                    elif key == 'slab_temp_diff':
                        value = float(value)
                    init_data[key] = value

        # updating the pressure unit in init file when a new path is created
        for j, entry in enumerate(init):
            if 'Input-arguments' in entry:
                redo_arguments = ', '.join(init_data['path_arguments'])
                redo_arguments = ["Input-arguments:", redo_arguments]
                redo_arguments = ''.join(redo_arguments)
                init[j] = redo_arguments
        # write new condition to init file
        redo_init = '\n'.join(init)
        with open(file_to_open, 'w') as file:
            file.write(redo_init)

        answer = int(init_data['path_arguments'][-1])
        init_data['path_arguments'] = init_data['path_arguments'][:-1]
        # answer = 2

        debugging_recorder.append("First init file block succesfully read.\n")

        breaks = []
        for i, entry in enumerate(init):
            if entry == '***':
                breaks.append(i)
        rock_dic = {}
        for i, item in enumerate(breaks):
            if item == breaks[-1]:
                continue
            rocktag = "rock" + f"{i:03}"
            rock_dic[rocktag] = init[item+1:breaks[i+1]]

        database = []
        bulk = []
        oxygen = []
        init_data['shear'] = []
        init_data['diffstress'] = []
        init_data['friction'] = []
        init_data['geometry'] = []
        init_data['Tensile strength'] = []
        init_data['Extraction scheme'] = []
        init_data['Min Permeability'] = []
        init_data['fluid_name_tag'] = []
        init_data['fluid_pressure'] = []
        init_data['fluid_volume'] = []
        init_data['extraction_percentage'] = []
        init_data['connectivity_fraction'] = []

        # Dictionary to map keys to their corresponding lists in init_data
        key_map = {
            'Database': 'database',
            'Bulk': 'bulk',
            'OxygenVal': 'oxygen',
            'Tensile strength': 'Tensile strength',
            'Geometry': 'geometry',
            'Friction': 'friction',
            'ShearStress': 'shear',
            'Diffential stress': 'diffstress',
            'Extraction scheme': 'Extraction scheme',
            'Minimum Permeability': 'Min Permeability',
            'Fluid phase name': 'fluid_name_tag',
            'Fluid pressure': 'fluid_pressure',
            'Fluid Volume Treshold[Vol%]': 'fluid_volume',
            'Extraction percentage[%]:': 'extraction_percentage',
            'Fluid Connectivity[Vol%]': 'connectivity_fraction'
            }



        # TODO add name from init file?
        for rock in rock_dic.keys():
            rock_init = rock_dic[rock]
            for entry in rock_init:
                if 'Database' in entry:
                    pos = entry.index(":")
                    db = entry[pos+1:].split('\t')[-1]
                    database.append(db + ".txt")
                if 'Bulk' in entry:
                    # read bulk composition and add hyrogen and carbon mole fraction
                    pos = entry.index(":")
                    bulkr = entry[pos+1:].split('\t')[-1]
                    bulkr = bulkr.split('\t')[-1]
                    bulkr = bulkr[1:-1].split(',')
                    for j, item in enumerate(bulkr):
                        bulkr[j] = float(item)
                    for inner_entry in rock_init:
                        if 'Mole of H' in inner_entry:
                            pos = inner_entry.index(":")
                            hydrogen_mole = np.float64(inner_entry[pos+1:].split('\t')[-1])
                        if 'Mole of C' in inner_entry:
                            pos = inner_entry.index(":")
                            carbon_mole = np.float64(inner_entry[pos+1:].split('\t')[-1])
                    bulkr.append(hydrogen_mole)
                    bulkr.append(carbon_mole)
                    bulk.append(bulkr)
                if 'OxygenVal' in entry:
                    pos = entry.index(":")
                    soxygen = entry[pos+1:].split('\t')[-1]
                    soxygen = float(soxygen)
                    oxygen.append(soxygen)
                if 'Tensile strength' in entry:
                    pos = entry.index(":")
                    tensile = entry[pos+1:].split('\t')[-1]
                    init_data['Tensile strength'].append(float(tensile))
                if 'Geometry' in entry:
                    pos = entry.index(":")
                    abc = entry[pos+1:].split('\t')[-1]
                    abc = abc[1:-1].split(',')
                    init_data['geometry'].append(abc)
                if 'Friction' in entry:
                    pos = entry.index(":")
                    friction = entry[pos+1:].split('\t')[-1]
                    init_data['friction'].append(float(friction))
                if 'ShearStress' in entry:
                    pos = entry.index(":")
                    shear = entry[pos+1:].split('\t')[-1]
                    init_data['shear'].append(float(shear))
                if 'Diffential stress' in entry:
                    pos = entry.index(":")
                    diff_stress = entry[pos+1:].split('\t')[-1]
                    init_data['diffstress'].append(float(diff_stress))
                if 'Extraction scheme' in entry:
                    pos = entry.index(":")
                    rock_mechanics = entry[pos+1:].split('\t')[-1]
                    init_data['Extraction scheme'].append(rock_mechanics)
                if 'Minimum Permeability' in entry:
                    pos = entry.index(":")
                    min_permeability = entry[pos+1:].split('\t')[-1]
                    init_data['Min Permeability'].append(float(min_permeability))
                if 'Fluid phase name' in entry:
                    pos = entry.index(":")
                    fluid_name = entry[pos+1:].split('\t')[-1]
                    init_data['fluid_name_tag'].append(fluid_name)
                if 'Fluid pressure' in entry:
                    pos = entry.index(":")
                    fluid_pressure = entry[pos+1:].split('\t')[-1]
                    init_data['fluid_pressure'].append(fluid_pressure.lower())
                if 'Extraction percentage[%]:' in entry:
                    pos = entry.index(":")
                    percentage = entry[pos+1:].split('\t')[-1]
                    init_data['extraction_percentage'].append(float(percentage) / 100)
                if "Fluid Volume Treshold[Vol%]" in entry:
                    pos = entry.index(":")
                    fluid_volume = entry[pos+1:].split('\t')[-1]
                    init_data['fluid_volume'].append(float(fluid_volume))
                if "Fluid Connectivity[Vol%]" in entry:
                    pos = entry.index(":")
                    fluid_volume = entry[pos+1:].split('\t')[-1]
                    init_data['connectivity_fraction'].append(float(fluid_volume))

        init_data['Database'] = database
        init_data['Path'] = init_data['path']
        init_data['Bulk'] = bulk
        init_data['Oxygen'] = oxygen

        debugging_recorder.append("Second init file block succesfully read. All rocks read and stored.\n")

        if __name__ == '__main__':
            print("File call __name__ is set to: {}" .format(__name__))
            from valhalla import Pathfinder
            from valhalla import routines_ThorPT
            from valhalla.tunorrad import run_theriak as test_theriak


        # test run for theriak
        test_output = test_theriak(init_data['theriak'], 'tc55.txt', 500.0, 20000.0, whole_rock="SI(7.9)AL(2.9)FE(0.8)MN(0.0)MG(1.7)CA(1.8)NA(0.7)TI(0.1)K(0.03)H(0.0)C(0.0)O(?)O(0.0)    * CalculatedBulk")
        if len(test_output) > 200:
            print("Theriak test run passed. Theriak is ready to use.")
        else:
            print("Theriak test run failed. Please check the theriak path in the init file.")
            print(test_output)
            print(init_data['theriak'])
            time.sleep(5)
            quit()

        debugging_recorder.append("Theriak test run passed. Theriak is ready to use.\n")

        # /////////////////////////////////////////////////////
        # Preparing input data for modelling routine

        # Set origin to file location
        set_origin()

        # Pre-routine activations from init inputs
        # database = init_data['Database']
        path = init_data['Path']
        rock_bulk = init_data['Bulk']
        oxygen = init_data['Oxygen']
        lowest_permeability = init_data['Min Permeability']

        # Deactivate grt fractionation
        if 'grt_frac_off' in init_data.keys():
            grt_frac = False
        else:
            grt_frac = True

        # /////////////////////////////////////////////////////
        # LINK P-T-t path selection
        # Choose P-T path scheme - True is active
        calc_path = path == "Calculus"
        vho = path == "Vho"
        pathfinder = path == "Finder"
        olivine = path == 'OlivineMod'

        path_arguments = init_data['path_arguments']
        path_increment = init_data['path_increment']

        # P-T- pathway calculation / preparation
        # Calling subduction module
        if calc_path is True:
            print("===== Pathfinder creator active =====")
            # Using pathfinder module to calculate P-T-t path
            function = Pathfinder.Pathfinder_calc(100_000, 100e6, 10, 100)
            function.calc_time_model()
            # Storing P-T-t from pathfinder - modul
            temperatures = np.array(function.T)
            pressures = np.array(function.P)/1e5
            track_time = np.array(function.t_start)

        elif pathfinder is True:
            # Calling Pathfinder module and executing digitization function
            nasa = Pathfinder.Pathfinder()
            # nasa.execute_digi()
            if path_arguments is False and path_increment is False:
                nasa.connect_extern()
            elif path_arguments is False:
                nasa.connect_extern(path_increment)
            elif path_increment is False:
                nasa.connect_extern(path_arguments)
            else:
                nasa.connect_extern(path_arguments=path_arguments, path_increment=path_increment)
            # Store P-T-t information
            temperatures = nasa.temperature
            pressures = nasa.pressure
            track_time = nasa.time
            depth = nasa.depth
            conv_speed = nasa.metadata['Convergence rate [cm/year]']
            angle = nasa.metadata['Burial angle [Degree]']


        elif vho is True:
            print("===== Vho P-T active =====")

            temperatures_start = [350, 375, 400, 450, 500,
                                550, 600, 700, 705, 715, 720, 650, 600, 550, 500]
            pressures_start = [13_000, 14_500, 16_000,
                            18_000, 20_000, 21_500, 23_000, 26_000,
                            25_850, 25_580, 25_000, 23_000, 21_500, 19_000, 18_000]

            # temperatures_start = [350, 375, 400, 450, 500, 550, 600, 700]
            # pressures_start = [13_000, 14_500, 16_000,
            #                    18_000, 20_000, 21_500, 23_000, 26_000]

            # Call pathfinder
            nasa = Pathfinder.call_Pathfinder(temp=temperatures_start,
                                pressure=pressures_start)
            nasa.execute_digi_mod()
            temperatures = nasa.temp
            pressures = nasa.pressure
            track_time = nasa.time_var
            depth = nasa.depth
            pathfinder = True

        else:
            print("No proper P-T-path is selected - no calculations are executed")
            answer2 = input("Do you want to continue? no P-T-path is selected!")
            if answer2.lower().startswith("n"):
                print("Pfff, as if you know the question to the answer 42.....")
                exit()
            if answer2.lower().startswith("y"):
                print("=============================")
                print(
                    "NOOOOOB --- Trying to calculate wihtout a P-T-path...In your dreams!!!")
                exit()

        # /////////////////////////////////////////////////////
        # REVIEW fluid release flag old version - deactivated
        """
        # Choose fluid fractionation scheme - set it to True
        factor_method = False
        steady_method = False
        dynamic_method = False
        coulomb = False
        coulomb_permea = False
        coulomb_permea2 = False
        if extraction == 'No extraction':
            factor_method = False
            steady_method = False
            dynamic_method = False
        if extraction == 'Always extraction':
            steady_method = True
        if extraction == 'Factor method':
            factor_method = True
        if extraction == 'Dynamic routine':
            dynamic_method = True
        if extraction == 'Mohr-Coulomb':
            coulomb = True
        if extraction == 'Mohr-Coulomb-Permea':
            coulomb_permea = True
        if extraction == 'Mohr-Coulomb-Permea2':
            coulomb_permea2 = True"""

        # /////////////////////////////////////////////////////
        # LINK rock directory
        # Setting up the main data directory for on-the-fly storage
        # each rock gets its own tag in master-rock dictionary
        # LINK rock directory
        master_rock = {}
        for i, item in enumerate(rock_bulk):
            tag = 'rock' + f"{i:03}"
            master_rock[tag] = {
                'count': 0,
                'bulk': item,
                'depth': depth,
                'database': database[i],
                'theriak_input_record': False,
                'df_var_dictionary': {},
                'df_h2o_content_dic': {},
                'df_element_total': pd.DataFrame(),
                'g_sys': [],
                'pot_data': [],
                'mica_K': [],
                'geometry': init_data['geometry'][i],
                'master_norm': [],
                'st_fluid_before': [],
                'st_fluid_after': [],
                'st_solid': [],
                'st_elements': pd.DataFrame(),
                'extracted_fluid_data': pd.DataFrame(),
                'fluid_hydrogen': [],
                'fluid_oxygen': [],
                'track_refolidv': [],
                'database_fluid_name': init_data['fluid_name_tag'][i],
                'fluid_pressure_mode': init_data['fluid_pressure'][i],
                'save_oxygen': [],
                'bulk_oxygen': oxygen[i],
                'save_bulk_oxygen_pre': [],
                'save_bulk_oxygen_post': [],
                'bulk_oxygen_before_influx': [],
                'bulk_oxygen_after_influx': [],
                'trace_element_data': {},
                'trace_element_bulk': {},
                'extr_time': [],
                'extr_svol': [],
                'tensile strength': init_data['Tensile strength'][i],
                'fracture bool': [],
                'save_factor': [],
                'friction': init_data['friction'][i],
                'Extraction scheme': init_data['Extraction scheme'][i],
                'failure module': [],
                'fluid_influx_data': pd.DataFrame(),
                'garnet': [],
                'garnet_check': [],
                'meta_grt_volume': [],
                'meta_grt_weight': [],
                'live_fluid-flux': [],
                'live_permeability': [],
                'reactivity': rockactivity(function='base', react=False) if tag == 'rock000' else rockactivity(function='stack', react=False)
            }
            if len(init_data['diffstress']) > 0:
                master_rock[tag]['diff. stress'] = init_data['diffstress'][i]
            if len(init_data['shear']) > 0:
                master_rock[tag]['shear'] = init_data['shear'][i]
            if 'fluid_volume' in init_data.keys():
                master_rock[tag]['extraction treshold'] = init_data['fluid_volume'][i] / 100
            else:
                master_rock[tag]['extraction treshold'] = False
            if 'extraction_percentage' in init_data.keys():
                master_rock[tag]['extraction percentage'] = init_data['extraction_percentage'][i]
            else:
                master_rock[tag]['extraction percentage'] = 1.0
            if 'connectivity_fraction' in init_data.keys():
                master_rock[tag]['fluid connectivity'] = init_data['connectivity_fraction'][i] / 100
            else:
                master_rock[tag]['fluid connectivity'] = 0.0

        # copy pf the master rock dictionary to save all original data before modification while the routine
        rock_origin = copy.deepcopy(master_rock)
        # format bulk rock entry to list to store each iteration
        for rocki in rock_origin.keys():
            entry = rock_origin[rocki]['bulk']
            rock_origin[rocki]['bulk'] = []
            rock_origin[rocki]['bulk'].append(entry)
            rock_origin[rocki]['df_element_total'] = []

        # read time and depth contrains calculated from Pathfinder (only set fpr digi_path, digi and pre_digi)
        if pathfinder is True:
            print("Pathfinder module option is True")
            track_depth = nasa.depth
            time_step = nasa.dt
        else:
            track_depth = [0]
            time_step = 0

        print('\n===================\nScript initialization passed\n======================')

        path_method = (calc_path, vho, pathfinder)

        """mechanical_methods = (
                factor_method,
                steady_method,
                dynamic_method,
                coulomb,
                coulomb_permea,
                coulomb_permea2
                )"""

        debugging_recorder.append("Initializing the routine.\n")

        # Initializing the main module for the routines and select the routine given by the anser in the init file
        if answer == 1:
            debugging_recorder.append("Routine 1 selected.\n")
            # routine 1 simulates every rock individually - no fluid transport between rocks

            """
            # metastable test
            temperatures = np.array([500.0, 550.0])
            pressures = np.array([20000.0, 21500.0])
            track_time = track_time[:2]
            track_depth[17:19]
            """

            ThorPT = routines_ThorPT.ThorPT_Routines(temperatures, pressures, master_rock, rock_origin,
                track_time, track_depth, grt_frac, path_method,
                lowest_permeability, conv_speed, angle, time_step, init_data['theriak'], debugging_recorder)
            ThorPT.unreactive_multi_rock()

        elif answer == 2:
            debugging_recorder.append("Routine 2 selected.\n")
            # routine 2 simulates every rock but allows fluid transport - the P-T path is the same for all rocks
            ThorPT = routines_ThorPT.ThorPT_Routines(temperatures, pressures, master_rock, rock_origin,
                track_time, track_depth, grt_frac, path_method,
                lowest_permeability, conv_speed, angle, time_step, init_data['theriak'], debugging_recorder)
            ThorPT.transmitting_multi_rock()

        elif answer == 3:
            debugging_recorder.append("Routine 3 selected.\n")
            # routine 3 simulates every rock and allows fluid transport - the P-T path is different for all rocks following a layerd scheme
            # layered model and get matrix of temperatures and pressures

            if 'slab_temp_diff' in init_data.keys():
                temperature_matrix, pressure_matrix = Pathfinder.layered_model_PTpatch(temperatures, pressures, init_data['geometry'], init_data['slab_temp_diff'])
            else:
                debugging_recorder.append("No slab temperature difference is set. Default of 100°C is default.\n")
                temperature_matrix, pressure_matrix = Pathfinder.layered_model_PTpatch(temperatures, pressures, init_data['geometry'], 100)

            ThorPT = routines_ThorPT.ThorPT_Routines(temperature_matrix, pressure_matrix, master_rock, rock_origin,
                track_time, track_depth, grt_frac, path_method,
                lowest_permeability, conv_speed, angle, time_step, init_data['theriak'], debugging_recorder)
            ThorPT.transmitting_multi_rock_altPT()

        else:
            print("Script is ending...\u03BA\u03B1\u03BB\u03B7\u03BD\u03C5\u03C7\u03C4\u03B1!")


        if answer in [1, 2, 3]:
            sound_flag = True
            if sound_flag:
                # ANCHOR Sound at end of routine run
                pass
                # playsound(r'C:/Users/Markmann/Downloads/Tequila.mp3')

            # Call the data reduction function
            debugging_recorder.append("Routine finished.\n")
            debugging_recorder.append("Data reduction is starting.\n")
            ThorPT.data_reduction(init_file)
            debugging_recorder.append("Data reduction finished.\n")

    # NOTE playsound to be fixed
    # Play a sound from file location and import a module
    dirname = os.path.dirname(os.path.abspath(__file__))
    new_folder_path = os.path.join(dirname, 'DataFiles')
    new_folder_path = os.path.join(dirname, 'sound')
    sound_01 = os.path.join(new_folder_path, 'wow.mp3')
    sound_02 = os.path.join(new_folder_path, 'Tequila.mp3')


    print("Script is ending...\u03BA\u03B1\u03BB\u03B7\u03BD\u03C5\u03C7\u03C4\u03B1!")
    time.sleep(1)

    """
    import pygame
    # Directing to sounds - play at end
    dirname = os.path.dirname(os.path.abspath(__file__))
    pygame.mixer.init()
    pygame.mixer.music.load(os.path.abspath(f'{dirname}/DataFiles/sound/wow.mp3'))
    pygame.mixer.music.play()
    while pygame.mixer.music.get_busy():
        pygame.time.Clock().tick(10)
    pygame.mixer.music.load(os.path.abspath(f'{dirname}/DataFiles/sound/tequila.mp3'))
    pygame.mixer.music.play()
    while pygame.mixer.music.get_busy():
        pygame.time.Clock().tick(10)
    """ 



if __name__ == '__main__':
    print("File call __name__ is set to: {}" .format(__name__))
    from valhalla import Pathfinder
    from valhalla import routines_ThorPT
    from valhalla.tunorrad import run_theriak as test_theriak

    run_main_routine()

else:
    # The script was imported as a module
    print("Import the modules of ThorPT during package import")
    from thorpt.valhalla import Pathfinder
    from thorpt.valhalla.routines_ThorPT import *
    from thorpt.valhalla.tunorrad import run_theriak as test_theriak
    print("Import finished")






