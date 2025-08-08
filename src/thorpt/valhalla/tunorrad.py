
"""
Written by
Thorsten Markmann
thorsten.markmann@unibe.ch
status: 16.07.2024
"""

# from subprocess import check_output
from pathlib import Path
from collections.abc import Iterable
from matplotlib import pyplot as plt
import pandas as pd
import numpy as np
import scipy
import subprocess
import time
from dataclasses import dataclass, field
from typing import List
from scipy.optimize import minimize
import copy
from collections.abc import Iterable
import keyboard
import os


@dataclass
class Phasecomp:
    """
    Represents a phase composition.

    Attributes:
        name (str): The name of the phase.
        temperature (float): The temperature in degree Celsius.
        pressure (float): The pressure in bar.
        moles (float): The number of moles.
        volume (float): The volume in cubic centimeters (ccm).
        volp (float): The volume percentage.
        mass (float): The mass in grams.
        massp (float): The mass percentage.
        density (float): The density in grams per cubic centimeter (g/ccm).
        elements (any): The elements present in the phase.
        volPmole (float): The volume per mole in cubic centimeters per mole (ccm/mol).
    """
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


@dataclass
class Phaseclass:
    """
    Represents a phase class.

    Attributes:
        comps (List[Phasecomp]): List of phase components.
    """
    comps: List[Phasecomp]

# collision of mohr circle with shear failure envelope


def checkCollision(a, b, c, x, y, radius):
    """
    Check if a point (x, y) is colliding with a line defined by the equation ax + by + c = 0.

    Parameters:
    a (float): The coefficient of x in the line equation.
    b (float): The coefficient of y in the line equation.
    c (float): The constant term in the line equation.
    x (float): The x-coordinate of the point.
    y (float): The y-coordinate of the point.
    radius (float): The radius of the collision detection.

    Returns:
    str: The collision status. Possible values are "Touch" if the point is on the line,
    "Intersect" if the point is inside the collision radius, and "Outside" if the point is outside the collision radius.
    """
    dist = ((abs(a * x + b * y + c)) /
            np.sqrt(a * a + b * b))
    if (radius == dist):
        output = "Touch"
    elif (radius > dist):
        output = "Intersect"
    else:
        output = "Outside"

    return output

def checkCollision_linear(a, b, c, x, y, radius):
    """
    Check if a point (x, y) lies within the collision radius of a linear equation ax + by + c = 0.

    Parameters:
    a (float): The coefficient of x in the linear equation.
    b (float): The coefficient of y in the linear equation.
    c (float): The constant term in the linear equation.
    x (float): The x-coordinate of the point.
    y (float): The y-coordinate of the point.
    radius (float): The collision radius.

    Returns:
    collision (bool): True if the point lies within the collision radius, False otherwise.
    dist (float): The distance of the point from the line defined by the linear equation.
    """
    # Finding the distance of line
    # from center.
    dist = ((abs(a * x + b * y + c)) /
            np.sqrt(a * a + b * b))

    # Checking if the distance is less
    # than, greater than or equal to radius.
    if (radius == dist):
        output = "Touch"
    elif (radius > dist):
        output = "Intersect"
    else:
        output = "Outside"

    if output == "Touch" or output == "Intersect":
        collision = True
    else:
        collision = False

    return collision, dist

def checkCollision_curve(pos, diff, tensile):
    """
    Checks if there is a collision between a circle and a curve.

    Args:
        pos (float): The x-coordinate of the circle center.
        diff (float): The diameter of the circle.
        tensile (float): The tensile strength of the curve.

    Returns:
        tuple: A tuple containing a boolean value indicating whether there is a collision,
               and the minimum distance between the circle center and the curve.
    """

    # Finding the minimum distance of circle center to curve
    # Defining the function for griffith failure
    def f(normal, T=tensile):
        return np.sqrt(4*normal*T +4*T**2)

    P = np.array([pos,0])

    def objective(X):
        X = np.array(X)
        return np.linalg.norm(X - P)

    def c1(X):
        x,y = X
        return f(x) - y

    sol = minimize(objective, x0=[P[0], f(P[0])], constraints={'type': 'eq', 'fun': c1})
    X = sol.x

    minimum = objective(X)
    r = diff/2

    if r >= minimum:
        collision = True
    else:
        collision = False

    return collision, minimum

def check_redo_bulk(bulk):
    """
    Check and update the bulk dictionary with missing elements and their default values.

    Args:
        bulk (dict): The dictionary representing the bulk composition.

    Returns:
        str: The updated bulk composition string.
    """
    el_list = bulk.keys()
    elements = ['H', 'C', 'MN', 'CA', 'TI', 'MG', 'NA', 'K']
    for element in elements:
        if element not in el_list:
            bulk[element] = 0

    new_bulk = (
        f"SI({bulk['SI']})AL({bulk['AL']})FE({bulk['FE']})"
        f"MN({bulk['MN']})MG({bulk['MG']})CA({bulk['CA']})"
        f"NA({bulk['NA']})TI({bulk['TI']})K({bulk['K']})"
        f"H({bulk['H']})C({bulk['C']})O(?)O(0)    * CalculatedBulk"
    )
    return new_bulk


def first_Occurrence_char(char, string):
    """Finds the position of the first occurrence of a character in a string.

    Args:
        char (str): The character to search for.
        string (str): The string to search in.

    Returns:
        int: The position of the first occurrence of the character in the string,
             or 'no' if the character is not found.
    """
    for num, item in enumerate(string):
        if string[num] == char:
            return num
    return 'no'


def remove_items(test_list, item):
    """
    Removes specified item from a list.

    Args:
        test_list (list): A list with items read from a txt file.
        item (str or int): The item to be removed from the list.

    Returns:
        list: A new list with the specified item removed.
    """
    # using list comprehension to perform the task
    res = [i for i in test_list if i != item]
    return res


def flatten(lis):
    """
    Cleaning nested lists to one level.

    Args:
        lis (nested list): Any nested list.

    Yields:
        list: Single level list.
    """
    for item in lis:
        if isinstance(item, Iterable) and not isinstance(item, str):
            for x in flatten(item):
                yield x
        else:
            yield item


def mineral_translation():
    """
    Read DB_database in DataFiles file and stores to list.
    This function translates between theriak and the DBOxygen dataset.

    Returns:
        Dictionary: Including 'DB_phases', 'Theriak_phases', 'SSolution_Check', 'SSolution_set', 'Oxygen_count'
    """

    # Checks for database used and reads the name, to be used for translating file
    with open('XBIN') as f:
        reading = f.readlines()
    database = reading[0].strip()
    name = "SSModel_" + database + ".txt"

    # selecting path and file depending on OS
    main_folder = Path.cwd()
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
    # FIXME - quick fix for spaces at the end of translation file
    # end_line = len(read_translation) - 16
    comp_list = read_translation[mineral_data_line+1:]
    for num, item in enumerate(comp_list):
        if '\n' in item:
            line = item.split('\n')
            line = remove_items(line, '')
            line = line[0].split('\t')
        else:
            line = item.split('\t')
        # function to remove all occurences of '' in list
        line = remove_items(line, '')
        # test for the length of the line - if it is 3, it is ok, if not, we need an additional symbol filter
        if len(line) == 3:
            pass
        else:
            line = remove_items(line, ' ')
            line = remove_items(line, '   ')
            line = remove_items(line, '\n')
            line = remove_items(line, '\t')
        # Store in comprehensive list
        # read first entry and erase whitespaces in the first entry
        line[0] = line[0].replace(' ', '')
        phases_DB.append(line[0])
        # read second entry and erase whitespaces in the second entry
        line[1] = line[1].replace(' ', '')
        phases_Ther.append(line[1])
        # read third entry
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


def decode_lines(line_number, line_file, number=11, one_element_row=True):
    """
    Decode lines from a file and return a dictionary.

    Args:
        line_number (int): The line number to decode.
        line_file (list): The list of lines from the file.
        number (int, optional): The number of elements in each line. Defaults to 11.
        one_element_row (bool, optional): Whether each line has only one element. Defaults to True.

    Returns:
        dict: A dictionary containing the decoded line.

    """
    if one_element_row is False:
        line_number = line_number + 1

    temp_list = line_file[line_number]
    list_decrypt = []
    temp_dictionary = {}

    temp_list = temp_list.split(' ')
    for word in temp_list:
        list_decrypt.append(word)

    while '' in list_decrypt:
        list_decrypt.remove('')

    if one_element_row is False:
        second_row = line_file[line_number+1]
        second_row = second_row.split(' ')
        while '' in second_row:
            second_row.remove('')

        for word in second_row:
            list_decrypt.append(word)

    while '' in list_decrypt:
        list_decrypt.remove('')

    temp_dictionary[list_decrypt[0]] = list_decrypt[1:]
    return temp_dictionary


def run_theriak(theriak_path, database, temperature, pressure, whole_rock, theriak_input_rock_before=False):
    """
    Function to run theriak with specified P-T condition and returns output as a list.
    it includes the path where theriak can be executed, writes the Therin file for the P-T
    conditions and uses a specific chemical composition (NannoOoze - Plank 2014) - date: february, 8 2021)

    Args:
        theriak_path (str): The path where theriak can be executed.
        database (str): The chemical composition database to be used by theriak.
        temperature (int): The value for the temperature condition.
        pressure (int): The value for the pressure condition.
        whole_rock (str): The whole rock composition.

    Returns:
        list: The output from theriak as a list of lines.

    Raises:
        FileNotFoundError: If the theriak executable file is not found.

    """

    ######################################################################
    # New way of calling theriak - now with input by user in init
    # initializes the list were the theriak output is stored
    therin_condition = '    ' + str(temperature) + '    ' + str(pressure)
    file_to_open = Path(theriak_path)
    whole_rock_write = "1   " + whole_rock
    # stores the momentary P, T condition passed to Theriak for calculation
    with open('THERIN', 'w') as file_object:
        file_object.write(therin_condition)
    # opens THERIN and writes new P,T condition
    # with open('THERIN', 'a') as file_object:
        file_object.write("\n")
        file_object.write(whole_rock_write)
    with open('XBIN', 'w') as file_object:
        file_object.write(database)
        file_object.write("\n")
        file_object.write("no")

    ######################################################################
    # Runs Theriak, saves output, strips it to list
    ######################################################################
    # # opens THERIN and writes more input parameters as elemental composition
    # Option 1 - Philips approach
    theriak_input_rock_before = False
    if theriak_input_rock_before == False:
        # Executing minimization for new bulk, P, T condition
        theriak_xbin_in = database + "\n" + "no\n"
        theriak_exe = file_to_open / "theriak"
        out = subprocess.run([theriak_exe],
                                input=theriak_xbin_in,
                                encoding="utf-8",
                                capture_output=True)

        theriak_output = out.stdout
        theriak_output = theriak_output.splitlines()

    elif theriak_input_rock_before['bulk'][-1] == whole_rock and theriak_input_rock_before['temperature'][-1] == temperature and theriak_input_rock_before['pressure'][-1] == pressure:
        with open('ThkRun.log', 'rb') as file:
            # read the txt file and write it into the list
            output = file.read()

        _thk_run_check = output.decode('utf-8').splitlines()
        if 'GARNET' in _thk_run_check[0]:
            # Executing minimization for new bulk, P, T condition
            theriak_xbin_in = database + "\n" + "no\n"
            theriak_exe = file_to_open / "theriak"
            out = subprocess.run([theriak_exe],
                                    input=theriak_xbin_in,
                                    encoding="utf-8",
                                    capture_output=True)

            theriak_output = out.stdout
            theriak_output = theriak_output.splitlines()

        else:
            print("no minimization - txt-read only")
            # this argument runs the copy script - avoiding running the same minimization twice
            # reading the already written theriakoutput
            # does not call the process but only reads ThkOut.txt into the 'theriak_in_lines'-list
            output = 1
            with open('ThkOut.txt', 'rb') as file:
                # read the txt file and write it into the list
                output = file.read()

            theriak_output = output.decode('utf-8').splitlines()

    else:
        # Executing minimization for new bulk, P, T condition
        theriak_xbin_in = database + "\n" + "no\n"
        theriak_exe = file_to_open / "theriak"
        out = subprocess.run([theriak_exe],
                                input=theriak_xbin_in,
                                encoding="utf-8",
                                capture_output=True)

        theriak_output = out.stdout
        theriak_output = theriak_output.splitlines()

    ####################################
    # Option 2 - Old approach
    """cmd = subprocess.Popen([file_to_open, 'XBIN', 'THERIN'],
                           stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    cmdout, cmderr = cmd.communicate()
    if cmderr != b'':
        print("CRITICAL ERROR:\nCmd out error - test input for theriak in thermo loop")
        print(cmderr)
        time.sleep(10)
    output = cmdout
    # t1 = time.time()
    # print(f"The time of option 2 is: {t1-t0}")

    theriak_in_lines = output.decode('utf-8').splitlines()"""


    
    return theriak_output


def read_to_dataframe(index_entry, df_index_as_list, input_in_lines, skip=2, number=200):
    """
    A file read from lines to a list is transformed to a DataFrame.

    Args:
        index_entry (int): Index number in the list that has been read from a txt file and checked for a specific entry.
        df_index_as_list (list): The header for the DataFrame (moles etc from theriak output).
        input_in_lines (list): The file from the txt file that has been read to a list by its lines.
        skip (int, optional): The number of lines to skip before the aimed entry is encountered. Defaults to 2.
        number (int, optional): The number of entries to process. Defaults to 200.

    Returns:
        DataFrame: A DataFrame with the defined headers and the values for each stable phase at the current P-T condition.
    """

    data = {}
    current_number = 0
    while current_number <= 12:
        try:
            entry_decoded = decode_lines(
                index_entry + skip + current_number, input_in_lines, number=number
            )
        # calls the decode line function that reads lines into object dictionaries
        # as long as the while loop runs
        except IndexError:
            break
        for item in entry_decoded.values():
            while '|' in item:
                item.remove('|')
            # removing separator
        for key in entry_decoded.keys():
            if key == '----------' or key == 'of' or key == 'total' or key == '' or key == '--------':
                break
            else:
                # NOTE - NaN bug fix substitute by 0.0 - working? Mistake?
                if 'NaN' in entry_decoded[key]:
                    entry_decoded[key] = [0.0 if value == 'NaN' else value for value in entry_decoded[key]]
                entries = pd.to_numeric(entry_decoded[key])
                # turns objects into numbers
                temp_series = pd.Series(entries)
                # stores value from dic key into series
                data[key] = temp_series
        current_number += 1
    df_data = pd.concat(data, axis=1)
    # writes the created 'data' dictionary which contains all
    # stable phases of solids into a DataFrame
    if len(df_data.index) == len(df_index_as_list):
        df_data.index = df_index_as_list
        return df_data
    else:
        pass


def read_theriak(theriak_path, database, temperature, pressure, whole_rock, theriak_input_rock_before):
    """
    Starts and reads the theriak minimization at specific P-T and bulk rock to Dataframes that
    contain the physical values for each stable phase,
    elemental distribution, solid solutions, and phases with water content.

    Args:
        theriak_path (str): The path to the theriak executable.
        database (str): The path to the thermodynamic database.
        temperature (int): The temperature at the given condition.
        pressure (int): The pressure at the given condition.
        whole_rock (list of floats): The calculated bulk rock for the system.

    Returns:
        dict: A dictionary that contains all the created dataframes.
    """

    # Call the connection to the console and execute theriak
    Data_dic = {}
    del Data_dic
    Data_dic = {}
    theriak_in_lines = run_theriak(theriak_path,
        database, temperature, pressure, whole_rock, theriak_input_rock_before)

    # #######################################################
    # Following reads passages from the theriak output
    # #######################################################


    # Snippet cutting method
    # ANCHOR time start Cavalaire
    TIL = theriak_in_lines.copy()
    # first important keyword to shorten list
    keyword = ' equilibrium assemblage:'
    snip_index = TIL.index(keyword)
    TIL = TIL[snip_index:]

    """
    # 2) Volume and densities from stable phases #####
    keyword = ' volumes and densities of stable phases:'
    keyword = (
        '  gases and fluids       N       '
        'volume/mol  volume[ccm]               wt/mol       '
        'wt [g]              density [g/ccm]'
    )
    # 3) H2O content in stable phases
    keyword = ' H2O content of stable phases:'

    # 4) Elements in stable phases #####
    #####################################
    keyword = ' elements in stable phases:'

    ####################################
    # 5) reading chemical potentials
    keyword = ' chemical potentials of components:'
    """

    # 1) Volume and densities from stable phases #####
    ###################################################
    data = {}
    fluid_content_solids_index = 0
    # read volume and densitiy line - theriak - stores values for fluids and solids
    keyword = ' volumes and densities of stable phases:'
    volume_density_index = int(TIL.index(keyword))
    # keyword_search function used for searching the line of index
    volume_density_values_solids = volume_density_index + 5
    # row index of solids and fluids in theriak output
    keyword = (
        '  gases and fluids       N       '
        'volume/mol  volume[ccm]               wt/mol       '
        'wt [g]              density [g/ccm]'
    )
    # try:
    if keyword in TIL:
        volume_density_values_fluids = int(TIL.index(keyword))
        gases_fluids = True
    # except ValueError:
    else:
        volume_density_values_fluids = 1
        gases_fluids = False
        del volume_density_values_fluids
        # print("=1= No gases and fluids in the system detected.")
    # checks if gases and fluids were stable in the system
    fluid_check = (
        '  gases and fluids       N       '
        'volume/mol  volume[ccm]               wt/mol       '
        'wt [g]              density [g/ccm]'
    )
    # theriak entry if fluids are stable - dont edit
    current_number = 0

    try:
        while current_number <= 10:
            try:
                vol_dens_solid = decode_lines(
                    volume_density_values_solids + current_number, TIL, number=200
                )
            except IndexError:
                break
            # decodes line specified by 'volume_density_values_solids' in variable 'theriak_in_lines'
            # calls the decode line function that reads lines into object dictionaries
            # as long as the while loop runs
            for item in vol_dens_solid.values():
                while '|' in item:
                    item.remove('|')
                    # removes symbol in dic
            for key in vol_dens_solid.keys():
                try:
                    if key == '----------' or key == 'of' or key == 'total':
                        break
                    else:
                        entries = pd.to_numeric(vol_dens_solid[key])
                        # converts object from list to int/float
                        temp_series = pd.Series(entries)
                        # stores value from dic key into series
                        entries = pd.to_numeric(vol_dens_solid[key])
                        # converts object from list to int/float
                        temp_series = pd.Series(entries)
                        # stores value from dic key into series
                        data[key] = temp_series
                except ValueError:
                    break
            current_number += 1
    except NameError:
        print("---ERROR:---What? No solid stable phases?")

    # test if data dictionary is empty
    if not data:
        df_Vol_Dens = pd.DataFrame(data)
        # print("---ERROR:---No solid stable phases detected.")
    else:
        df_Vol_Dens = pd.concat(data, axis=1)
    # writes the created 'data' dictionary which contains all
    # stable phases of solids into a DataFrame
    try:
        df_Vol_Dens.index = [
            'N', 'volume/mol', 'volume[ccm]', 'vol%',
            'wt/mol', 'wt[g]', 'wt%', 'density[g/ccm]'
        ]
    except ValueError:
        print("---ERROR:---Value Error with Volume and Density index")

    try:
        TIL[volume_density_values_fluids] == fluid_check
        index_list = ['N', 'volume/mol', 'volume[ccm]',
                      'wt/mol', 'wt[g]', 'density[g/ccm]']
        df_Vol_Dens_fluids = read_to_dataframe(
            index_entry=volume_density_values_fluids, df_index_as_list=index_list, input_in_lines=TIL)
        # decodes output theriakouput saved as a list at given index, index=Volume and densities from stable phases
        df_Vol_Dens = pd.concat([df_Vol_Dens, df_Vol_Dens_fluids], axis=1)
        # merges Dataframes
    except NameError:
        pass

    # 2) H2O content in stable phases #####
    ########################################
    # solid phases
    # keyword = 'H2O content of stable phases:'   ---- error?
    keyword = ' H2O content of stable phases:'
    try:
        fluid_content_solids_index = int(TIL.index(keyword)) + 3
        hydrous_solid = True
        # print("=2== Solids that bound water detected")
    except ValueError:
        # print("No H2O content in solid phases")
        hydrous_solid = False
        del fluid_content_solids_index

    try:
        index_list = ['N', 'H2O[pfu]', 'H2O[mol]', 'H2O[g]',
                      'wt%_phase', 'wt%_solids', 'wt%_H2O.solid']
        df_H2O_content_solids = read_to_dataframe(
            index_entry=fluid_content_solids_index, df_index_as_list=index_list, input_in_lines=TIL)
        # decodes output theriakouput saved as a list at given index, index=H2O content of solid phases
        test_hydsolid = True
        # old test if hydrous minerals are present because of oxxygen isotope issue
        # if temperature > 533:
        #     print(df_H2O_content_solids)
    except NameError:
        test_hydsolid = False
        pass
    # gases and fluids
    keyword = '  gases and fluids       N      H2O[pfu]    H2O[mol]     H2O [g]  |   phase'
    try:
        h2o_content_fluids_index = int(TIL.index(keyword))
        # print("Index fluid - true")
    except ValueError:
        keyword = '  solid phases           N      H2O[pfu]    H2O[mol]     H2O [g]  |   phase'
        try:
            h2o_content_fluids_index = int(TIL.index(keyword))
            test_fluidfluid = True
            print("???solids with water detected")
        except ValueError:
            test_fluidfluid = False
            h2o_content_fluids_index = 1
            # print("optional path fluid not true")
            del h2o_content_fluids_index

    try:
        index_list = ['N', 'H2O[pfu]', 'H2O[mol]', 'H2O[g]', 'wt%_phase']
        df_H2O_content_fluids = read_to_dataframe(
            index_entry=h2o_content_fluids_index, df_index_as_list=index_list, input_in_lines=TIL)
    except NameError:
        df_H2O_content_fluids = pd.DataFrame()
        pass
    if gases_fluids is True:
        if test_hydsolid is False and test_fluidfluid is False:
            pass
        else:
            df_h2o_content = pd.concat(
                [df_H2O_content_solids, df_H2O_content_fluids], axis=1)
        # print("Fluid+solid h2o cont true")
    elif hydrous_solid is True:
        df_h2o_content = pd.concat(
            [df_H2O_content_solids], axis=1)
        # print("Only solid h2o cont true")
    else:
        # print("---EXCEPTION--- no h2o content true")
        pass

    # 3) Elements in stable phases #####
    #####################################

    keyword = ' elements in stable phases:'
    try:
        elements_stable_index = int(TIL.index(keyword))
    except ValueError:
        print("---ERROR:---Where are my elements???")
    # need to define index list by decoding the output first

    temp = TIL[elements_stable_index +
                            3] + TIL[elements_stable_index+4]

    element_list = temp.split()
    one_element_row = True
    element_list = element_list[:element_list.index('E')+1]
    if len(element_list) > 10:
        # print(f"More than one element row - lenght: {len(element_list)}")
        one_element_row = False
    data = {}
    current_number = 0

    TIL_elements = TIL[elements_stable_index:int(TIL.index(' elements per formula unit:'))-1]
    # Read the TIL lines for each phase when more than one element row is present (>10 elements)
    if one_element_row is False:
        # Read every phase from TIL_elements to a temporaray dictionary
        # Each phase has two lines of entry
        num_entry_lines = 2
        # first 5 rows are not needed
        start_index = 5
        TIL_elements = TIL_elements[start_index:]
        phase_temp = []
        # put every two lines into one entry of a list
        i = 0
        while i < len(TIL_elements):
            phase_temp.append(TIL_elements[i:i+num_entry_lines])
            i += num_entry_lines
        # read each item from phase_temp and put it into the data dictionary
        for item in phase_temp:
            tempdata = item[0].split() + item[1].split()
            data[tempdata[0]] = tempdata[1:]

        # test if entry in tempdata[1:] is longer than 8 digits. If yes, split it into two entries
        for key in data.keys():
            tempdata = data[key]
            tempdata_new = []
            for entry in tempdata:
                    # test if length of entry is longer than 8 digits and if two "." are present
                    if len(entry) > 8 and entry.count(".") == 2 :
                        tempdata_new.append(entry[:8])
                        tempdata_new.append(entry[8:])
                    else:
                        tempdata_new.append(entry)
            # convert each entry in tempdata to numeric
            tempdata_new = [float(num) for num in tempdata_new]
            data[key] = tempdata_new

        # write the dictionary to a dataframe
        df_elements_in_phases = pd.DataFrame(data)
        df_elements_in_phases.index = element_list

        # old way of reading the TIL lines for each phase when more than one element row is present (>10 elements)
        """while True:
            if TIL[elements_stable_index + 4 + current_number] == ' elements per formula unit:':
                break
            elements_in_phases = decode_lines(
                elements_stable_index + 4 + current_number,
                TIL,
                number=10, one_element_row=one_element_row
            )
            print(elements_in_phases)
            for key in elements_in_phases.keys():
                if key == ' elements per formula unit:':
                    break
                else:
                    try:
                        entries = pd.to_numeric(elements_in_phases[key])
                        temp_series = pd.Series(entries)
                        data[key] = temp_series
                    except ValueError:
                        pass
                    # print(data[key])
            current_number += 1
            if one_element_row is False:
                current_number += 1"""
    # Read the TIL lines for each phase when only one element row is present (<=10 elements)
    if one_element_row is True:
        count = 0
        while True:
            if 'total' in TIL[elements_stable_index+4+count]:
                el_phase = TIL[elements_stable_index +
                                            4+count].split()
                for i, item in enumerate(el_phase[1:]):
                    el_phase[i+1] = float(item)
                data[el_phase[0]] = pd.Series(el_phase[1:])
                break
            el_phase = TIL[elements_stable_index+4+count].split()
            for i, item in enumerate(el_phase[1:]):
                el_phase[i+1] = float(item)
            data[el_phase[0]] = pd.Series(el_phase[1:])
            count += 1

        df_elements_in_phases = pd.concat(data, axis=1)
        df_elements_in_phases.index = element_list

    Data_dic['df_Vol_Dens'] = df_Vol_Dens
    try:
        Data_dic['df_h2o_content'] = df_h2o_content
    except UnboundLocalError:
        pass
    Data_dic['df_elements_in_phases'] = df_elements_in_phases

    ####################################
    # 4) reading equilibrium assemblage
    keyword = ' equilibrium assemblage:'
    eq_ass_index = int(TIL.index(keyword)) + 12
    stop = volume_density_index
    phase_names = []
    endmember_names = []
    eq_dic = {}
    phase_pos = []
    phase_mol = []
    phase_molper = []
    solid_S = []
    elem_comp = []

    numb_p = int(TIL[int(
        TIL.index(keyword)) + 5].split()[2])
    temp2 = TIL[TIL.index(
        keyword) + 7].split('     ')
    # read the free gibbs system energy
    try:
        g_sys = float(temp2[0].split()[-1])
    except:
        if temp2[1] == '':
            g_sys = float(temp2[0].split()[2])
        elif '*' in temp2[0].split()[-1]:
            print(" * in system energy - failure?")
            g_sys = np.nan
        else:
            g_sys = float(temp2[1])

    while eq_ass_index < stop:

        temp = TIL[eq_ass_index].split()

        if len(temp) > 4:
            try:
                # print(temp)
                float(temp[0])
                float(temp[1])
                res = True

            except ValueError:
                '''entries are not a float'''
                res = False
            if res is True:
                phase_pos.append(eq_ass_index)

        eq_ass_index = eq_ass_index + 1

    for num, item in enumerate(phase_pos):
        # iterating through phase_pos items which are the phases with or without endmembers
        # test item line, if its longer than 5 it has endmembers
        temp = TIL[item].split()
        if temp[2] == 'GARNET':
            break
        if '**' in temp:
            temp.remove('**')

        phase_names.append(temp[2])
        phase_mol.append(temp[3])
        phase_molper.append(temp[4])

        if len(temp) > 5:
            solid_S.append(True)
            endmember_names.append(phase_names[-1])

            elem_comp.append(temp[6])

            check_solid_solutions = True
            phases_endm_temp = []
            elem_comp_temp = []
            endmember = 0
            while check_solid_solutions is True:
                endmember += 1
                # print(TIL[item+endmember])
                ass_cache = TIL[item+endmember].split()
                if '**' in ass_cache:
                    ass_cache.remove('**')
                if not ass_cache:
                    break
                if len(ass_cache) != 5:
                    break

                phases_endm_temp.append(ass_cache[0])
                elem_comp_temp.append(ass_cache[1])

            phases_endm_temp.insert(0, temp[5])
            elem_comp_temp.insert(0, elem_comp[-1])
            for i, ele in enumerate(elem_comp_temp):
                if '*' in ele:
                    elem_comp_temp[i] = 0.0000001
            # old backup
            # elem_comp_temp = [float(ele) for ele in elem_comp_temp]
            # new reading because of odd '00-0.9' value
            # REVIEW odd shit value
            ttemp = []
            for ele in elem_comp_temp:
                if isinstance(ele, float):
                    ele = str(ele)
                else:
                    pass
                if ele[:3] == '00-':
                    ttemp.append(float(ele[2:]))
                else:
                    ttemp.append(float(ele))
            elem_comp_temp = ttemp

            endmember_names[num] = phases_endm_temp
            elem_comp[num] = elem_comp_temp

        else:
            solid_S.append(False)
            endmember_names.append([])

    ####################################
    # 5) reading chemical potentials
    keyword = ' chemical potentials of components:'
    i_base = int(TIL.index(keyword))
    n_elem = True
    pot_index = i_base + 12

    # get number of potential variables in output to define important lines
    i = 0
    while n_elem is True:
        test = TIL[pot_index+i]
        if test == ' ':
            n_elem = i
        i += 1

    pot_col = TIL[pot_index].split()
    pot_d = TIL[pot_index+2:pot_index+n_elem]
    pot_el = []
    pot_val = []
    for i, item in enumerate(pot_d):
        item = item.split()
        pot_el.append(item[0].replace('"', ''))
        pot_val.append(float(item[1]))

    pot_frame = pd.DataFrame(pot_val, index=pot_el, columns=[temperature])

    # ===================================================================
    solid_solution_dic = {}
    solid_solution_dic = {'Name': phase_names, 'Moles': phase_mol, 'Molper': phase_molper,
                          'Memb_Names': endmember_names, 'Comp': elem_comp}
    eq_dic['Names'] = phase_names

    Data_dic['solid_solution_dic'] = solid_solution_dic

    # ANCHOR time stop Cavalaire
    # store iteration end timestamp
    end = time.time()
    # show time of execution per iteration
    # print(f"Script Cavalaire - Time taken: {(end-start)*10**3:.03f}ms")

    return Data_dic, g_sys, pot_frame


def boron_fraction(fluidV, rockV, conc_TE_rock, Frac_fac_whitemica=1.4, conc_fluid=300, open=True):
    """
    First try to implement a boron fractionation between mica and water
    - a lot to do here


    Args:
        fluidV ([type]): [description]
        rockV ([type]): [description]
        conc_TE_rock ([type]): [description]
        Frac_fac_whitemica (float, optional): [description]. Defaults to 1.4.
        conc_fluid (int, optional): [description]. Defaults to 300.
        open (bool, optional): [description]. Defaults to True.

    Returns:
        [type]: [description]
    """
    Ratio_FR = fluidV/rockV * np.arange(1, 100, 10)
    conc_TE_rock = conc_fluid/Frac_fac_whitemica - \
        (conc_fluid/Frac_fac_whitemica - conc_TE_rock) * \
        np.exp(-Ratio_FR*Frac_fac_whitemica)

    return conc_TE_rock, Ratio_FR


def read_trace_element_content():
    # get script file path
    # Get the script file location
    script_path = os.path.dirname(os.path.abspath(__file__))
    # Enter DataFiles folder
    new_folder_path = os.path.join(script_path, 'DataFiles')
    # read txt file distribution_coeff_tracers.txt with pandas
    file_to_open = os.path.join(new_folder_path, 'distribution_coeff_tracers.txt')
    
    # read the txt file with pandas, only the first 7 rows are needed, space as separator
    df = pd.read_csv(file_to_open, sep=', ', nrows=6, engine='python')
    df.index = df['Ratio']
    df = df.iloc[:,2:]

    return df

    


def garnet_bulk_from_dataframe(frame, moles, volumePermole, volume):
    """
    Generate a bulk formula string for garnet from a dataframe.

    Args:
        frame (pandas.DataFrame): The dataframe containing the element concentrations.
        moles (float): The number of moles of the bulk.
        volumePermole (float): The volume per mole.
        volume (float): The total volume.

    Returns:
        str: The bulk formula string for garnet.
    """
    new_bulk1 = []

    # normalization to 1 mol
    rec_moles = volume/volumePermole
    frame = frame * 1/rec_moles
    for el in frame.index:
        if el == 'E':
            pass
        else:
            element_val_str = frame.loc[el][0]
            element_val_str = '{:f}'.format(element_val_str)
            new_bulk1.append(el+'('+element_val_str+')')
            # old version
            # new_bulk1.append(el+'('+str(float(frame.loc[el]))+')')
    new_bulk = ''.join(new_bulk1) + "    " + "* grt extract bulk"

    return new_bulk

def recalculate_FeWtPerc_speciation(total_feo_value, x_fe3):
    oxide_molar_weight = [60.08, 79.87, 101.96, 71.85, 159.69, 70.94, 40.30, 56.08, 61.98, 94.20] # g/mol
    cation_molar_weight = [28.08, 47.87, 26.98, 55.85, 55.85, 54.94, 24.30, 40.08, 22.98, 39.09] # g/mol
    oxide_number_cations = [1, 1, 2, 1, 2, 1, 1, 1, 2, 2]
    oxide_number_oxygen = [2, 2, 3, 1, 3, 1, 1, 1, 1, 1]


    fe_mol = total_feo_value / (oxide_molar_weight[4] * oxide_number_cations[4])
    # fe_mol = n_fe2 + n_fe3 = (1-0.15)*fe_mol + 0.15 * fe_mol

    n_fe3 = fe_mol * x_fe3
    n_fe2 = fe_mol * (1-x_fe3)

    wt_FeO = n_fe2 / oxide_number_cations[4] * oxide_molar_weight[4]
    wt_Fe2O3 = n_fe3 / oxide_number_cations[5] * oxide_molar_weight[5]

    # print("wt Fe2 ", wt_FeO, "\n", "wt Fe3 ", wt_Fe2O3, "\n")

    return wt_FeO, wt_Fe2O3


class Therm_dyn_ther_looper:
    """
    A class representing a thermodynamic looping station for the Tunorrad object.

    Attributes:
        theriak_path (str): The path to the Theriak software.
        database (str): The name of the thermodynamic database.
        bulk_rock (str): The composition of the bulk rock.
        temperature (float): The temperature in Kelvin.
        pressure (float): The pressure in bars.
        df_var_dictionary (DataFrame): The variable dictionary.
        df_hydrous_data_dic (DataFrame): The hydrous data dictionary.
        df_all_elements (DataFrame): The dataframe containing all elements.
        num (int): The looping step.
        fluid_name_tag (str): The tag for the fluid name based on the database used.

    Methods:
        __init__(self, theriak_path, database, bulk_rock, temperature, pressure, df_var_dictionary, df_hydrous_data_dic, df_all_elements, num, fluid_name_tag):
            Initializes the Therm_dyn_ther_looper object.

        thermodynamic_looping_station(self, marco=False):
            Calls theriak and passes T, P and bulk rock. Resulting data is read and formatted.
            First check for 'water.fluid' is included.
            Prints messages and gives Volume, Density, Mol and Weight of present fluid as feedback.

        merge_dataframe_dic(self):
            Merges the calculated values from theriak output with the existing dataframes.

        step_on_water(self):
            Collecting data about solids and fluid from t=0 (new calculation) and t=-1 (previous calculation).
    """

    def __init__(
            self, theriak_path, database, bulk_rock,
            temperature, pressure, df_var_dictionary,
            df_hydrous_data_dic, df_all_elements, num, fluid_name_tag):
        """
        Initialize the Tunorrad object.

        Args:
            theriak_path (str): The path to the Theriak software.
            database (str): The name of the thermodynamic database.
            bulk_rock (str): The composition of the bulk rock.
            temperature (float): The temperature in Kelvin.
            pressure (float): The pressure in bars.
            df_var_dictionary (DataFrame): The variable dictionary.
            df_hydrous_data_dic (DataFrame): The hydrous data dictionary.
            df_all_elements (DataFrame): The dataframe containing all elements.
            num (int): The looping step.
            fluid_name_tag (str): The tag for the fluid name based on the database used.
        """
        self.theriak_path = theriak_path
        self.database = database
        # store variables from thermodynamic read
        self.df_phase_data = 0
        self.df_hydrous_data = 0
        self.phase_cache = {}
        self.hydrous_cache = {}
        self.frac_system = {}
        self.sol_sol_base = 0

        self.bulk_rock = bulk_rock

        self.df_all_elements = df_all_elements
        self.df_var_dictionary = df_var_dictionary
        self.df_hydrous_data_dic = df_hydrous_data_dic
        self.new_fluid_Vol = 0
        self.new_fluid_dens = 0
        self.new_fluid_N = 0
        self.new_fluid_weight = 0
        self.temperature = temperature
        self.pressure = pressure

        self.sys_vol_no_extraction = 0
        self.free_water_before = 0
        self.solid_vol_before = 0
        self.new_water = 0
        self.new_solid_Vol = 0

        self.g_sys = 0
        self.pot_frame = 0

        # garnet fraction separate
        self.separate = 0

        self.num = num  # looping step

        # passing the fluid name tag based on the database used
        self.fluid_name_tag = fluid_name_tag

    def thermodynamic_looping_station(
            self, marco=False, oversaturation=False,
                        theriak_input_rock_before=False
                        ):
        """
        Calls theriak and passes T, P and bulk rock. Resulting data is read and formatted.
        First check for 'water.fluid' is included.
        Prints messages and gives Volume, Density, Mol and Weight of present fluid as feedback.

        Parameters:
            marco (bool): Flag indicating whether to perform recalculation for oversaturation of water. Default is False.

        Returns:
            None
        """

        # Flag: is set to True if water.fluid is stable in the system for further calculations
        content_h2o = True

        # Calls function to read the theriak output
        # (function also runs theriak and passes temperature, pressure and whole rock)
        # - returns theriak_data (which is a dictionary with Dataframes)

        theriak_data, g_sys, pot_frame = read_theriak(
            self.theriak_path, self.database, self.temperature, self.pressure, self.bulk_rock,
            theriak_input_rock_before=theriak_input_rock_before)

        # recalculation to delete for oversaturation of water
        if marco is False and oversaturation is False:
            if self.num < 1:
                if self.fluid_name_tag in theriak_data['df_elements_in_phases'].columns:
                    # print("333333 - First step free fluid detected - recalc initialized")
                    new_bulk = {}
                    water_H = theriak_data['df_elements_in_phases'][self.fluid_name_tag]['H']
                    total = theriak_data['df_elements_in_phases']['total:']
                    # NOTE - rounding water hydrogen here - undersaturation now?
                    total['H'] = total['H'] - round(water_H, 5)
                    for el in total.index:
                        if el == 'O' or el == 'E':
                            pass
                        else:
                            new_bulk[el] = total[el]
                    reset_bulk = check_redo_bulk(new_bulk)
                    self.bulk_rock = reset_bulk
                    theriak_data, g_sys, pot_frame = read_theriak(
                        self.theriak_path, self.database, self.temperature, self.pressure, self.bulk_rock, 
                        theriak_input_rock_before=theriak_input_rock_before)

        self.g_sys = g_sys
        self.pot_frame = pot_frame

        # Stores returned data from "read_theriak" to variables -
        self.df_phase_data = theriak_data['df_Vol_Dens']
        self.df_all_elements = theriak_data['df_elements_in_phases']
        self.sol_sol_base = theriak_data['solid_solution_dic']

        # if self.temperature > 507.43:
        # print(self.temperature)
        # print("=3=== Test water presence in system...")
        if 'df_h2o_content' in theriak_data.keys():
            self.df_hydrous_data = theriak_data['df_h2o_content']
            if self.fluid_name_tag in self.df_hydrous_data.columns:
                content_h2o = True
                # print("=3=== free water is present.")
            else:
                content_h2o = False
                # print("3=== only hydrous solid phases, no free water.")
        else:
            content_h2o = False
            # Option for - No free water present in calculation
            self.df_hydrous_data = pd.DataFrame(
                {'empty': [0, 0, 0, 0, 0, 0, 0]})
            print("Status Quo = ########### No free H2O to extract #############")
        # print(self.df_hydrous_data)

        # Redefining volume and weight percentages if fluid phase is present
        # theriak only calculates Vol% for all soldi phases, no fluid phase involved

        if content_h2o is True:
            for phase in self.df_phase_data.columns:
                self.df_phase_data.loc['vol%', phase] = self.df_phase_data.loc[
                    'volume[ccm]', phase] / sum(
                    self.df_phase_data.loc['volume[ccm]', :]) * 100
                self.df_phase_data.loc['wt%', phase] = self.df_phase_data.loc[
                    'wt[g]', phase] / sum(
                    self.df_phase_data.loc['wt[g]', :]) * 100

        # Volume and Density ouput - Dataframes (df_N, df_Vol% etc)
        for variable in list(self.df_phase_data.index):
            self.phase_cache['df_' + str(variable)] = pd.DataFrame()

        # Dataframes for physical variables of hydrous phases
        if 'df_h2o_content' in theriak_data.keys():
            water_cont_ind = ["N", "H2O[pfu]", "H2O[mol]",
                              "H2O[g]", "wt%_phase", "wt%_solids", "wt%_H2O.solid"]
            for variable in water_cont_ind:
                self.hydrous_cache['df_' + str(variable)] = pd.DataFrame()

        # If fluid component is present loop is activated. Stores pyhsical properties of the fluid and the system.

        if content_h2o is True:
            # Created free water volume in system at P-T condition:
            self.new_fluid_Vol = self.df_phase_data.loc['volume[ccm]', self.fluid_name_tag]
            self.new_fluid_dens = self.df_phase_data.loc['density[g/ccm]', self.fluid_name_tag]
            self.new_fluid_N = self.df_phase_data.loc['N', self.fluid_name_tag]
            self.new_fluid_weight = self.df_phase_data.loc['wt[g]',
                                                           self.fluid_name_tag]

        # print("\t Information on the fluid:")
        # print('\t Vol:{} Dens:{} N:{} Weight:{}'.format(self.new_fluid_Vol,
        #       self.new_fluid_dens, self.new_fluid_N, self.new_fluid_weight))

    def merge_dataframe_dic(self):
        """
        Merges the calculated values from theriak output with the existing dataframes.

        This method merges the newly calculated phase data from the phase cache with the
        existing dataframe for each physical property ('N', 'volume', 'density', etc).
        It also merges the hydrous solids and fluid phases from the hydrous cache with
        the existing hydrous data dataframe. Finally, it merges the element data with
        the existing dataframe of all elements.

        Returns:
            None
        """
        # Merging calculated values from theriak output with concat into (first empty) dataframe
        # This is added each P-T-step the new calculated values to the
        # dataframe of the physical property and the stable phases

        for num, key in enumerate(self.phase_cache.keys()):
            self.df_var_dictionary[key] = pd.concat(
                [self.df_var_dictionary[key], self.df_phase_data.iloc[num, :]], axis=1)

        # Merging hydrous solids and fluid phases
        for num, key in enumerate(self.hydrous_cache.keys()):
            if len(self.df_hydrous_data.index) == len(self.hydrous_cache.keys()):
                self.df_hydrous_data_dic[key] = pd.concat(
                    [self.df_hydrous_data_dic[key], self.df_hydrous_data.iloc[num, :]], axis=1)
            else:
                if key == 'df_N':
                    print("\t Merging dataframe information:")
                    # print("\t Phases not saved:{}".format(self.df_hydrous_data.columns))
                    # print("\t Array of zeros will be added to existing data")
                # creating a zero DataFrame similar to last condition to fill gap
                d_zero = pd.DataFrame(
                    np.zeros(len(self.df_hydrous_data_dic[key].index)))
                d_zero.index = self.df_hydrous_data_dic[key].index
                self.df_hydrous_data_dic[key] = pd.concat(
                    [self.df_hydrous_data_dic[key], d_zero], axis=1)

        # merging element data
        self.df_all_elements = pd.concat(
            [self.df_all_elements, self.df_all_elements.loc[:, 'total:']], axis=1)

    def step_on_water(self):
        """
        Collecting data about solids and fluid from t=0 (new calculation) and t=-1 (previous calculation)

        This method calculates and collects data related to solids and fluid at two time points: t=0 (new calculation) and t=-1 (previous calculation).
        It performs the following steps:
        1. Calculates the total volume of the system without extraction at t=0.
        2. Retrieves the volume data before the current time step from the variable dictionary.
        3. Compares the water data between the current time step and the previous time step.
        4. Determines the volume of free water before the current time step.
        5. Calculates the volume of solids before the current time step.
        6. Calculates the new volume of solids by subtracting the volume of free water from the total volume of the system at t=0.
        """
        self.sys_vol_no_extraction = self.df_phase_data.loc['volume[ccm]'].sum(
        )
        try:
            volume_data_before = np.array(
                self.df_var_dictionary['df_volume[ccm]'].iloc[:, -2])
            volume_data_before = np.nan_to_num(volume_data_before)
        except IndexError:
            volume_data_before = np.array(
                self.df_var_dictionary['df_volume[ccm]'].iloc[:, -1])
            volume_data_before = np.nan_to_num(volume_data_before)

        # comparing water data (momentary step and previous step)
        # important for extraction steps (factor, volume changes, porosity)
        if self.fluid_name_tag in list(self.df_all_elements.columns):
            fluid_volumes = np.array(
                self.df_var_dictionary['df_volume[ccm]'].loc[self.fluid_name_tag])
            try:
                # it is essential to state the "== False" otherwise the value is not read
                if np.isnan(fluid_volumes[-2]) == False:
                    # print("cache water[-2] is a number")
                    self.free_water_before = fluid_volumes[-2]
                else:
                    # print("cache water[-2] is NaN")
                    self.free_water_before = 0
            # no previous step possible - first oocurence of free water
            except IndexError:
                # print("first round fluid vol")
                self.free_water_before = fluid_volumes[-1]

            if np.isnan(self.free_water_before) is True:
                self.free_water_before = 0
        else:
            self.free_water_before = 0

        self.solid_vol_before = volume_data_before.sum() - self.free_water_before
        self.new_solid_Vol = sum(
            list(self.df_phase_data.loc['volume[ccm]'])) - self.new_fluid_Vol

    def system_condi(self):
        """
        Collecting the system data from calculations

        This method collects various system data from calculations, such as volume, density, weight, number of moles,
        volume percentage, and weight percentage. It then performs some data cleaning by filling any missing values with 0.

        Finally, it calculates the total volume of stable phases for each temperature (and pressure) step, as well as the
        solid phase volume (excluding the fluid phase if present) for each temperature (and pressure) step.
        """

        sys_vol_data = self.df_phase_data.loc['volume[ccm]']
        sys_dens_data = self.df_phase_data.loc['density[g/ccm]']
        sys_weight_data = self.df_var_dictionary['wt[g]']
        sys_mol_data = self.df_phase_data.loc['N']
        sys_volperc_data = self.df_phase_data.loc['vol%']
        sys_weightperc_data = self.df_phase_data.loc['wt%']

        sys_vol_data = sys_vol_data.fillna(0, inplace=True)
        sys_dens_data = sys_dens_data.fillna(0, inplace=True)
        sys_weight_data = sys_weight_data.fillna(0, inplace=True)
        sys_mol_data = sys_mol_data.fillna(0, inplace=True)
        sys_volperc_data = sys_volperc_data.fillna(0, inplace=True)
        sys_weightperc_data = sys_weightperc_data.fillna(0, inplace=True)

        # system volume conditions
        for temperature in list(self.sys_vol_data.index):
            # Sums up total volume of stable phases for each temperature (and pressure) step
            tot_volume = self.sys_vol_data.loc[temperature, :].sum()
            self.sys_vol_pre.append(tot_volume)

            if self.fluid_name_tag in list(self.sys_vol_data.columns):
                # Sums up solid phase volume for each temperature (and pressure) step
                solid_volume = self.sys_vol_data.loc[temperature, :].sum(
                ) - self.sys_vol_data.loc[temperature, self.fluid_name_tag]
                self.solid_vol_recalc.append(solid_volume)
            else:
                solid_volume = self.sys_vol_data.loc[temperature, :].sum()
                self.solid_vol_recalc.append(solid_volume)

    def mineral_fractionation(self, oxygen_data, name):
        """
        Calculates the fractionation of a mineral from the bulk composition.

        Parameters:
        oxygen_data (dict): Dictionary containing oxygen data.
        name (str): Name of the mineral to be fractionated.

        Returns:
        float: The new oxygen bulk composition after fractionation.
        """

        # diminish double total: in dataframe
        if self.df_all_elements.columns[-1] == self.df_all_elements.columns[-2]:
            self.df_all_elements = self.df_all_elements.iloc[:, :-1]
        phase_list = list(self.df_all_elements.columns)
        corr_list = []
        # test stable phases for mineral to be fractionated from bulk
        for phase in phase_list:
            if '_' in phase:
                pos = phase.index('_')
                phase = phase[:pos]
            corr_list.append(phase)

        # fractionation in the case the phase is present
        if name in corr_list:
            # postion to read the to be fractionated phase from DataFrame
            min_pos = corr_list.index(name)

            # store phase composition before fractionation
            el_raw_data = [np.array(self.df_all_elements[phase_list[min_pos]]),
                           self.df_all_elements[phase_list[min_pos]].index]
            # store to dataclass
            self.separate = Phasecomp(phase_list[min_pos], self.temperature, self.pressure,
                                      self.df_phase_data[phase_list[min_pos]]['N'],
                                      self.df_phase_data[phase_list[min_pos]]['volume[ccm]'],
                                      self.df_phase_data[phase_list[min_pos]]['vol%'],
                                      self.df_phase_data[phase_list[min_pos]]['wt[g]'],
                                      self.df_phase_data[phase_list[min_pos]]['wt%'],
                                      self.df_phase_data[phase_list[min_pos]]['density[g/ccm]'],
                                      el_raw_data,
                                      self.df_phase_data[phase_list[min_pos]]['volume/mol']
                                      )

            # element data of phase to be fractionated
            frac_mineral_el = self.df_all_elements[phase_list[min_pos]]
            # total element bulk
            active_bulk = self.df_all_elements['total:']
            self.frac_system[self.temperature] = active_bulk
            # elements substracted
            minus_element = active_bulk - frac_mineral_el
            # update the system bulk used for next bulk calculation
            self.df_all_elements['total:'] = minus_element
            # print(f"Fractionated phase:\n{phase_list[min_pos]}")
            # update the bulk rock oxygen
            collect_phases = []
            for phase in oxygen_data['Phases']:
                collect_phases.append(self.df_all_elements.loc['O'][phase])

            phase_oxygen = np.array(collect_phases)
            phase_oxygen[min_pos] = 0
            # ANCHOR syntax fix pandas - caviat from jupyter notebook
            self.df_all_elements.loc['O', phase_list[min_pos]] = 0
            # self.df_all_elements.loc['O'][min_pos] = 0
            # self.df_all_elements.loc['O'].iloc[min_pos] = 0
            phase_doxy = oxygen_data['delta_O']

            new_O_bulk = sum(phase_oxygen*phase_doxy / sum(phase_oxygen))

            return new_O_bulk



class ThermodynamicPressureSolver:
    """
    Solves the coupled system of equations (5), (6), and (9) to find
    thermodynamically consistent pressure evolution.
    """
    
    def __init__(self, thermodynamic_data_t0, thermodynamic_data_t1, dt=1.0):
        """
        Initialize with thermodynamic data from two Gibbs minimizations.
        
        Parameters:
        -----------
        thermodynamic_data_t0 : dict
            Data at t=0 containing: rho_f, rho_s, X_h, T, P
        thermodynamic_data_t1 : dict  
            Data at t=1 containing: rho_f, rho_s, X_h, T, P
        dt : float
            Time step
        """
        self.data_t0 = thermodynamic_data_t0
        self.data_t1 = thermodynamic_data_t1
        self.dt = dt
        
        # Extract values at t=0
        self.rho_f_t0 = thermodynamic_data_t0['rho_f']
        self.rho_s_t0 = thermodynamic_data_t0['rho_s']
        self.X_h_t0 = thermodynamic_data_t0['X_h']
        self.T_t0 = thermodynamic_data_t0['T']
        self.P_t0 = thermodynamic_data_t0['P']
        
        # Extract values at t+dt (these are from thermodynamic equilibrium)
        self.rho_f_t1 = thermodynamic_data_t1['rho_f']
        self.rho_s_t1 = thermodynamic_data_t1['rho_s']
        self.X_h_t1 = thermodynamic_data_t1['X_h']
        self.T_t1 = thermodynamic_data_t1['T']
        self.P_t1_thermo = thermodynamic_data_t1['P']  # This is from your thermo solver
        
        print("=== THERMODYNAMIC PRESSURE SOLVER ===")
        print(f"Initial state (t=0): P={self.P_t0:.3f}, T={self.T_t0:.1f}, ρf={self.rho_f_t0:.3f}, ρs={self.rho_s_t0:.3f}, Xh={self.X_h_t0:.6f}")
        print(f"Final state (t+dt): P={self.P_t1_thermo:.3f}, T={self.T_t1:.1f}, ρf={self.rho_f_t1:.3f}, ρs={self.rho_s_t1:.3f}, Xh={self.X_h_t1:.6f}")
    
    def calculate_porosity(self, rho_s, X_h, rho_f):
        """Calculate porosity from equation (3): φf = (ρs×Xh)/[(1-Xh)×ρf]"""
        return (rho_s * X_h) / ((1 - X_h) * rho_f)

    def solve_simplified_approach(self):
        """
        Simplified approach using direct thermodynamic data.
        
        Since you have thermodynamic equilibrium data at both time points,
        you can use the thermodynamic consistency condition.
        """
        print(f"\n=== SIMPLIFIED APPROACH ===")
        
        # Method 1: Direct thermodynamic result
        P_direct = self.P_t1_thermo
        print(f"Method 1 (Direct thermo): P = {P_direct:.6f}")
        
        # Method 2: Pressure term evolution
        pressure_term_t0 = (self.rho_s_t0 * self.X_h_t0) / (1 - self.X_h_t0)
        pressure_term_t1 = (self.rho_s_t1 * self.X_h_t1) / (1 - self.X_h_t1)
        pressure_term_ratio = pressure_term_t1 / pressure_term_t0
        
        P_evolution = self.P_t0 * pressure_term_ratio
        print(f"Method 2 (Pressure term): P = {P_evolution:.6f} (ratio = {pressure_term_ratio:.6f})")
        
        # Method 3: Density-weighted average
        density_ratio = self.rho_f_t1 / self.rho_f_t0
        weight_ratio = self.X_h_t1 / self.X_h_t0
        
        combined_ratio = (pressure_term_ratio * density_ratio * weight_ratio)**(1/3)
        P_combined = self.P_t0 * combined_ratio
        print(f"Method 3 (Combined): P = {P_combined:.6f} (ratio = {combined_ratio:.6f})")
        
        return P_direct, P_evolution, P_combined
    
    def validate_solution(self, P_solution):
        """Validate the solution by checking thermodynamic consistency."""
        print(f"\n=== SOLUTION VALIDATION ===")
        
        # Calculate porosity at both times
        phi_t0 = self.calculate_porosity(self.rho_s_t0, self.X_h_t0, self.rho_f_t0)
        phi_t1 = self.calculate_porosity(self.rho_s_t1, self.X_h_t1, self.rho_f_t1)
        
        print(f"Porosity evolution: {phi_t0:.6f} → {phi_t1:.6f}")
        
        # Check mass conservation
        mass_term_t0 = self.rho_s_t0 * (1 - self.X_h_t0) * (1 - phi_t0)
        mass_term_t1 = self.rho_s_t1 * (1 - self.X_h_t1) * (1 - phi_t1)
        mass_error = abs(mass_term_t1 - mass_term_t0)
        
        print(f"Mass conservation: {mass_term_t0:.6f} → {mass_term_t1:.6f}")
        print(f"Mass conservation error: {mass_error:.6f}")
        
        # Check pressure evolution
        pressure_change = P_solution - self.P_t0
        relative_change = (P_solution / self.P_t0 - 1) * 100
        
        print(f"Pressure evolution: {self.P_t0:.6f} → {P_solution:.6f}")
        print(f"Pressure change: {pressure_change:.6f} ({relative_change:.2f}%)")
        
        return {
            'porosity_t0': phi_t0,
            'porosity_t1': phi_t1,
            'mass_conservation_error': mass_error,
            'pressure_change': pressure_change,
            'relative_change': relative_change
        }



class Ext_method_master:
    """
    Module to calculate the factor important for the extraction of the water from the system

    Attributes:
        rock_item (int): The tag for the rock item.
        fluid_t1 (int): The fluid volume from the new P-T step.
        fluid_t0 (int): The fluid volume from the previous P-T step.
        solid_t1 (int): The solids volume from the new P-T step.
        solid_t0 (int): The solids volume from the previous P-T step.
        master_norm (list): The normalization of the values due to the slope of P-T steps.
        save_factor (list): The factor regulating the extraction - saved to a list.
        phase_data (dict): The phase data.
        unlock_freewater (bool): Flag indicating if the freewater is unlocked.
        pressure (float): The pressure.
        tensile_strength (float): The tensile strength.
        fracture (bool): Flag indicating if there is fracture.
        shear_stress (float): The shear stress.
        frac_respo (int): The response to fracturing.
        angle (float): The subduction angle.
        diff_stress (float): The differential stress.
        friction (float): The friction.
        failure_dictionary (dict): The failure dictionary.
        fluid_name_tag (str): The fluid name tag.
        fluid_pressure_mode (str): The fluid pressure mode.

    Methods:
        __init__: Initialize all the values and data necessary for calculations.
        couloumb_method: Extraction of water/fluid follows the overstepping of a certain factor.
        couloumb_method2: Calculate the Coulomb test method 2 for rock failure.
    """

    def __init__(
            self, pressure, temperature,
            moles_vol_fluid,
            fluid_volume_before, fluid_volume_new,
            solid_volume_before, solid_volume_new,
            save_factor, master_norm, phase_data,
            tensile_s, differential_stress, friction, subduction_angle,
            reviewer_mode, phase_data_complete, hydrous_data_complete,
            pressure_before,
            fluid_pressure_mode, fluid_name_tag, extraction_connectivity=0.0, extraction_threshold=False, rock_item_tag=0
            ):
        """
        Initialize all the values and data necessary for calculations

        Args:
            fluid_volume_before (int): value of the fluid volume from previous P-T step
            fluid_volume_new (int): value of the fluid volume from new P-T step
            solid_volume_before (int): value of the solids volume from previous P-T step
            solid_volume_new (int): value of the solids volume from new P-T step
            save_factor (list): factor regulating the extraction - saved to a list
            master_norm ([type]): normalization  of the values due to the slope of P-T steps
        """
        # REVIEW fluid t0 and t1 where twisted?
        self.rock_item = rock_item_tag
        self.fluid_t1 = fluid_volume_new
        self.fluid_t0 = fluid_volume_before
        self.solid_t1 = solid_volume_new
        self.solid_t0 = solid_volume_before
        self.master_norm = master_norm
        self.save_factor = save_factor
        self.phase_data = phase_data
        self.unlock_freewater = False
        self.pressure = pressure
        self.temperature = temperature
        self.pressure_before = pressure_before
        self.moles_vol_fluid = moles_vol_fluid
        # assumption and connected to main code - should be sent by input
        self.tensile_strength = tensile_s
        self.fracture = False
        self.shear_stress = 0
        self.frac_respo = 0
        self.angle = subduction_angle
        self.diff_stress = differential_stress
        self.friction = friction
        self.failure_dictionary = {}
        self.fluid_name_tag = fluid_name_tag
        self.fluid_pressure_mode = fluid_pressure_mode
        self.extraction_threshold = extraction_threshold
        self.extraction_connectivity = extraction_connectivity
        self.reviewer_mode = reviewer_mode
        self.phase_data_complete = phase_data_complete
        self.hydrous_data_complete = hydrous_data_complete
        self.pressure_results = []

    def couloumb_method(self, t_ref_solid, tensile=20):
        """
        Extraction of water/fluid follows the overstepping of a certain factor - calculated similar to Etheridge (2020)

        Args:
            t_ref_solid (float): Reference temperature for solid.
            tensile (float, optional): Tensile strength. Defaults to 20.

        Returns:
            None
        """

        print('\tsolid_t0:{} solid_t-1:{} fluid_t0:{} fluid_t-1:{}'.format(self.solid_t0,
              self.solid_t1, self.fluid_t0, self.fluid_t1))
        print('\t-----------------')

        #################################
        # get system conditions at present step and previous condition
        vol_t0 = self.solid_t0 + self.fluid_t0
        vol_new = self.solid_t1 + self.fluid_t1
        # define lithostatic pressure and sys pressure - convert Bar to MPa = one magnitude smaller
        litho = self.pressure/10
        rock = litho * vol_new/vol_t0
        print('\tP_Lith:{} P_rocksys:{}'.format(litho, rock))
        # test whether is sigma 1 or sigma 3, it defines the stress regime
        if rock > litho:
            sig1 = rock
            sig3 = litho
        elif litho > rock:
            # the tensile overpressure here
            # hydraulic extension fracture when diff. stress < 4T
            # shear failure when diff. stress > 4T and radius*2theta of mohr circle equal or larger than envelope (TETAkrit)
            sig1 = litho
            sig3 = rock
        else:
            # probably not happening - to be tested - no fracturing?
            sig1 = litho
            sig3 = rock

        # Mohr Circle radius and centre
        diff_stress = sig1-sig3
        self.diff_stress = diff_stress
        r = abs(sig1-sig3)/2
        center = sig1 - r
        pos = center - rock

        # linear quation of the critical line
        # 50 for y axis intercept: large depth with 200MPa < normal stress < 2000 MPa
        # slope of 0.6
        cohesion = 50
        # FIXME - internal friction coefficient modified
        # internal_friction = 0.6
        internal_friction = 0.06
        a = -1
        b = 1/internal_friction
        c = -cohesion/internal_friction

        # criterion for failure via extensional fracturing or shear failure - what is hitted earlier?
        output = checkCollision(a, b, c, x=pos, y=0, radius=r)
        self.frac_respo = 0
        if diff_stress > 0:
            print("\tHydrofracture test...")

            # Hydrofracture criterion - double checked if envelope was hit earlier than tensile strength
            if (sig3-rock) <= -self.tensile_strength:
                tpos = -self.tensile_strength+abs(pos)
                # Test if envelope collision was ealier
                outputb = checkCollision(a, b, c, x=tpos, y=0, radius=r)
                print(f"\t-->{outputb}")
                if outputb == 'Intersect':
                    fracturing = True
                    self.frac_respo = 2
                    print(f"\t-->Shear-failure")
                elif outputb == 'Touch':
                    fracturing = True
                    self.frac_respo = 3
                    print(f"\t-->Ultimative failure - extension plus shear!?!?")
                # if all is False, there is extensional failure
                else:
                    fracturing = True
                    self.frac_respo = 1
                    print(f"\t-->Extensional fracturing")

            # The shear failure because the circle did not cross the tensile strength limit
            elif output == 'Touch' or output == 'Intersect':
                fracturing = True
                self.frac_respo = 2
                print(f"\t-->Shear-failure")

            # All failure criterions are False - there is no fracturing
            else:
                print("\t...nothing happened...")
                fracturing = False
                self.frac_respo = 0

        # No differential stress - there is no fracturing
        else:
            print("\t...no diff. stress, nothing happened...")
            fracturing = False
            self.frac_respo = 0

        # Update system condition of fracturing
        self.fracture = fracturing
        print("End fracture modul")

    def couloumb_method2(self, shear_stress, friction, cohesion):
        """
        Calculate the Coulomb test method 2 for rock failure.

        Args:
            shear_stress (float): The shear stress applied to the rock.
            friction (float): The internal friction angle of the rock.
            cohesion (float): The cohesion of the rock.

        Returns:
            None

        Raises:
            None

        Notes:
            - This method calculates the critical fluid pressure and checks for brittle shear failure.
            - It also performs a double check for Mohr-Coulomb failure, considering both brittle extension and brittle shear.
            - The results are printed to the console.
        """

        print("Coulomb test method 2 answer:")
        print('\tsolid_t0:{}\n\tsolid_t1:{}\n\tfluid_t0:{}\n\tfluid_t1:{}'.format(self.solid_t0,
              self.solid_t1, self.fluid_t0, self.fluid_t1))

        # test before executing module - phase data fluid volume should be equal the self.fluid_t1
        if self.phase_data[self.fluid_name_tag]['volume[ccm]'] != self.fluid_t1:
            print("Inconsistency in new fluid volume")
            # keyboard.wait('esc')

        #################################
        # Mohr-Coulomb rock yield
        # linear equation of the critical line
        # 50 for y axis intercept: large depth with 200MPa < normal stress < 2000 MPa
        # slope of 0.6
        # LINK Mohr-Coulomb slope for failure
        cohesion = 50
        internal_friction = 0.7

        cohesion = cohesion
        internal_friction = friction
        a = -1
        b = 1/internal_friction
        c = -cohesion/internal_friction

        # REVIEW - static fix to 45°
        self.angle = 45

        # #########################################
        # New method 16.02.2023 - brittle shear failure - Cox et al. 2010
        # define: lithostatic pressure, differential stress, sigma1, sigma3, normal stress
        litho = self.pressure/10 # convert Bar to MPa

        # differential stress from shear stress input, recasted after Cox et al. 2010
        self.diff_stress = 2*shear_stress/np.sin(2*self.angle*np.pi/180)

        # NOTE - setting lithostatic pressure for sig1
        # sigma 1 or sigma 3, it defines the stress regime
        sig1 = litho
        sig3 = litho - self.diff_stress

        # normal stress after Cox et al. 2010
        # normal stress = lithostatic when theta is 45°
        normal_stress = ((sig1+sig3)/2) - ((sig1-sig3)/2) * np.cos(2*self.angle*np.pi/180)

        # REVIEW sig1-sig3 version 18.06.2023 before revision from workshop
        """# sigma 1 or sigma 3, it defines the stress regime
        sig1 = litho + diff_stress/2
        sig3 = litho - diff_stress/2

        # normal stress after Cox et al. 2010
        # normal stress = lithostatic when theta is 45°
        normal_stress = ((sig1+sig3)/2) - ((sig1-sig3)/2) * np.cos(2*self.angle*np.pi/180)

        # reevaluate sig1 and sig3 after assigning normal_stress
        # remain same values if theta is 45°
        sig1 = normal_stress + diff_stress/2
        sig3 = normal_stress - diff_stress/2"""

        # Critical fluid pressure
        crit_fluid_pressure = cohesion/internal_friction + ((sig1+sig3)/2) - (
                (sig1-sig3)/2)*np.cos(2*self.angle*np.pi/180) - ((sig1-sig3)/2/internal_friction)*np.sin(2*self.angle*np.pi/180)

        # #########################################
        # Fluid pressure
        # get system conditions at present step and previous step
        vol_t0 = self.solid_t0 + self.fluid_t0
        vol_new = self.solid_t1 + self.fluid_t1

        # Fluid pressure calculation
        # Duesterhoft 2019 method
        hydro2 = normal_stress + normal_stress/vol_t0 * (vol_t0-(vol_new-vol_t0))-normal_stress
        # purely volume ratio method
        hydro = normal_stress * vol_new/vol_t0
        # NOTE fluid pressure after Duesterhoft 2019
        hydro = hydro2

        # brittle shear failure check
        print("Critical fluid pressure:{:.3f} Calc.-fluid pressure:{:.3f} difference:{:.3f}".format(
                        crit_fluid_pressure, hydro, crit_fluid_pressure-hydro))
        print("\tBrittle shear failure test...")
        if hydro >= crit_fluid_pressure:
            fracturing = True
            self.frac_respo = 2
            print(f"\t-->Brittle-Shear-Failure")
        # critical fluid pressure is not reached
        else:
            print("\t...not reaching critical value...nothing happens...")
            fracturing = False
            self.frac_respo = 0

        # #########################################
        # Mohr-Coulomb failure double check - brittle extension or brittle shear

        # Mohr Circle radius and centre
        self.shear_stress = shear_stress # in MPa
        r = abs(sig1-sig3)/2
        center = normal_stress
        pos = center - hydro

        # macroscopic griffith criterion

        normal_stress_line = np.linspace(-cohesion/2, 0, 100)
        fail_tau = np.sqrt(4*(normal_stress_line)*(cohesion/2)+4*(cohesion/2)**2)

        pf_crit_griffith = (8*self.tensile_strength*(sig1+sig3)-((sig1-sig3)**2))/(16*self.tensile_strength)

        # Test possible extensional fracturing
        output = checkCollision(a, b, c, x=pos, y=0, radius=r)

        # Condition 4T >= differential stress for extensional brittle failure
        if self.tensile_strength*4 > self.diff_stress:
            print("\tRe-evaluate failure test circle...")

            # Hydrofracture criterion - double checked if envelope was hit earlier than tensile strength
            if (sig3-hydro) <= -self.tensile_strength:

                # Test if shear envelope collision was ealier
                tpos = -self.tensile_strength+abs(pos)
                outputb = checkCollision(a, b, c, x=tpos, y=0, radius=r)
                print(f"\t-->{outputb}")

                # Shear
                if outputb == 'Intersect':
                    fracturing = True
                    self.frac_respo = 2
                    print(f"\t-->Shear-failure")

                # Hybrid shear
                elif outputb == 'Touch':
                    fracturing = True
                    self.frac_respo = 3
                    print(f"\t-->Ultimative failure - extension plus shear!?!?")

                # Brittle extensional
                else:
                    fracturing = True
                    self.frac_respo = 1
                    print(f"\t-->Extensional fracturing")

        # NOTE testing with plot
        # arrays for plotting the Mohr-Coloumb diagramm
        normal_stress_line = np.linspace(-60, 2000, 100)
        tkrit = cohesion + internal_friction*normal_stress_line
        tau_griffith = np.sqrt(4*(cohesion/2)*(normal_stress_line+(cohesion/2)))
        # mohr circle
        theta = np.linspace(0, 2*np.pi, 100)
        x1f = r*np.cos(theta) + pos
        x2f = r*np.sin(theta)
        x1f2 = r*np.cos(theta) + center
        plt.figure(1003)
        plt.plot(normal_stress_line, tkrit, 'r-', x1f, x2f, 'b--', x1f2, x2f, 'g-')
        plt.plot(normal_stress_line, tau_griffith, 'r--')
        label = "{:.2f}".format(self.diff_stress)
        plt.annotate(label, (sig3-sig1, 0),
                    textcoords="offset points", xytext=(0, 10), ha='center')
        plt.axvline(color='red', x=-cohesion/2)
        plt.axvline(color='black', x=0)
        plt.axhline(color='black', y=0)
        plt.xlabel(r"$\sigma\ MPa$")
        plt.ylabel(r"$\tau\ MPa$")
        plt.ylim(-1, 100)
        plt.xlim(-125, 200)
        plt.show()
        """
        if self.rock_item == 'rock1':
            # arrays for plotting the Mohr-Coloumb diagramm
            normal_stress_line = np.linspace(-60, 2000, 100)
            tkrit = cohesion + internal_friction*normal_stress_line
            # mohr circle
            theta = np.linspace(0, 2*np.pi, 100)
            x1f = r*np.cos(theta) + pos
            x2f = r*np.sin(theta)
            x1f2 = r*np.cos(theta) + center
            plt.figure(1003)
            plt.plot(normal_stress_line, tkrit, 'r-', x1f, x2f, 'b--', x1f2, x2f, 'g-')
            label = "{:.2f}".format(self.diff_stress)
            plt.annotate(label, (sig3-sig1, 0),
                        textcoords="offset points", xytext=(0, 10), ha='center')
            plt.axvline(color='red', x=-self.tensile_strength)
            plt.axvline(color='black', x=0)
            plt.axhline(color='black', y=0)
            plt.xlabel(r"$\sigma\ MPa$")
            plt.ylabel(r"$\tau\ MPa$")
            plt.ylim(-1, 100)
            plt.xlim(-125, 200)
            # plt.show()
        """

        # Update system condition of fracturing
        self.fracture = fracturing
        print("End fracture modul")

    def mohr_cloulomb_griffith(self, shear_stress=False):
        """
        Estimates the differential stress, tensile strength, and fluid pressure
        based on the Mohr-Coulomb failure criterion and Griffith's criterion.

        Args:
            shear_stress (bool): Flag indicating whether to use shear stress as input.
                If False, the differential stress is directly input.

        Returns:
            None
        """

        # 26.06.2023
        # Estimates on
        # 1) differential stress (sig1 is lithostatic, shear resolves sig3),
        # 2) tensile strength (estimates from literature and experiments therein),
        # 3) fluid pressure from volume change (Duesterhoft 2019)
        # Blanpied et al. 1992 state that low effective stresses can lead to slide failure
        # or even extensional failure eventhough the system has low differential stresses.
        # Extraction of water/fluid when the shear envelope or tensile strength are intersectedf
        # - idea similar to Etheridge (2020) - no fluid factor approach


        # NOTE testing with plot - Not working after Cassis update (23.01.2024)
        # arrays for plotting the Mohr-Coloumb diagramm
        def mcg_plot(cohesion, internal_friction, diff_stress, shear_stress, sig1, hydro):

            # ##########################################
            # Mohr Circle
            r = diff_stress/2
            center = sig1-r
            pos = center - hydro

            # plotting conditions
            stress_line = np.linspace(-60, 2000, 10000)
            tkrit = cohesion + internal_friction*stress_line
            tau_griffith = np.sqrt(4*(cohesion/2)*(stress_line+(cohesion/2)))
            # mohr circle
            theta = np.linspace(0, 2*np.pi, 100)
            x1f = r*np.cos(theta) + pos
            x2f = r*np.sin(theta)
            x1f2 = r*np.cos(theta) + center
            plt.figure(10001)
            plt.plot(stress_line, tkrit, 'r-', x1f, x2f, 'b--', x1f2, x2f, 'g-')
            plt.plot(stress_line, tau_griffith, 'r--')
            label = "{:.2f}".format(self.diff_stress)
            plt.annotate(label, (sig3-sig1, 0),
                        textcoords="offset points", xytext=(0, 10), ha='center')
            plt.axvline(color='red', x=-cohesion/2)
            plt.axvline(color='black', x=0)
            plt.axhline(color='black', y=0)
            plt.xlabel(r"$\sigma\ MPa$")
            plt.ylabel(r"$\tau\ MPa$")
            plt.ylim(-1, 100)
            plt.xlim(-60, 200)
            plt.show()


        print("M-C-G test module answer:")
        # print('\tsolid_t0:{}\n\tsolid_t1:{}\n\tfluid_t0:{}\n\tfluid_t1:{}'.format(self.solid_t0,
        #       self.solid_t1, self.fluid_t0, self.fluid_t1))

        # test before executing module - phase data fluid volume should be equal the self.fluid_t1
        if self.fluid_name_tag in self.phase_data.columns:
            if self.phase_data[self.fluid_name_tag]['volume[ccm]'] != self.fluid_t1:
                print("Inconsistency in new fluid volume")
                keyboard.wait('esc')
        else:
            print("No fluid phase data available")

        #################################
        # Mohr-Coulomb rock yield
        # linear equation of the critical line
        # 50 for y axis intercept: large depth with 200MPa < normal stress < 2000 MPa
        # slope of 0.6
        # LINK Mohr-Coulomb slope for failure
        cohesion = 2* self.tensile_strength

        # REVIEW - static fix to 27° for the optimal angle of failure after Cox 2010 and Sibson 2000
        # 45 degree gives that diff stress is two time shear stress (and not more)
        # friction of 0.75 gives ~26.6°
        self.angle = np.round(0.5 * np.arctan(1/self.friction) *180/np.pi,1)

        theta = self.angle*np.pi/180
        # #########################################
        # New method 16.02.2023 - brittle shear failure - Cox et al. 2010
        # define: lithostatic pressure, differential stress, sigma1, sigma3, normal stress
        litho = self.pressure/10 # convert Bar to MPa

        # NOTE - differential stress direct input from init file (23.01.2024)
        # input of shear stress is deactivated, direct input of differential stress (15.01.2024)
        # differential stress is taken as 2*shear stress - fix before proper diff stress input from input file
        # self.diff_stress = 2*shear_stress
        shear_stress = self.diff_stress*np.sin(2*theta)/2
        """
        # differential stress from shear stress input, recasted after Cox et al. 2010
        # 45 degree gives that diff stress is two time shear stress (and not more)
        self.diff_stress = 2*shear_stress/np.sin(2*theta)"""

        # NOTE - setting lithostatic pressure for sig3
        # sigma 1 or sigma 3, it defines the stress regime
        sig3 = litho
        sig1 = litho + self.diff_stress

        # #########################################
        # Normal stress of the system defined after Cox et al 2010
        # normal stress after Cox et al. 2010
        # normal stress = lithostatic when theta is 90°
        normal_stress = ((sig1+sig3)/2) - ((sig1-sig3)/2) * np.cos(2*theta)
        mean_stress = litho + (sig1-sig3)/2
        # #########################################
        # Fluid pressure
        # get system conditions at present step and previous step
        vol_t0 = self.solid_t0 + self.fluid_t0
        vol_new = self.solid_t1 + self.fluid_t1

        # ####################################
        # CORK-real - Calculate P2/P1

        R = 0.083144621  # cm³kbarK⁻¹mol⁻¹ (gas constant)
        T = 298+self.temperature  # K (temperature)
        #n = self.moles_fluid  # moles
        #V1 = self.fluid_t1  # cm^3 (initial volume)
        Vm1 = self.moles_vol_fluid

        # Cork parameters - values from Holland and Powell 1991 -
        # A Compensated-Redlich-Kwong (CORK) equation for volumes and fugacities of CO2 and H2O in the range 1 bar to 50 kbar and 100-1600°C
        a0 = 113.4
        a4 = -0.22291
        a5 = -3.8022 * 10**-4
        a6 = 1.7791 * 10**-4
        a = a0 + a4 * (T-673) + a5 * (T-673)**2 + a6 * (T-673)**3

        b = 1.465  # kJ kbar^-1 mol^-1, converted to L/mol

        c0 = -3.025650 * 10**-2
        c1 = -5.343144*10**-6
        c = c0 + c1 * T

        d0 = -3.2297554*10**-3
        d1 = 2.2215221*10**-6
        d = d0 + d1 * T

        """Vm1 = Vm1 + c * np.sqrt(P-P0) + d*(P-P0)"""

        if self.solid_t1 == vol_t0 and self.fluid_t1 == 0:
            V2 = self.fluid_t1
            p2_p1_cork_real = 0
        elif self.fluid_t1 == 0:
            V2 = vol_t0-self.solid_t1
            p2_p1_cork_real = 0
        else:
            comp_term = c * np.sqrt(litho-0.2) + d*(litho-0.2)
            # Assuming constant number of moles: n = V0 / Vm1 = V1 / Vm2
            # Solving for Vm2 (V2 here), based on current fluid volume
            # litho - 0.2 is the compensation term and 0.2 is pressure in bar as the value when MRK is deviating at high pressures
            V2 = 1/(self.fluid_t1/(vol_t0-self.solid_t1)) * (Vm1) # - comp_term
            
            # check the mass conservation
            moles_0 = self.fluid_t1 / Vm1
            moles_1 = (vol_t0-self.solid_t1) / V2
            # assert when mass is not conserved
            # NOTE - this is not working for the case when fluid_t1 = 0
            assert np.isclose(moles_0, moles_1, rtol=1e-4), "Mass not conserved!"
            # statement when conserved
            print(f"Mass conserved: {np.isclose(moles_0, moles_1, rtol=1e-4)}")

            # Calculate the real fluid factor with CORK EOS
            p2_p1_cork_real = ( (R * T) / (V2 - b) - (a) / (V2 * (V2 + b) * np.sqrt(T)) ) / \
                ( (R * T) / (Vm1 - b) - (a) / (Vm1 * (Vm1 + b) * np.sqrt(T)) )


        # Fluid pressure calculation
        # fluid pressure close to mean stress
        if self.fluid_pressure_mode == 'mean stress':
            # NOTE Modification here after revision 14.05.2025
            # hydro = mean_stress/vol_t0 * (vol_t0+(vol_new-vol_t0))
            # NOTE new equation here
            # ideal fluid
            hydro = litho * self.fluid_t1/(vol_t0-self.solid_t1)
            # real fluid factor
            # hydro = litho * self.fluid_t1/(vol_t0-self.solid_t1) * (0.0651 * self.fluid_t1/(vol_t0-self.solid_t1) + 0.936)
            hydro_CORK = litho * p2_p1_cork_real

            if self.reviewer_mode:
                # read weight fraction of fluid from phase data

                from scipy.optimize import root_scalar

                # From thermodynamic output

                # Calculate fluid and solid weights and densities for current step (t=1)
                fluid_weight_bound = self.hydrous_data_complete['df_H2O[g]'].iloc[:, -1].sum() / 1000 - self.hydrous_data_complete['df_H2O[g]'].loc[self.fluid_name_tag].iloc[-1] / 1000
                weight_total = self.phase_data.loc['wt[g]'].sum() / 1000
                if self.fluid_name_tag in self.phase_data.columns:
                    fluid_weight = self.phase_data[self.fluid_name_tag]['wt[g]'] / 1000
                    # Fluid density calculation
                    rho_f_1 = self.phase_data[self.fluid_name_tag]['density[g/ccm]'] * 1000  # g/ccm to kg/m^3
                    fluid_volume = self.phase_data[self.fluid_name_tag]['volume[ccm]'] / 1_000_000  # ccm to m^3
                else:
                    fluid_weight = 0.0
                    rho_f_1 = 0.0
                    fluid_volume = 0.0
                solid_weight = weight_total - fluid_weight

                # Volume calculations (convert ccm to m^3)
                solid_volume_total = self.phase_data.loc['volume[ccm]'].sum() / 1_000_000  # ccm to m^3
                solid_volume = solid_volume_total - fluid_volume

                # Solid density calculation
                rho_s_1 = solid_weight / solid_volume  # kg/m^3

                # Previous calculation step (t=0)
                fluid_weight_bound_before = self.hydrous_data_complete['df_H2O[g]'].iloc[:, -2].sum() / 1000 - self.hydrous_data_complete['df_H2O[g]'].loc[self.fluid_name_tag].iloc[-2] / 1000
                weight_total_before = self.phase_data_complete["df_wt[g]"].iloc[:, -2].sum() / 1000
                fluid_weight_before = self.phase_data_complete["df_wt[g]"].loc[self.fluid_name_tag].iloc[-2] / 1000
                solid_weight_before = weight_total_before - fluid_weight_before
                rho_f_0 = self.phase_data_complete["df_density[g/ccm]"].loc[self.fluid_name_tag].iloc[-2] * 1000  # g/ccm to kg/m^3

                # Volume calculations for previous step
                solid_volume_total_before = self.phase_data_complete["df_volume[ccm]"].iloc[:, -2].sum() / 1_000_000
                fluid_volume_before = self.phase_data_complete["df_volume[ccm]"].loc[self.fluid_name_tag].iloc[-2] / 1_000_000
                solid_volume_before = solid_volume_total_before - fluid_volume_before
                rho_s_0 = solid_weight_before / solid_volume_before

                # VOLUME CONSERVATION CONSTRAINT
                # If total volume remains constant, fluid must expand to fill space left by solid shrinkage
                V_total_const = solid_volume_total_before  # Constant total volume constraint

                # From thermodynamic output
                Xh_0 = fluid_weight_bound_before / solid_weight_before
                Xh_1 = fluid_weight_bound / solid_weight
                rho_bulk_0 = weight_total_before / solid_volume_total_before

                # Known pressures (Pa)
                Pf_0 = self.pressure_before * 1e5          # Bar → Pa
                Pf_1 = litho * 1e6         # litho is in MPa → convert to Pa

                # Linear interpolation function
                def interp(val_0, val_1, Pf_trial):
                    return val_0 + (val_1 - val_0) * (Pf_trial - Pf_0) / (Pf_1 - Pf_0)

                # Thermo property functions that vary with trial pressure
                rho_f_func = lambda Pf: interp(rho_f_0, rho_f_1, Pf)
                rho_s_func = lambda Pf: interp(rho_s_0, rho_s_1, Pf)
                Xh_func    = lambda Pf: interp(Xh_0, Xh_1, Pf)

                # Volume-conserving residual function
                def residual_volume_conserving(Pf_trial):
                    rho_f = rho_f_func(Pf_trial)
                    rho_s = rho_s_func(Pf_trial)
                    Xh = Xh_func(Pf_trial)

                    # Calculate solid volume at trial pressure (interpolated)
                    V_s_0 = solid_volume_before
                    V_s_1 = solid_volume
                    V_s_trial = interp(V_s_0, V_s_1, Pf_trial)

                    # Volume conservation: V_total = V_solid + V_fluid = constant
                    V_f_required = V_total_const - V_s_trial

                    # Calculate fluid mass at trial conditions
                    m_f_0 = fluid_weight_before
                    m_f_1 = fluid_weight
                    m_f_trial = interp(m_f_0, m_f_1, Pf_trial)

                    # Required fluid density to fit in available volume
                    rho_f_required = m_f_trial / V_f_required

                    # Residual: difference between thermodynamic density and required density
                    residual = rho_f - rho_f_required

                    # print(f"Pf = {Pf_trial/1e6:.1f} MPa → V_s = {V_s_trial:.9f} m³, V_f_req = {V_f_required:.9f} m³")
                    # print(f"    rho_f_thermo = {rho_f:.2f} kg/m³, rho_f_required = {rho_f_required:.2f} kg/m³, residual = {residual:.2f}")

                    return residual

                # Alternative: Bulk density conservation with volume constraint
                def residual_bulk_density_volume_conserving(Pf_trial):
                    rho_f = rho_f_func(Pf_trial)
                    rho_s = rho_s_func(Pf_trial)
                    Xh = Xh_func(Pf_trial)

                    # Calculate solid volume at trial pressure
                    V_s_0 = solid_volume_before
                    V_s_1 = solid_volume
                    V_s_trial = interp(V_s_0, V_s_1, Pf_trial)

                    # Volume conservation constraint
                    V_f_trial = V_total_const - V_s_trial

                    # Calculate masses at trial conditions
                    m_f_0 = fluid_weight_before
                    m_f_1 = fluid_weight
                    m_f_trial = interp(m_f_0, m_f_1, Pf_trial)

                    m_s_0 = solid_weight_before
                    m_s_1 = solid_weight
                    m_s_trial = interp(m_s_0, m_s_1, Pf_trial)

                    # Calculate bulk density with volume constraint
                    rho_bulk_trial = (m_f_trial + m_s_trial) / V_total_const

                    # Residual: bulk density should remain constant
                    residual = rho_bulk_trial - rho_bulk_0

                    # print(f"Pf = {Pf_trial/1e6:.1f} MPa → V_s = {V_s_trial:.9f}, V_f = {V_f_trial:.9f} m³")
                    # print(f"    m_f = {m_f_trial:.4f}, m_s = {m_s_trial:.4f} kg, rho_bulk = {rho_bulk_trial:.2f} kg/m³, residual = {residual:.2f}")

                    return residual

                # Root-finding wrapper
                def solve_pressure_volume_conserving(Pf_range=(1000e6, 5000e6)):
                    a, b = Pf_range
                    res_a = residual_func(a)
                    res_b = residual_func(b)
                    # print(f"Residual at {a/1e6:.1f} MPa: {res_a:.4f}")
                    # print(f"Residual at {b/1e6:.1f} MPa: {res_b:.4f}")
                    if res_a * res_b > 0:
                        raise ValueError("Residuals at both ends have the same sign — adjust the bracket!")

                    result = root_scalar(residual_func, bracket=Pf_range, method='brentq', xtol=1e-5)
                    if result.converged:
                        return result.root
                    else:
                        raise RuntimeError("Root-finding did not converge.")

                def molar_volume_cork(P_kbar, T_C, a_params, b):
                    """
                    Solve CORK EOS for molar volume (Vm) at given pressure (kbar) and temperature (°C).
                    Returns Vm in cm³/mol.
                    """
                    import numpy as np
                    from scipy.optimize import root_scalar

                    R = 83.14472  # cm³·bar/mol·K
                    T_K = T_C + 273.15
                    a, a1, a2, a3 = a_params

                    def cork_residual(Vm):
                        # CORK EOS residual: f(Vm) = LHS - RHS = 0
                        lhs = (R * T_K) / (Vm - b) - (a + a1 * T_K + a2 * T_K**2 + a3 * T_K**3) / \
                            (Vm * (Vm + b) * np.sqrt(T_K))
                        rhs = P_kbar * 1000  # Convert kbar → bar
                        return lhs - rhs

                    # Solve: Vm typically between 10 and 100 cm³/mol
                    sol = root_scalar(cork_residual, bracket=[10.0, 100.0], method='brentq')
                    if sol.converged:
                        return sol.root
                    else:
                        raise RuntimeError("CORK EOS did not converge for molar volume.")

                def residual_volume_conserving_cork(Pf_trial):
                    # Constants and CORK parameters  # in °C
                    a_params = (113.4, -0.22291, -3.8022e-4, 1.7791e-4)  # Holland & Powell (1991)
                    b = 1.465  # L/mol
                    M_H2O = 18.01528 / 1000  # kg/mol

                    # Convert Pf to kbar for CORK EOS
                    Pf_kbar = Pf_trial / 1e8  # Pa → kbar

                    try:
                        Vm_cm3mol = molar_volume_cork(Pf_kbar, T, a_params, b)  # cm³/mol
                        Vm_m3mol = Vm_cm3mol * 1e-6  # → m³/mol
                        rho_f_cork = M_H2O / Vm_m3mol  # kg/m³
                    except RuntimeError:
                        return np.nan

                    # Interpolate solid volume at trial pressure
                    V_s_trial = interp(solid_volume_before, solid_volume, Pf_trial)
                    V_f_trial = V_total_const - V_s_trial
                    if V_f_trial <= 0:
                        return np.nan

                    # Interpolate fluid mass
                    m_f_trial = interp(fluid_weight_before, fluid_weight, Pf_trial)

                    # Required fluid density to match volume constraint
                    rho_f_required = m_f_trial / V_f_trial

                    # Residual = (EOS-predicted density) - (required density for conservation)
                    residual = rho_f_cork - rho_f_required
                    return residual

                # Fixed bracket finding function
                def find_valid_bracket(start=10e6, stop=50000e6, step=1000e6):
                    Pf_values = np.arange(start, stop + step, step)
                    prev_Pf = Pf_values[0]
                    prev_res = residual_func(prev_Pf)

                    for Pf in Pf_values[1:]:
                        curr_res = residual_func(Pf)
                        if prev_res * curr_res < 0:
                            print(f"Valid bracket found: {prev_Pf/1e6:.0f} to {Pf/1e6:.0f} MPa")
                            return prev_Pf, Pf
                        prev_Pf, prev_res = Pf, curr_res

                    raise ValueError("No valid pressure bracket found — residual doesn't cross zero.")

                # Choose which residual function to use
                # Option 1: Volume-conserving fluid density constraint
                residual_func = residual_volume_conserving

                # Option 2: Bulk density conservation with volume constraint
                # residual_func = residual_bulk_density_volume_conserving

                # Solve and print for volume-conserving pressure Reviewer version
                if fluid_volume > 0.0:
                    print(f"\n=== SOLVING FOR VOLUME-CONSERVING PRESSURE ===")
                    try:
                        Pf_range = find_valid_bracket(start=10e6, stop=50000e7, step=1000e6)
                        Pf_corrected = solve_pressure_volume_conserving(Pf_range=Pf_range)
                        print(f"\nCorrected fluid pressure: {Pf_corrected / 1e6:.2f} MPa")
                        # Verify the solution
                        print(f"\n=== SOLUTION VERIFICATION ===")
                        final_residual = residual_func(Pf_corrected)
                        print(f"Final residual: {final_residual:.6f}")

                    except ValueError as e:
                        print(f"Error: {e}")
                        print("Trying alternative approach with bulk density conservation...")

                        # Switch to alternative residual function
                        residual_func = residual_bulk_density_volume_conserving
                        try:
                            Pf_range = find_valid_bracket(start=10e6, stop=50000e6, step=1000e6)
                            Pf_corrected = solve_pressure_volume_conserving(Pf_range=Pf_range)
                            print(f"\nCorrected fluid pressure: {Pf_corrected / 1e6:.2f} MPa")
                        except ValueError as e2:
                            print(f"Both approaches failed: {e2}")
                            print("Consider checking physical constraints or expanding pressure range.")
                            Pf_corrected = np.nan # Set to NaN if no solution found
                else:
                    print("Fluid volume is zero, cannot solve for pressure.")
                    Pf_corrected = np.nan







                """residual_func = residual_volume_conserving_cork
                #cork version
                # Solve and print for volume-conserving pressure Reviewer+CORK version
                if fluid_volume > 0.0:
                    print(f"\n=== SOLVING FOR VOLUME-CONSERVING PRESSURE ===")
                    try:
                        Pf_range = find_valid_bracket(start=10e6, stop=50000e7, step=1000e6)
                        Pf_corrected_cork = solve_pressure_volume_conserving(Pf_range=Pf_range)
                        print(f"\nCorrected fluid pressure: {Pf_corrected / 1e6:.2f} MPa")
                        # Verify the solution
                        print(f"\n=== SOLUTION VERIFICATION ===")
                        final_residual = residual_func(Pf_corrected_cork)
                        print(f"Final residual: {final_residual:.6f}")

                    except ValueError as e:
                        print(f"Error: {e}")
                        print("Trying alternative approach with bulk density conservation...")

                        # Switch to alternative residual function
                        residual_func = residual_bulk_density_volume_conserving
                        try:
                            Pf_range = find_valid_bracket(start=10e6, stop=50000e7, step=1000e6)
                            Pf_corrected_cork = solve_pressure_volume_conserving(Pf_range=Pf_range)
                            print(f"\nCorrected fluid pressure: {Pf_corrected_cork / 1e6:.2f} MPa")
                        except ValueError as e2:
                            print(f"Both approaches failed: {e2}")
                            print("Consider checking physical constraints or expanding pressure range.")
                            Pf_corrected_cork = np.nan # Set to NaN if no solution found
                else:
                    print("Fluid volume is zero, cannot solve for pressure.")
                    Pf_corrected_cork = np.nan"""



                hydro = Pf_corrected / 1e6  # Convert Pa to MPa

                print(f"{'Pressure type':<20} | {'Value (MPa)':>15}")
                print("-" * 38)
                print(f"{'Hydrostatic':<20} | {litho:15.3f}")
                print(f"{'CORK':<20} | {hydro_CORK:15.3f}")
                print(f"{'Reviewer':<20} | {Pf_corrected / 1e6:15.3f}")
                # print(f"{'Reviewer+CORK':<20} | {Pf_corrected_cork / 1e6:15.3f}")
            else:
                hydro = hydro_CORK
            hydro = hydro_CORK


            sig3 = litho - self.diff_stress/2
            sig1 = litho + self.diff_stress/2
            # hydro = mean_stress * (self.solid_t0+(self.solid_t1-self.solid_t0))/self.solid_t0 * (vol_new-vol_t0)/self.fluid_t0
            delta_p = litho - hydro

        elif self.fluid_pressure_mode == 'sig3':
            # NOTE Modification here after revision 14.05.2025
            # hydro = mean_stress/vol_t0 * (vol_t0+(vol_new-vol_t0))
            # NOTE new equation here
            hydro = sig3 * self.fluid_t1/(vol_t0-self.solid_t1)
            # hydro = sig3 * self.fluid_t1/(vol_t0-self.solid_t1)
            # real fluid factor
            hydro = sig3 * self.fluid_t1/(vol_t0-self.solid_t1) * (0.0651 * self.fluid_t1/(vol_t0-self.solid_t1) + 0.936)
            hydro_CORK = sig3 * p2_p1_cork_real

            hydro = hydro_CORK

            # hydro = mean_stress * (self.solid_t0+(self.solid_t1-self.solid_t0))/self.solid_t0 * (vol_new-vol_t0)/self.fluid_t0
            delta_p = sig3 - hydro

        elif self.fluid_pressure_mode == 'normal stress':
            # Duesterhoft 2019 method
            hydro = normal_stress/vol_t0 * (vol_t0+(vol_new-vol_t0))
            delta_p = normal_stress - hydro
        # if fluid_pressure_mode is a string but not mean_stress or normal_stress
        elif isinstance(self.fluid_pressure_mode, str):
            # NOTE - this is a very specific solution for the input of the fluid pressure
            # read the string and extract value "mean_stress-10" -> 10
            s = self.fluid_pressure_mode
            value = s.split("-")[1]
            new_pressure = mean_stress - float(value)
            hydro = new_pressure/vol_t0 * (vol_t0+(vol_new-vol_t0))
            delta_p = new_pressure - hydro
        else:
            hydro = mean_stress
            delta_p = mean_stress - hydro

        if self.fluid_name_tag not in self.phase_data.columns:
            print("No fluid phase data available")
            hydro = 0

        # Mohr circle arguments
        r = self.diff_stress/2      # radius of the circle
        center = sig1-r             # centre of the cricel
        pos = center - hydro        # position shift of centre due to fluid pressure

        # Coulomb failure envelope
        a = -1
        b = 1/self.friction
        c = -cohesion/self.friction

        # Critical fluid pressures
        crit_fluid_pressure = cohesion/self.friction + ((sig1+sig3)/2) - (
                (sig1-sig3)/2)*np.cos(2*theta) - ((sig1-sig3)/2/self.friction)*np.sin(2*theta)
        pf_crit_griffith = (8*self.tensile_strength*(sig1+sig3)-((sig1-sig3)**2))/(16*self.tensile_strength)

        # ##########################################
        # Failure envelope test
        thresh_value = self.extraction_connectivity
        # mcg_plot(cohesion, self.friction, self.diff_stress, self.shear_stress, sig1, hydro)
        if (vol_t0-self.solid_t1)/vol_t0 >= thresh_value:
            if self.diff_stress >= self.tensile_strength*5.66:

                print("Compressive shear failure test")
                # print(f"Diff.-stress is {self.diff_stress} and > than T*5.66")

                # Test possible extensional fracturing
                output, minimum = checkCollision_linear(a, b, c, x=pos, y=0, radius=r)
                # print(f"Maximum differential stress calculated as {minimum*2} MPa, but is {self.diff_stress}.")

                # Critical fluid pressure
                # print(f"Difference of fluid pressure and critical fluid pressure is:\n{crit_fluid_pressure-hydro} (Pf_crit - Pf)")

                if output is True:
                    print("---> Failure due to compressive shear.")
                    self.frac_respo = 3
                else:
                    print("---> No failure")

                # mcg_plot(cohesion, internal_friction, self.diff_stress, self.shear_stress, sig1, hydro)

            elif self.diff_stress < self.tensile_strength*5.66 and self.diff_stress > 4*self.tensile_strength:
                print("Extensional shear failure test")
                #print(f"Diff.-stress is {self.diff_stress} and < than T*5.66 and > T*4")
                output, minimum = checkCollision_curve(pos=pos, diff=self.diff_stress, tensile=self.tensile_strength)
                #print(f"Maximum differential stress calculated as {minimum*2} MPa, but is {self.diff_stress}.")
                #print(f"Difference of fluid pressure and critical fluid pressure is:\n{pf_crit_griffith-hydro} (Pf_crit - Pf)")

                if output is True:
                    print("---> Failure due to extensional shear.")
                    self.frac_respo = 2
                else:
                    print("---> No failure")

                # mcg_plot(cohesion, internal_friction, self.diff_stress, self.shear_stress, sig1, hydro)

            elif self.diff_stress <= 4*self.tensile_strength:
                print("Pure extensional failure test")
                # print(f"Diff.-stress is {self.diff_stress} and < T*4")

                if sig3-hydro <= -cohesion/2:
                    print("---> Failure due to pure extensional fail.")
                    self.frac_respo = 1
                else:
                    print("---> No failure")

                # mcg_plot(cohesion, internal_friction, self.diff_stress, self.shear_stress, sig1, hydro)

            else:
                print("Differential stress seems to have a problem. We decided to test the value. The value is:\n")
                print(self.diff_stress)
                self.frac_respo = 0
        else:
            print("Fluid-filled porosity to small to interconnect.")
            self.frac_respo = 0

        # NOTE Treshold sequence
        if self.frac_respo == 0:
            print("No mechanical failure detected.")
            print("Testing porosity and interconnectivity.")
            # print("P_f factor is {:.3f}".format(hydro/litho), "and porosity is {:.3f}".format(self.fluid_t1/(vol_t0-self.solid_t1)))
            # if hydro/litho > 0.9 and (vol_t0-self.solid_t1)/vol_t0 > 0.002:
            if (vol_t0-self.solid_t1)/vol_t0 >= self.extraction_threshold and self.fluid_t1 > 0:
                self.frac_respo = 5
            elif (vol_t0-self.solid_t1)/vol_t0 >= thresh_value:
                pass
            else:
                pass

        self.failure_dictionary = {
            "sigma 1":copy.deepcopy(sig1),
            "sigma 3":copy.deepcopy(sig3),
            "fluid pressure": copy.deepcopy(hydro),
            "delta p": copy.deepcopy(delta_p),
            "normal pressure": copy.deepcopy(normal_stress),
            "mean_stress": copy.deepcopy(mean_stress),
            "diff stress":(self.diff_stress),
            "fracture key": copy.deepcopy(self.frac_respo),
            "tensile strength": copy.deepcopy(self.tensile_strength),
            "shear stress": shear_stress,
            "friction coeff": copy.deepcopy(self.friction),
            "cohesion": (cohesion),
            "angle": copy.deepcopy(self.angle),
            "vol t0":copy.deepcopy(vol_t0),
            "solid t0":copy.deepcopy(self.solid_t0),
            "fluid t0":copy.deepcopy(self.fluid_t0),
            "vol t1":copy.deepcopy(vol_new),
            "solid t1":copy.deepcopy(self.solid_t1),
            "fluid t1":copy.deepcopy(self.fluid_t1),
            "radius": copy.deepcopy(r),
            "center": copy.deepcopy(center),
            "effective position": copy.deepcopy(pos),
            "critical fluid pressure normal": copy.deepcopy(crit_fluid_pressure),
            "critical fluid pressure griffith": copy.deepcopy(pf_crit_griffith),
            "pressure_reviewer": copy.deepcopy(Pf_corrected / 1e6 if self.reviewer_mode else np.nan),
            "Density fluid t0": rho_f_0 if 'rho_f_0' in locals() else np.nan,
            "Density fluid t1": rho_f_1 if 'rho_f_1' in locals() else np.nan,
            "Density solid t0": rho_s_0 if 'rho_s_0' in locals() else np.nan,
            "Density solid t1": rho_s_1 if 'rho_s_1' in locals() else np.nan,
            "Xh t0": Xh_0 if 'Xh_0' in locals() else np.nan,
            "Xh t1": Xh_1 if 'Xh_1' in locals() else np.nan
            }

        if self.frac_respo > 0:
            fracturing = True
        elif self.frac_respo == 0:
            fracturing = False
        else:
            print("The fracturing response failed. Thorsten has to work more here.")

        # !!!!
        #REVIEW - fluid is pressureized or expands into new volume
        # self.fluid_t1 = vol_t0 - self.solid_t1


        # Update system condition of fracturing
        self.fracture = fracturing
        print("End fracture modul")

    def factor_method(self):
        """
        Extraction of water/fluid follows the overstepping of a certain factor - calculated similar to Etheridge (2020)
        """

        print('\n')
        print('-----------------')
        print('solid_t0:{} solid_t-1:{} fluid_t0:{} fluid_t-1:{}'.format(self.solid_t0,
              self.solid_t1, self.fluid_t0, self.fluid_t1))
        print('-----------------')
        print('\n')

        diff_V_solids = (self.solid_t0 - self.solid_t1)
        diff_V_free_water = (self.fluid_t0 - self.fluid_t1)
        print(f"Diff in solid volume:\t\t {diff_V_solids}")
        print(f"Diff in fluid volume:\t\t {diff_V_free_water}")

        # Option 1: Calculating a fluid factor after Etheridge 2020
        fluid_factor_1 = (
            self.fluid_t0 + diff_V_free_water + diff_V_solids
        ) / self.fluid_t0
        print(f"Value of fluid factor 1:\t\t {fluid_factor_1}")
        # Option 2: Calculating a fluid factor 2 - in-house equation
        fluid_factor_2 = (self.solid_t0 + self.fluid_t0)/self.solid_t1
        print(f"Value of fluid factor 2: \t\t{fluid_factor_2}")
        # listing each fluid factor - saves it

        # Solution for new factor calculation (Etheridge):
        sys_vol_t0 = self.solid_t0 + self.fluid_t0
        sys_vol_t1 = self.solid_t1 + self.fluid_t1
        diff_sys = sys_vol_t0 - sys_vol_t1
        print(f"Diff in Sys:\t\t {diff_sys}")
        fluid_factor_ether = (self.pressure + self.pressure *
                              diff_sys/self.fluid_t0) / self.pressure
        print(f"Value of fluid factor Ether: \t\t{fluid_factor_ether}")

        # Alternative method - own idea:
        fluid_factor_mod = (self.fluid_t0/diff_V_solids*self.pressure /
                            10 - self.tensile_strength) / (self.pressure/10)
        print(f"Value of fluid factor Mod: \t\t{fluid_factor_mod}")

        # store selected fluid factor to be compared with extensional failure factor
        self.save_factor.append(abs(fluid_factor_ether))
        # self.save_factor.append(abs(fluid_factor_mod))


class Fluid_master():
    """
    Module to extract (modify the H content) the fluid in the system
        - different ways are possible
    """

    def __init__(self, phase_data, ext_data, temperature, new_fluid_V, sys_H, element_frame, st_fluid_post, fluid_name_tag):
        """
        Initialize the Tunorrad class.

        Args:
            phase_data (dict): Phase data for the fluid.
            ext_data (dict): External data.
            temperature (float): Temperature value.
            new_fluid_V (float): New fluid volume.
            sys_H (float): System H value.
            element_frame (DataFrame): Element frame.
            st_fluid_post (str): Fluid post.
            fluid_name_tag (str): Fluid name tag.
        """
        self.phase_data_fluid = phase_data
        self.ext_data = ext_data
        self.temp = temperature
        self.new_fluid_V = new_fluid_V
        self.sys_H = sys_H
        self.element_frame = element_frame.copy()
        self.st_fluid_post = st_fluid_post
        self.fluid_name_tag = fluid_name_tag

    def hydrogen_ext_all(self, extraction_percentage):
        """
        Recalculates hydrogen content due to fluid extraction.

        This method updates the hydrogen content in the system after fluid extraction.
        It subtracts the hydrogen content of the extracted fluid from the total system,
        and sets the remaining hydrogen content as the new bulk hydrogen content.
        The fluid hydrogen content is set to zero after this step.

        Parameters:
        - None

        Returns:
        - None
        """
        new_extraction = self.phase_data_fluid
        self.ext_data = pd.concat([self.ext_data, new_extraction], axis=1)
        self.ext_data = self.ext_data.rename(
            columns={self.fluid_name_tag: self.temp})

        # diminish double total: in dataframe
        if self.element_frame.columns[-1] == self.element_frame.columns[-2]:
            self.element_frame = self.element_frame.iloc[:, :-1]

        # settings some important values from element data frame
        extracted_fluid_vol = self.phase_data_fluid['volume[ccm]']
        fluid_H = self.element_frame.loc['H', self.fluid_name_tag]
        sys_H = self.element_frame.loc['H', 'total:']
        fluid_O = self.element_frame.loc['O', self.fluid_name_tag]
        sys_O = self.element_frame.loc['O','total:']

        if extraction_percentage == 1.0:
            # remaining H and O in bulk
            total_hydrogen = sys_H - fluid_H
            total_oxygen = sys_O - fluid_O

            # IMPORTANT step - substracting H and O of fluid from total system - total system will be new bulk
            self.element_frame.loc['H', 'total:'] = total_hydrogen
            self.element_frame.loc['O', 'total:'] = total_oxygen

            # set fluid O and H to zero after this step
            self.element_frame.loc['O', self.fluid_name_tag] = 0.0
            self.element_frame.loc['H', self.fluid_name_tag] = 0.0
            self.st_fluid_post[-1] = self.st_fluid_post[-1] - extracted_fluid_vol
            # print("=====================================")
            print("======= Fluid fractionation 100% ========") 
        else:
            if new_extraction['vol%']/100 < extraction_percentage:
                fraction = 0.0
            else:
                fraction = (new_extraction['vol%']/100 - extraction_percentage) / (new_extraction['vol%']/100)

            # remaining H and O in bulk
            total_hydrogen = sys_H - fluid_H * fraction
            total_oxygen = sys_O - fluid_O * fraction

            # IMPORTANT step - substracting H and O of fluid from total system - total system will be new bulk
            self.element_frame.loc['H', 'total:'] = total_hydrogen
            self.element_frame.loc['O', 'total:'] = total_oxygen

            # set fluid O and H to zero after this step
            self.element_frame.loc['O', self.fluid_name_tag] = fluid_O - fluid_O * fraction
            self.element_frame.loc['H', self.fluid_name_tag] = fluid_H - fluid_H * fraction
            self.st_fluid_post[-1] = self.st_fluid_post[-1] - extracted_fluid_vol * fraction
            # print("=====================================")
            print(f"======= Fluid fractionation {fraction*100}% ========")

        print("\t    \N{BLACK DROPLET}")
        print("\t   \N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}")
        print("\t \N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}")
        print("\t\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}")
        print("\t\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}")
        print("\t \N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}")
        print("\t   \N{BLACK DROPLET}\N{BLACK DROPLET}\N{BLACK DROPLET}")
        # print("=====================================")

    def hydrogen_partial_ext(self, fraction_value):
        new_extraction = self.phase_data_fluid
        self.ext_data = pd.concat([self.ext_data, new_extraction], axis=1)
        self.ext_data = self.ext_data.rename(
            columns={self.fluid_name_tag: self.temp})

        # diminish double total: in dataframe
        if self.element_frame.columns[-1] == self.element_frame.columns[-2]:
            self.element_frame = self.element_frame.iloc[:, :-1]

        # settings some important values from element data frame
        extracted_fluid_vol = self.phase_data_fluid['volume[ccm]']
        fluid_H = self.element_frame.loc['H', self.fluid_name_tag]
        sys_H = self.element_frame.loc['H', 'total:']
        fluid_O = self.element_frame.loc['O', self.fluid_name_tag]
        sys_O = self.element_frame.loc['O','total:']

        # remaining H and O in bulk
        fractionation_value = (new_extraction['vol%']/100 - fraction_value) / (new_extraction['vol%']/100)
        total_hydrogen = sys_H - fluid_H * fractionation_value
        total_oxygen = sys_O - fluid_O * fractionation_value

        # IMPORTANT step - substracting H and O of fluid from total system - total system will be new bulk
        self.element_frame.loc['H', 'total:'] = total_hydrogen
        self.element_frame.loc['O', 'total:'] = total_oxygen

        # set fluid O and H to zero after this step
        self.element_frame.loc['O', self.fluid_name_tag] = fluid_O - fluid_O * fractionation_value
        self.element_frame.loc['H', self.fluid_name_tag] = fluid_H - fluid_H * fractionation_value
        self.st_fluid_post[-1] = self.st_fluid_post[-1] - extracted_fluid_vol * fractionation_value

        print(f"======= Fluid fractionation {fractionation_value*100}% ========")

class System_status:
    """
    Compiling the system conditions for the P-T-t path. Ready for plotting and investigations.
    """

    def __init__(self, df_var_dictionary, df_h2o_content_dic, df_element_total, st_elements):
        """
        Initialize the Tunorrad class.

        Args:
            df_var_dictionary (dict): A dictionary containing thermodynamic data.
            df_h2o_content_dic (dict): A dictionary containing hydrate data.
            df_element_total (DataFrame): A DataFrame containing total element data.
            st_elements (str): A string representing the elements.

        Returns:
            None
        """
        self.therm_data = df_var_dictionary
        self.hydrate_data = df_h2o_content_dic
        self.df_element_total = df_element_total
        self.st_elements = st_elements
        self.sys_dat = 0

    def formatting_data(self, temperatures, st_solid, st_fluid_before, st_fluid_after, extracted_fluid_data, line=0):
        """
        Re-organizing the Dataframes for finalizing the data

        Args:
            temperatures (list): all the temperature values
            st_solid (list): list of the solid volumina over P-T
            st_fluid_before (list): list of the fluid volumina before extraction over P-T
            st_fluid_after ([type]): list of the fluid volumina after extraction over P-T
            extracted_fluid_data (Dataframe): Data of the fluid
        """

        # edit 220406 iterating number is now index instead of temperature
        # based on len of temperature

        for key in self.therm_data.keys():
            self.therm_data[key] = self.therm_data[key].T
            # self.therm_data[key].index = temperatures
            self.therm_data[key].index = line
        for key in self.hydrate_data.keys():
            self.hydrate_data[key] = self.hydrate_data[key].T
            # self.hydrate_data[key].index = temperatures
            if len(line) > len(self.hydrate_data[key].index):
                # FIXME line list of hydrated phases is wrong
                print("line list of hydrated phases is wrong")
                self.hydrate_data[key].index = np.arange(0, len(self.hydrate_data[key].index),1)
            else:
                self.hydrate_data[key].index = line

        self.df_element_total = self.df_element_total.T

        self.st_elements.columns = line

        # system volumina
        system_vol_pre = np.array(st_solid) + np.array(st_fluid_before)
        system_vol_post = np.array(st_solid) + np.array(st_fluid_after)

        # system weights
        system_weight_pre = self.therm_data['df_wt[g]'].T.sum()
        chache1 = self.therm_data['df_wt[g]'].T.sum()
        if 'wt[g]' in extracted_fluid_data.index:
            chache2 = extracted_fluid_data.loc['wt[g]']
            chache1[list(chache2.index)] = chache1[list(chache2.index)] - chache2

        system_weight_post = chache1

        # system densities
        system_density_pre = system_weight_pre/system_vol_pre
        system_density_post = system_weight_post/system_vol_post

        self.sys_dat = {'system_vol_pre': system_vol_pre, 'system_vol_post': system_vol_post,
                        'system_weight_pre': system_weight_pre, 'system_weight_post': system_weight_post,
                        'system_density_pre': system_density_pre, 'system_density_post': system_density_post}


class Isotope_calc():
    """
    Module to calculate the oxygen fractionation between the stable phases for each P-T step
    """

    def __init__(self, data_mol_from_frame, eq_data, element_data, oxygen_signature):
        """
        Initilization of datasets, frames and values

        Args:
            data_mol_from_frame (Dataframe): Data with information about the moles for each phase
            eq_data (Dataframe): Information about the phases and sthe solid solutions
            element_data (Dataframe): Distribution of oxygen in the stable phases
            oxygen_signature (int): Bulk rock oxygen value at given condition
                - changes over P-T and especially extraction
        """
        # open to write sth
        self.oxygen_dic = {}
        # stable phase data (df_var_dic/df_N)
        self.phase_data = data_mol_from_frame
        # minimization.sol_sol_base
        self.eq_data = eq_data
        self.start_bulk = oxygen_signature
        self.element_oxy = element_data.loc['O']

    def frac_oxygen(self, temperature):
        """
        Checks for the stable phases and including solid solution.
        Recalculates the fractions.
        Reads fractionation factors.
        Then applies the fractionation of oxygen between all phases by using a linear equation system.
        Saves results to "stable phases" and "oxygen values"

        Args:
            temperature (int): Temperature condition for P-T step
        """
        phase_translation = mineral_translation()

        stable_phases = list(self.element_oxy.index[:-1])

        pos_index_list = []
        enabled_phase = []
        layered_name = []
        for i, phase in enumerate(stable_phases):
            # print(phase)
            try:
                # Catching stable phases and checking with name pattern and solid-solution entry
                if phase in phase_translation['Theriak_phases']:
                    # position from already clean naming pattern
                    pos_index = phase_translation['Theriak_phases'].index(
                        phase)
                    pos_index_list.append(pos_index)
                    # collect detected phase present in model
                    enabled_phase.append(phase)
                    # save enabled position
                    layered_name.append(True)
                # In the case the name looks like "phase_spec"
                elif '_' in phase:
                    dunder_pos = phase.index('_')
                    phaseOK = phase[0:dunder_pos]
                    pos_index = phase_translation['Theriak_phases'].index(
                        phaseOK)
                    pos_index_list.append(pos_index)
                    # collect detected phase present in model
                    enabled_phase.append(phase)
                    layered_name.append(True)
                else:
                    layered_name.append(False)
            except ValueError:
                print(
                    f"---ERROR:---Searching for stable phase in translation failed. No such phase: {phase} found!")
                for name in self.eq_data['Memb_Names'][i]:
                    if name in phase_translation['Theriak_phases']:
                        print(f"---But {name} is found instead.")
                        pos_index = phase_translation['Theriak_phases'].index(
                            name)
                        pos_index_list.append(pos_index)
                        enabled_phase.append(name)
                        layered_name.append(True)
                else:
                    layered_name.append(False)

        # chooses detected phase present in model and stable at P-T
        stable_phases = enabled_phase

        # disabled phases (like "Si", ...)
        all_phases = self.eq_data['Name']
        excluded_phases = []
        excluded_phases_index = []
        for phase in all_phases:
            if phase not in stable_phases:
                excluded_phases.append(phase)
                # get index of excluded phases
                pos_index = self.eq_data['Name'].index(phase)
                excluded_phases_index.append(pos_index)

        # deep copy of moles for the stable phases to avoid changing the original data (some values might be excluded)
        moles_phases_copy = copy.deepcopy(self.eq_data['Moles'])
        moles_phases_copy = np.delete(moles_phases_copy,excluded_phases_index,0)
        if len(moles_phases_copy) != len(layered_name):
            layered_name = list(np.delete(layered_name, excluded_phases_index, 0))

        phases_X = []
        database_names = []
        ii = 0
        # collects the translation names and the composition of all the sable phases and their solid solution members
        for i, pos in enumerate(pos_index_list):
            # testing if entry for position is a solid solution
            if phase_translation['SSolution_Check'][pos] is True:

                # selecting composition for solid solution members
                solsol_frac_temp = []
                solsol_memb = phase_translation['SSolution_set'][pos][1]
                # pos_index = self.eq_data['Name'].index(phase)
                j = 0
                # test indexing i using layered_name - False marks an offset
                # (if a phase such as "CRPH_mcar" is not in database, there is an offset)
                if layered_name[i] is False:
                    ii += 1
                    # print(f"=====\nOffset in oxygen database reader no\ni is {i} and name {self.eq_data['Memb_Names'][i]}\nfollow name is {self.eq_data['Memb_Names'][i+1]}=====")
                # checks for solid-solution endmembers
                for phase in solsol_memb:
                    if phase in self.eq_data['Memb_Names'][i+ii]:
                        memb_ind = self.eq_data['Memb_Names'][i +
                                                              ii].index(phase)
                        solsol_frac_temp.append(
                            self.eq_data['Comp'][i+ii][memb_ind])
                        j += 1
                    else:
                        print(
                            f"=4== OxyFracModul reports: \n\tEndmembers are {self.eq_data['Memb_Names'][i+ii]} -- {phase} not stable.")
                phases_X.append(solsol_frac_temp)

                # appends solid solution list to nested list - depending on which solidsol member is stable
                #   (maybe just three of four garnet endmembers are stable)
                database_names.append(
                    phase_translation['SSolution_set'][pos][0][0:j])

            if phase_translation['SSolution_Check'][pos] is False:
                # get the name of the phase from the oxygen database
                database_names.append(
                    phase_translation['DB_phases'][pos].split())
                # append the moles of phase X to the list
                phases_X.append(float(moles_phases_copy[i+ii]))

                # selecting composition for pure phases
                # mono = phase_translation['Theriak_phases'][pos]

        # undo nested list into simple list
        # phases_frac = list(flatten(phases_X))

        # read oxygen isotope frac database

        main_folder = Path.cwd()
        file_to_open = main_folder / "DataFiles" / "DO18db2.0.3.dat"

        df = pd.read_csv(file_to_open, sep=" ", header=None)
        df.index = list(df[0])
        df.columns = ['Name', 'A', 'B', 'C', 'a', 'b', 'c', 'R2']
        df.drop('Name', axis=1)
        df_ABC = df.drop(['a', 'b', 'c', 'R2'], axis=1)
        df_ABC = df_ABC.T

        # frame matrix for thermodynamic fractionation calculation
        ma = []
        TK = temperature + 273.15

        phase = list(flatten(database_names))[0]
        norm_ABC = np.array(
            [df_ABC[phase].iloc[1], df_ABC[phase].iloc[2], df_ABC[phase].iloc[3]])

        for phase in list(flatten(database_names)):
            a = df_ABC[phase].iloc[1] - norm_ABC[0]
            b = df_ABC[phase].iloc[2] - norm_ABC[1]
            c = df_ABC[phase].iloc[3] - norm_ABC[2]
            calc = a*10**6/TK**2 + b*10**3/TK + c
            ma.append(calc)

        # Using themo frac and oxygen data for endmembers to calculate isotope values for phases
        # Sum part of equation (7) from Vho 2020

        # preparing fractionation for LES (linear equation system)
        ma = ma[1:]
        ma.append(0)
        ma.append(self.start_bulk)
        # print(f"Bulk oxygen is {self.start_bulk}")

        # deep copy moles of oxygen in the phase to not modifiy data
        oxygen_moles = copy.deepcopy(np.asarray(self.element_oxy))
        # delete the excluded values from phases such as ("SI", ...)
        oxygen_moles = np.delete(oxygen_moles, excluded_phases_index, 0)

        # collect the fraction of stable phases depending on a solid solution or a monotopic phase
        phases_fraction = []
        moles_frac = []
        for i, phase in enumerate(phases_X):
            pos = pos_index_list[i]
            if phase_translation['SSolution_Check'][pos] is True:
                for member in phase:
                    frac_temp = member / sum(np.array(phase))
                    phases_fraction.append(frac_temp)
                    mole_temp = frac_temp*oxygen_moles[i]
                    moles_frac.append(mole_temp)
            else:
                phases_fraction.append(phase)
                mole_temp = oxygen_moles[i]
                moles_frac.append(mole_temp)

        # creating fractionation matrix for LES
        mat = np.ones(len(list(flatten(database_names)))-1)*-1
        d = np.diag(mat)
        n, m = d.shape
        n0 = np.ones((n, 1))
        d = np.hstack((n0, d))
        d = np.hstack((d, np.zeros((n, 1))))
        # moles_frac[-2] = 0
        moles_frac = list(np.nan_to_num(np.asarray(moles_frac), 0))
        if len(list(flatten(database_names))) == len(moles_frac):
            moles_frac.append(-sum(moles_frac))
        else:
            moles_frac.append(-sum(moles_frac[:-1]))
        d = np.vstack((d, moles_frac))
        d = np.vstack((d, np.zeros(len(moles_frac))))
        d[-1, -1] = 1

        # solving equation system
        X = scipy.linalg.solve(d, ma)

        oxygen_values = []
        count = 0
        # Summing up the endmembers oxygen isotope comp into one value for all stable phases (check with stable_phases)
        for i, phase in enumerate(database_names):
            pos = pos_index_list[i]
            # in the case of more than one endmember for a phase
            if phase_translation['SSolution_Check'][pos] is True:
                # print(phase)
                length = len(database_names[i])
                # print('X = {} fract = {}'.format(X[count:count+length], phases_fraction[count:count+length]))
                result = sum(np.array(X[count:count+length])
                             * np.array(phases_fraction[count:count+length]))
                oxygen_values.append(result)
                count += length
            # case if phase is not a solid solution
            else:
                # print('X = {} fract = {}'.format(X[count], phases_fraction[count]))
                result = X[count]
                oxygen_values.append(result)
                count += 1
            # print('phase:{} Oxy-Val:{}'.format(phase, result))
        self.oxygen_dic['Phases'] = stable_phases
        self.oxygen_dic['delta_O'] = oxygen_values
        #print("Oxygen fractionation calculation done.")


def phases_translated(database, phases):
    """
    Returns a list of phase names and their corresponding colors based on the given database and phases.

    Parameters:
    - database (str): The name of the database.
    - phases (list): A list of phase names.

    Returns:
    - phase_set (list): A list of phase names.
    - color_set (list): A list of RGB color tuples corresponding to each phase.

    """

    # Translation file - database to XMT names
    script_folder = Path(__file__).parent.absolute()
    dot_indi = database.index('.')
    file_to_open = script_folder / "DataFiles" / \
        f"MINERAL_NAMES_{database[:dot_indi]}_to_TE.txt"
    min_translation = pd.read_csv(file_to_open, delimiter='\t')

    # Iterating through phases 
    phase_set = []
    phase_original = []
    z = 0
    for mini in phases:
        if mini == 'fluid' or mini == 'H2O.liq' or mini == 'Liqtc6_H2Ol' or mini == 'H2O' or mini == 'water.fluid':
            phase_set.append('Fluid')
            phase_original.append(mini)
        elif mini in list(min_translation['Database-name']):
            phase_set.append(min_translation[
                min_translation['Database-name']== mini]['XMT-name'].iloc[0])
            phase_original.append(mini)
        elif '_' in mini:
            mini_inid = mini.index('_')
            base_name = mini[:mini_inid]
            if base_name in list(min_translation['Database-name']):
                phase_set.append(min_translation[min_translation['Database-name']
                                           == base_name]['XMT-name'].iloc[0])
                phase_original.append(mini)

    return phase_set, phase_original

# Function to check if a variable is np.float64 and NaN
def is_float64_nan(value):
    return isinstance(value, np.float64) and np.isnan(value)

class TraceElementDistribution():

    def __init__(self, data_mol_from_frame, eq_data, element_data, bulk_trace_elements, database):
        """
        Initilization of datasets, frames and values

        Args:
            data_mol_from_frame (Dataframe): Data with information about the moles for each phase
            eq_data (Dataframe): Information about the phases and sthe solid solutions
            element_data (Dataframe): Element distribution in moles for the stable phases
            bulk_trace_elements (int): Bulk rock trace element value at input condition
        """

        self.phase_data = data_mol_from_frame
        self.eq_data = eq_data
        self.element_data = element_data
        self.start_bulk = bulk_trace_elements
        self.distribution_coefficients = read_trace_element_content()
        self.database = database

        if 'tc55' in database:
            self.phase_set, self.phase_original = phases_translated('tc55.txt', self.phase_data.index)
        else:
            self.phase_set, self.phase_original = phases_translated(self.database, self.phase_data.index)

    def distribute_tracers(self, temperature, pressure, iteration):
        """
        runs the distribution of trace elements between the stable phases
        """

        # name list is self.distribution_coefficients.index but clip /Grt from each name
        name_list = [name[:-4] for name in self.distribution_coefficients.index]
        name_list.append('Garnet')
        phase_list = self.phase_data.index
        element_list = self.distribution_coefficients.columns

        # Assuming self.distribution_coefficients is a DataFrame
        # add 'Grt' to distribution_coefficients with array of ones
        distribution_coefficients = self.distribution_coefficients
        distribution_coefficients.loc['Garnet'] = np.ones(len(distribution_coefficients.columns))
        distribution_coefficients.index = name_list

        # collecting stable phases moles and names
        assembled_distribution_coefficients = pd.DataFrame([])
        for name in name_list:
            if name in self.phase_set:
                assembled_distribution_coefficients = pd.concat([
                        assembled_distribution_coefficients, 
                        distribution_coefficients.loc[name]], axis=1)

        # final assembled distribution coefficients based on predicted stable phases
        assembled_distribution_coefficients = assembled_distribution_coefficients.T
        # update name_list to the names of the phases that are present
        name_list = assembled_distribution_coefficients.index
        
        # Calculating mineral/matrix distribution coefficients 
        # (D_min-grt_e / sum(D_min-grt_e) = D_min-mtrx_e)
        # -----------------------------------------------------------

        # Convert the DataFrame to a numpy array for faster operations
        coeff_array = assembled_distribution_coefficients.values

        # Calculate the sum of each column
        column_sums = np.sum(coeff_array, axis=0)

        # Calculate the coefficient values
        _m_min_matrix_coeff = coeff_array / column_sums.T
        #TODO - check if this is correct if not all mineral phases are present

        _m_min_matrix_coeff_df = pd.DataFrame(
            _m_min_matrix_coeff,
            index=name_list,
            columns=element_list)

        # Calculating k-matrix coefficients
        # (c_bulk_e / sum(moles_min *D_min-mtrx_e) = K_e)
        # -----------------------------------------------------------
        # create empty dataframe to be filled with k values
        mass_balanced_trace_elements = pd.DataFrame([])
        selected_phase_moles = []
        assembled_matrix_coeff = pd.DataFrame([])
        selected_phase_name_original = []
        #assembling the trace element distribution matrix
        for i, name in enumerate(self.phase_set):
            if name in name_list:
                if is_float64_nan(self.phase_data.iloc[i,-1]):
                    pass
                else:
                    selected_phase_moles.append(self.phase_data.iloc[i,-1])
                    selected_phase_name_original.append(self.phase_original[i])
                    assembled_matrix_coeff = pd.concat([
                        assembled_matrix_coeff, _m_min_matrix_coeff_df.loc[name]], axis=1)

        assembled_matrix_coeff = assembled_matrix_coeff.T
        selected_phase_moles = np.array(selected_phase_moles)

        # Calculating the mass balanced distribution coefficients
        # mass_balanced_trace_elements = selected_phase_moles.values * assembled_matrix_coeff.values
        mass_balanced_trace_elements = selected_phase_moles[:, np.newaxis] * assembled_matrix_coeff.values
        # mass_balanced_trace_elements = pd.DataFrame(selected_phase_moles.values * assembled_matrix_coeff.values)
        # mass_balanced_trace_elements.columns = element_list
        # mass_balanced_trace_elements.index = selected_phase_moles.index
        
        # Calculate the sum of each column
        column_sums = np.sum(mass_balanced_trace_elements, axis=0)
        
        # Calculate the coefficient values by dividing the bulk values by the sum of the column
        if len(column_sums) == 0:
            
            column_sums = np.nan
            k_factor = self.start_bulk.values[0] / column_sums

            content =  k_factor
            selected_phase_name_original = 'None'
            content_df = pd.DataFrame(content)
            content_df.index = element_list
            #content_df.columns = selected_phase_name_original
        else:
            k_factor = self.start_bulk.values[0] / column_sums

            # Calculating the content of element e in mineral
            # content = K_e * moles_min * D_min-mtrx_e
            # -----------------------------------------------------------

            #content =  k_factor * selected_phase_moles.values  * assembled_matrix_coeff.values
            content =  k_factor * selected_phase_moles[:, np.newaxis]  * assembled_matrix_coeff.values
            
            #build dataframe of content
            #content_df = pd.DataFrame(content, index=selected_phase_moles.index, columns=element_list)
            #content_df = pd.DataFrame(content, index=selected_phase_moles.index, columns=element_list)
            content_df = pd.DataFrame(content)
            content_df.columns = element_list
            content_df.index = selected_phase_name_original
        
            """
            norming = np.array([
                0.3670, 0.9570, 0.1370, 0.7110, 0.2310, 0.0870, 0.3060, 
                0.0580, 0.3810, 0.0851, 0.2490, 0.0356, 0.2480, 0.0381])
            test = content_df/norming

            # test plot in log scale on the y axis, y_label is REE/Chondrite and x_label is elements
            for phase in test.index:
                plt.plot(test.loc[phase], 'D--' , label=phase)
            plt.yscale('log')
            plt.xlabel('Elements')
            plt.ylabel('REE/Chondrite')
            plt.legend()
            """
        return content_df


class Garnet_recalc():
    """
    Class for performing recalculation of garnets.

    Attributes:
        theriak (str): The path to the Theriak software.
        dataclass (list): List of garnet data objects.
        temperature (float): The temperature in Kelvin.
        pressure (float): The pressure in bars.
        recalc_volume (float): The recalculated volume of garnets.
        recalc_weight (float): The recalculated weight of garnets.
    """

    def __init__(self, theriak, dataclass, temperature, pressure):
        """
        Initialize a Tunorrad object.

        Args:
            theriak (str): The path to the Theriak executable.
            dataclass: The data class for the mineral.
            temperature (float): The temperature in Kelvin.
            pressure (float): The pressure in bars.
        """
        self.mineral = dataclass
        self.recalc_volume = 0
        self.recalc_weight = 0
        self.temperature = temperature
        self.pressure = pressure
        self.theriak_path = theriak

    def recalculation_of_garnets(self, database, garnet_name):
        """
        Recalculates the volume and weight of garnets based on the provided data.

        Returns:
            None
        """
        for garnet in self.mineral:
            vals = garnet.elements[0]
            index = garnet.elements[1]
            relements = pd.DataFrame(vals, index=index)

            # create the bulk from element entry - normalized to 1 mol for theriak
            # forward - volume needs to be back-converted to the moles of the shell (garnet.moles)
            bulk = garnet_bulk_from_dataframe(relements, garnet.moles, garnet.volPmole, garnet.volume)

            # FIXME static database tc55 for garnet recalculation
            db = f"{database}    {garnet_name}"
            Data_dic, g_sys, pot_frame = read_theriak(
                self.theriak_path,
                database=db, temperature=self.temperature,
                pressure=self.pressure, whole_rock=bulk,
                theriak_input_rock_before=False)
            grt = Data_dic['df_Vol_Dens'].columns[0]
            phase_list = list(Data_dic['df_Vol_Dens'].columns)

            # check for couble garnet
            double_grt_check = []
            selected_garnet_name = False
            for item in phase_list:
                if garnet_name in item:
                    selected_garnet_name = item
                    double_grt_check.append(item)

            if selected_garnet_name is False:
                print("WARNING: Garnet increment no longer stable. Press ESC to continue.")
                print(f"Bulk is: {bulk}")
                # keyboard.wait('esc')
                v = 0
                g = 0
            else:
                # backward - Volume of the garnet shell to be added (back-normlaized from 1 mole to x-mol (garnet.moles))
                rec_mole = garnet.volume/garnet.volPmole
                v = Data_dic['df_Vol_Dens'][grt]['volume/mol']*rec_mole
                g = Data_dic['df_Vol_Dens'][grt]['wt/mol']*rec_mole

            self.recalc_volume += v
            self.recalc_weight += g


class Metastable_phase_recalculator():
    """
    Class for performing recalculation of garnets.

    Attributes:
        theriak (str): The path to the Theriak software.
        dataclass (list): List of garnet data objects.
        temperature (float): The temperature in Kelvin.
        pressure (float): The pressure in bars.
        recalc_volume (float): The recalculated volume of garnets.
        recalc_weight (float): The recalculated weight of garnets.
    """

    def __init__(self, theriak, dataclass, temperature, pressure):
        """
        Initialize a Tunorrad object.

        Args:
            theriak (str): The path to the Theriak executable.
            dataclass: The data class for the mineral.
            temperature (float): The temperature in Kelvin.
            pressure (float): The pressure in bars.
        """
        self.mineral = dataclass
        self.recalc_volume = 0
        self.recalc_weight = 0
        self.temperature = temperature
        self.pressure = pressure
        self.theriak_path = theriak

    def recalculation_of_garnets(self, database, garnet_name):
        """
        Recalculates the volume and weight of garnets based on the provided data.

        Returns:
            None
        """
        for garnet in self.mineral:
            vals = garnet.elements[0]
            index = garnet.elements[1]
            relements = pd.DataFrame(vals, index=index)

            # create the bulk from element entry - normalized to 1 mol for theriak
            # forward - volume needs to be back-converted to the moles of the shell (garnet.moles)
            bulk = garnet_bulk_from_dataframe(relements, garnet.moles, garnet.volPmole, garnet.volume)

            # FIXME static database tc55 for garnet recalculation
            db = f"{database}    {garnet_name}"
            Data_dic, g_sys, pot_frame = read_theriak(
                self.theriak_path,
                database=db, temperature=self.temperature,
                pressure=self.pressure, whole_rock=bulk,
                theriak_input_rock_before=False)
            grt = Data_dic['df_Vol_Dens'].columns[0]
            phase_list = list(Data_dic['df_Vol_Dens'].columns)

            # check for couble garnet
            double_grt_check = []
            selected_garnet_name = False
            for item in phase_list:
                if garnet_name in item:
                    selected_garnet_name = item
                    double_grt_check.append(item)

            if selected_garnet_name is False:
                print("WARNING: Garnet increment no longer stable. Press ESC to continue.")
                print(f"Bulk is: {bulk}")
                # keyboard.wait('esc')
                v = 0
                g = 0
            else:
                # backward - Volume of the garnet shell to be added (back-normlaized from 1 mole to x-mol (garnet.moles))
                rec_mole = garnet.volume/garnet.volPmole
                v = Data_dic['df_Vol_Dens'][grt]['volume/mol']*rec_mole
                g = Data_dic['df_Vol_Dens'][grt]['wt/mol']*rec_mole

            self.recalc_volume += v
            self.recalc_weight += g





