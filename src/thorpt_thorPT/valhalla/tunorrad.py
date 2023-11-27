
"""
Written by
Thorsten Markmann
thorsten.markmann@geo.unibe.ch
status: 17.02.2023
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


@dataclass
class Phaseclass:
    comps: List[Phasecomp]

# collision of mohr circle with shear failure envelope


def checkCollision(a, b, c, x, y, radius):
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

    return output

def checkCollision_linear(a, b, c, x, y, radius):
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
    bulk = bulk
    el_list = bulk.keys()
    if 'H' not in el_list:
        bulk['H'] = 0
    if 'C' not in el_list:
        bulk['C'] = 0
    if 'MN' not in el_list:
        bulk['MN'] = 0
    if 'CA' not in el_list:
        bulk['CA'] = 0
    if 'TI' not in el_list:
        bulk['TI'] = 0
    if 'MG' not in el_list:
        bulk['MG'] = 0
    if 'NA' not in el_list:
        bulk['NA'] = 0
    if 'K' not in el_list:
        bulk['K'] = 0

    new_bulk = (
        f"SI({bulk['SI']})AL({bulk['AL']})FE({bulk['FE']})"
        f"MN({bulk['MN']})MG({bulk['MG']})CA({bulk['CA']})"
        f"NA({bulk['NA']})TI({bulk['TI']})K({bulk['K']})"
        f"H({bulk['H']})C({bulk['C']})O(?)O(0)    * CalculatedBulk"
    )
    return new_bulk


def first_Occurrence_char(char, string):
    """finds position of first occurence

    Args:
        char ([string]): [can be a letter or anything]
        string ([string]): [can be a word]

    Returns:
        [type]: [description]
    """
    # print(string)
    for num, item in enumerate(string):
        if(string[num] == char):
            return num
    return 'no'


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


def flatten(lis):
    """
    Cleaning nested lists to on level

    Args:
        lis ([nested list]): [any nested list]

    Yields:
        [list]: [single level list]
    """
    for item in lis:
        if isinstance(item, Iterable) and not isinstance(item, str):
            for x in flatten(item):
                yield x
        else:
            yield item


def mineral_translation():
    """
    read DB_database in DataFiles file and stores to list
    - is the translation between theriak nad the DBOxygen dataset

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


def decode_lines(line_number, line_file, number=11, one_element_row=True):
    """
    Decodes the output from theriak from a line into dictionary.
        - Especially the whitespaces are annoying

    Args:
        line_number (int): Line number of an output or txt file read as a list -
        def "keyword search" find the line for specific keywords
        line_file (type): declares the variable name in which the file
            (eg. the theriak output) is stored
        number (int, optional): Variable to handle the lines with less or more
        than 10 entries (Theriak switches to the next line e.g.,
        if more than 10 elements are in the system). Defaults to 10.
    """
    if one_element_row is False:
        line_number = line_number + 1

    temp_list = line_file[line_number]
    list_decrypt = []
    temp_dictionary = {}
    # print(f"1) {temp_list}\n")
    temp_list = temp_list.split(' ')
    # print(f"2) {temp_list}\n")
    for word in temp_list:
        list_decrypt.append(word)
        # takes entry in selected line (temp_list) from line_file
        # and splits the entries from the whitespaces. Then appends them to the decrpyt_list
    while '' in list_decrypt:
        list_decrypt.remove('')

    if one_element_row is False:
        second_row = line_file[line_number+1]
        second_row = second_row.split(' ')
        while '' in second_row:
            second_row.remove('')
        # if any(chr.isdigit() for chr in second_row[0]) == False:
        #    pass
        # else:
        for word in second_row:
            list_decrypt.append(word)
    while '' in list_decrypt:
        list_decrypt.remove('')

    temp_dictionary[list_decrypt[0]] = list_decrypt[1:]
    return temp_dictionary


def run_theriak(theriak_path, database, temperature, pressure, whole_rock):
    """
    Function to run theriak with specified P-T condition and returns output as a list.
    it includes the path where theriak can be executed, writes the Therin file for the P-T
    conditions and uses a specific chemical composition (NannoOoze - Plank 2014) - date: february, 8 2021)

    Args:
        temperature (int): value for temperature condition
        pressure (int): value for pressure condition

    Returns:
        list: output (theriak in lines) from theriak as a list - ready to read and process
    """
    # Old way of calling theriak - now with input by user in init
    """
    if platform.system() == 'Windows':
        data_folder = Path(theriak_path)
        file_to_open = data_folder / "theriak.exe"
    else:
        data_folder = Path(theriak_path)
        file_to_open = data_folder / "theriak"
    """

    """
    main_folder = Path(__file__).parent.absolute()
    if platform.system() == 'Windows':
        data_folder = main_folder / "_theriak" / "WIN" / "Programs"
        file_to_open = data_folder / "theriak.exe"
    else:
        data_folder = main_folder / "_theriak" / "MAC" / "Programs"
        file_to_open = data_folder / "theriak"
    """
    """path_split = theriak_path.split("\\")
    print(path_split)
    pos = path_split.index("_theriak")
    path_split = '\\'.join(path_split[:pos])

    # define THERIN and XBIN location
    therin = Path(path_split) / 'THERIN'
    xbin = Path(path_split) / 'XBIN'"""

    """# stores the momentary P, T condition passed to Theriak for calculation
    with open(therin, 'w') as file_object:
        file_object.write(therin_condition)
    # opens THERIN and writes new P,T condition
    # with open('THERIN', 'a') as file_object:
        file_object.write("\n")
        file_object.write(whole_rock_write)
    with open(xbin, 'w') as file_object:
        file_object.write(database)
        file_object.write("\n")
        file_object.write("no")"""
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


    # reads each line in the output.txt and apends it to the 'theriak_in_lines'-list
    return theriak_output


def read_to_dataframe(index_entry, df_index_as_list, input_in_lines, skip=2, number=200):
    """
    A file read from lines to a list is transformed to a Dataframe

    Args:
        index_entry (int): index number in list that has been read from a txt file and got checked for a specific entry
        df_index_as_list (list): passing the header for the DataFrame (moles etc from theriak output)
        input_in_lines (list): file from txt file has been read to a list by its lines
        skip (int, optional): Sometimes there are lines to skip before the aimed entry is coming. Defaults to 2.
        number (int, optional):
            defines how long the loop is running for the entries
            - sometimes not necessary - could be done smarter. Defaults to 200.

    Returns:
        Dataframe: A dataframe with the defined headers and the values
            for each stable phase at the momentary P-T condition
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


def read_theriak(theriak_path, database, temperature, pressure, whole_rock):
    """
    Starts and read the theriak minimization at specific P-T and bulk rock to Dataframes that
    contain the physical values for each stable phase,
    elemental distribution, solid solutions and phases with water content.

    Args:
        temperature (int): temperature at given condition
        pressure (int): pressure at given condition
        whole_rock (list of floats): calculated bulk rock for the system

    Returns:
        Dictionary: A dictionary that contains all the created dataframes
    """
    # Print the calculation inforamtion
    # print(
    #     f"--------Calculation for P = {pressure} Bar and T = {temperature} °C-------"
    #     "With whole rock:"
    # )
    # print(whole_rock)
    # Call the connection to the console and execute theriak
    Data_dic = {}
    del Data_dic
    Data_dic = {}
    theriak_in_lines = run_theriak(theriak_path,
        database, temperature, pressure, whole_rock)

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
        print("=1= No gases and fluids in the system detected.")
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
        print("---ERROR:---No solid stable phases detected.")
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
        print("---ERROR:---Value Error with Volume and Desity index")

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
        print("---EXCEPTION--- no h2o content true")
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


def garnet_bulk_from_dataframe(frame, moles, volumePermole, volume):
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


class Therm_dyn_ther_looper:
    """
        Class for thermodynamic minimization for P-T-t steps and includes data
        reduction scheme. Data can be read after calling different functions.
        Order to call:
            1. initialize 'Therm_dyn_ther_looper'
            2. execute minimization with 'thermodynamic_looping_station'
            3. data reduction with 'merge_dataframe_dic' to incorporate dataframe into dictionaries
            4. Checking for volume changes and making variables easy to access via 'step_on_water'
    """

    def __init__(
            self, theriak_path, database, bulk_rock,
            temperature, pressure, df_var_dictionary,
            df_hydrous_data_dic, df_all_elements, num):
        """
        Initializes variables for thermodynamic minimization and data formation.
        -Bulk rock in normalized in mol to 1 kg of rock.
        -Temperature and pressure as a value from list.
        -Saving data in dictionaries and dataframe.

        Args:
            bulk_rock (list): is a list of values for each element as mol normalized to 1 kg of rock.
                        weight percent values are in the first looping round calculated with
                        'whole_rock_to_weight_normalizer' to mol per 1kg of rock.

            temperature (int): one temperature value in °C from list from P-T-t path
            pressure (type):
                one pressure value in bar from list from P-T-t path
            df_var_dictionary (dictionary): At start empty. Stores for each stable phase
            its physical properties (density, weight) a dataframe compiled in this dic
            df_hydrous_data_dic (dictionary): At start empty. Stores for each hydrous
            stable phase its physical properties (weight, H2O content) a dataframe compiled in this dic
            df_all_elements (dataframe): Frame that stores for each phase plus the
            total the respective element amount (O, C, H etc.)
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

    def thermodynamic_looping_station(self, marco=False):
        """
        Calls theriak and passes T, P and bulk rock. Resulting data is read and formatted.
        First check for 'water.fluid' is included.
        Prints messages and gives Volume, Density, Mol and Weight of present fluid as feedback.
        """

        # Flag: is set to True if water.fluid is stable in the system for further calculations
        content_h2o = True

        # Calls function to read the theriak output
        # (function also runs theriak and passes temperature, pressure and whole rock)
        # - returns theriak_data (which is a dictionary with Dataframes)

        theriak_data, g_sys, pot_frame = read_theriak(
            self.theriak_path, self.database, self.temperature, self.pressure, self.bulk_rock)

        # recalculation to delete for oversaturation of water
        if marco is False:
            if self.num < 1:
                if 'water.fluid' in theriak_data['df_elements_in_phases'].columns:
                    # print("333333 - First step free fluid detected - recalc initialized")
                    new_bulk = {}
                    water_H = theriak_data['df_elements_in_phases']['water.fluid']['H']
                    total = theriak_data['df_elements_in_phases']['total:']
                    total['H'] = total['H'] - water_H
                    for el in total.index:
                        if el == 'O' or el == 'E':
                            pass
                        else:
                            new_bulk[el] = total[el]
                    reset_bulk = check_redo_bulk(new_bulk)
                    self.bulk_rock = reset_bulk
                    theriak_data, g_sys, pot_frame = read_theriak(
                        self.theriak_path, self.database, self.temperature, self.pressure, self.bulk_rock)

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
            if 'water.fluid' in self.df_hydrous_data.columns:
                content_h2o = True
                # print("=3=== free water is present.")
            else:
                content_h2o = False
                print("3=== only hydrous solid phases, no free water.")
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
            self.new_fluid_Vol = self.df_phase_data.loc['volume[ccm]', 'water.fluid']
            self.new_fluid_dens = self.df_phase_data.loc['density[g/ccm]', 'water.fluid']
            self.new_fluid_N = self.df_phase_data.loc['N', 'water.fluid']
            self.new_fluid_weight = self.df_phase_data.loc['wt[g]',
                                                           'water.fluid']

        # print("\t Information on the fluid:")
        # print('\t Vol:{} Dens:{} N:{} Weight:{}'.format(self.new_fluid_Vol,
        #       self.new_fluid_dens, self.new_fluid_N, self.new_fluid_weight))

    def merge_dataframe_dic(self):
        """
        phase chache contains newly calculated phase data and is concenated
        with "df_var_dictionary" for each physical property ('N', 'volume', 'density' etc)"
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
        if 'water.fluid' in list(self.df_all_elements.columns):
            fluid_volumes = np.array(
                self.df_var_dictionary['df_volume[ccm]'].loc['water.fluid'])
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

            if 'water.fluid' in list(self.sys_vol_data.columns):
                # Sums up solid phase volume for each temperature (and pressure) step
                solid_volume = self.sys_vol_data.loc[temperature, :].sum(
                ) - self.sys_vol_data.loc[temperature, 'water.fluid']
                self.solid_vol_recalc.append(solid_volume)
            else:
                solid_volume = self.sys_vol_data.loc[temperature, :].sum()
                self.solid_vol_recalc.append(solid_volume)

    def mineral_fractionation(self, oxygen_data, name):
        """
        uses stable phases and checks for garnet
        we need a fractionation for garnet so we substract the element amount of garnet from the bulk
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


class Ext_method_master:
    """
    Module to calculate the factor important for the extraction of the water from the system
    """

    def __init__(
            self, pressure,
            fluid_volume_before, fluid_volume_new,
            solid_volume_before, solid_volume_new,
            save_factor, master_norm, phase_data,
            tensile_s, subduction_angle, rock_item_tag=0):
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
        # assumption and connected to main code - should be sent by input
        self.tensile_strength = tensile_s
        self.fracture = False
        self.shear_stress = 0
        self.frac_respo = 0
        self.angle = subduction_angle
        self.diff_stress = 0
        self.failure_dictionary = {}

    def couloumb_method(self, t_ref_solid, tensile=20):
        """
        Extraction of water/fluid follows the overstepping of a certain factor - calculated similar to Etheridge (2020)
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
        06.02.2023
        Assume
        1) differential stress (as inout from init file),
        2) tensile strength,
        3) fluid pressure from volume change
        Blanpied et al. 1992 state that low effective stresses can lead to slide failure
        or even extensional failure eventhough the system has low differential stresses.
        Extraction of water/fluid when the shear envelope or tensile strength are intersectedf
        - idea similar to Etheridge (2020) - no fluid factor approach
        """
        print("Coulomb test method 2 answer:")
        print('\tsolid_t0:{}\n\tsolid_t1:{}\n\tfluid_t0:{}\n\tfluid_t1:{}'.format(self.solid_t0,
              self.solid_t1, self.fluid_t0, self.fluid_t1))

        # test before executing module - phase data fluid volume should be equal the self.fluid_t1
        if self.phase_data['water.fluid']['volume[ccm]'] != self.fluid_t1:
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

    def mohr_cloulomb_griffith(self, shear_stress, friction, cohesion):
        """
        026.06.2023
        Estimates on
        1) differential stress (sig1 is lithostatic, shear resolves sig3),
        2) tensile strength (estimates from literature and experiments therein),
        3) fluid pressure from volume change (Duesterhoft 2019)
        Blanpied et al. 1992 state that low effective stresses can lead to slide failure
        or even extensional failure eventhough the system has low differential stresses.
        Extraction of water/fluid when the shear envelope or tensile strength are intersectedf
        - idea similar to Etheridge (2020) - no fluid factor approach
        """

        # NOTE testing with plot
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
        print('\tsolid_t0:{}\n\tsolid_t1:{}\n\tfluid_t0:{}\n\tfluid_t1:{}'.format(self.solid_t0,
              self.solid_t1, self.fluid_t0, self.fluid_t1))

        # test before executing module - phase data fluid volume should be equal the self.fluid_t1
        if self.phase_data['water.fluid']['volume[ccm]'] != self.fluid_t1:
            print("Inconsistency in new fluid volume")
            # keyboard.wait('esc')

        #################################
        # Mohr-Coulomb rock yield
        # linear equation of the critical line
        # 50 for y axis intercept: large depth with 200MPa < normal stress < 2000 MPa
        # slope of 0.6
        # LINK Mohr-Coulomb slope for failure
        cohesion = cohesion
        internal_friction = friction

        # REVIEW - static fix to 45°
        # 45 degree gives that diff stress is two time shear stress (and not more)
        self.angle = 27

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

        # #########################################
        # Normal stress of the system defined after Cox et al 2010
        # normal stress after Cox et al. 2010
        # normal stress = lithostatic when theta is 90°
        normal_stress = ((sig1+sig3)/2) - ((sig1-sig3)/2) * np.cos(2*self.angle*np.pi/180)

        # #########################################
        # Fluid pressure
        # get system conditions at present step and previous step
        vol_t0 = self.solid_t0 + self.fluid_t0
        vol_new = self.solid_t1 + self.fluid_t1
        # Fluid pressure calculation
        # Duesterhoft 2019 method
        hydro = normal_stress + normal_stress/vol_t0 * (vol_t0-(vol_new-vol_t0))-normal_stress

        # Mohr circle arguments
        r = self.diff_stress/2      # radius of the circle
        center = sig1-r             # centre of the cricel
        pos = center - hydro        # position shift of centre due to fluid pressure

        # Coulomb failure envelope
        a = -1
        b = 1/internal_friction
        c = -cohesion/internal_friction

        # Critical fluid pressures
        crit_fluid_pressure = cohesion/internal_friction + ((sig1+sig3)/2) - (
                (sig1-sig3)/2)*np.cos(2*self.angle*np.pi/180) - ((sig1-sig3)/2/internal_friction)*np.sin(2*self.angle*np.pi/180)
        pf_crit_griffith = (8*self.tensile_strength*(sig1+sig3)-((sig1-sig3)**2))/(16*self.tensile_strength)

        # ##########################################
        # Failure envelope test
        if self.diff_stress >= self.tensile_strength*5.66:
            print("Compressive shear failure test")
            print(f"Diff.-stress is {self.diff_stress} and > than T*5.66")

            # Test possible extensional fracturing
            output, minimum = checkCollision_linear(a, b, c, x=pos, y=0, radius=r)
            print(f"Maximum differential stress calculated as {minimum*2} MPa, but is {self.diff_stress}.")

            # Critical fluid pressure
            print(f"Difference of fluid pressure and critical fluid pressure is:\n{crit_fluid_pressure-hydro} (Pf_crit - Pf)")

            if output is True:
                print("---> Failure due to compressive shear.")
                self.frac_respo = 3
            else:
                print("---> No failure")

            # mcg_plot(cohesion, internal_friction, self.diff_stress, self.shear_stress, sig1, hydro)

        elif self.diff_stress < self.tensile_strength*5.66 and self.diff_stress > 4*self.tensile_strength:
            print("Extensional shear failure test")
            print(f"Diff.-stress is {self.diff_stress} and < than T*5.66 and > T*4")
            output, minimum = checkCollision_curve(pos=pos, diff=self.diff_stress, tensile=cohesion/2)
            print(f"Maximum differential stress calculated as {minimum*2} MPa, but is {self.diff_stress}.")
            print(f"Difference of fluid pressure and critical fluid pressure is:\n{pf_crit_griffith-hydro} (Pf_crit - Pf)")

            if output is True:
                print("---> Failure due to extensional shear.")
                self.frac_respo = 2
            else:
                print("---> No failure")

            # mcg_plot(cohesion, internal_friction, self.diff_stress, self.shear_stress, sig1, hydro)

        elif self.diff_stress <= 4*self.tensile_strength:
            print("Pure extensional failure test")
            print(f"Diff.-stress is {self.diff_stress} and < T*4")

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

        self.failure_dictionary = {
            "friction coeff": copy.deepcopy(internal_friction),
            "cohesion": (cohesion),
            "shear stress": shear_stress,
            "diff stress":(self.diff_stress),
            "radius": copy.deepcopy(r),
            "center": copy.deepcopy(center),
            "effective position": copy.deepcopy(pos),
            "critical fluid pressure normal": copy.deepcopy(crit_fluid_pressure),
            "critical fluid pressure griffith": copy.deepcopy(pf_crit_griffith),
            "fluid pressure": copy.deepcopy(hydro),
            "normal pressure": copy.deepcopy(normal_stress),
            "vol t0":copy.deepcopy(vol_t0),
            "solid t0":copy.deepcopy(self.solid_t0),
            "fluid t0":copy.deepcopy(self.fluid_t0),
            "vol t1":copy.deepcopy(vol_new),
            "solid t1":copy.deepcopy(self.solid_t1),
            "fluid t1":copy.deepcopy(self.fluid_t1),
            "sigma 1":copy.deepcopy(sig1),
            "sigma 3":copy.deepcopy(sig3),
            "angle": copy.deepcopy(self.angle),
            "fracture key": copy.deepcopy(self.frac_respo)
            }

        if self.frac_respo > 0:
            fracturing = True
        elif self.frac_respo == 0:
            fracturing = False
        else:
            print("The fracturing response failed. Thorsten has to work more here.")


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

    def __init__(self, phase_data, ext_data, temperature, new_fluid_V, sys_H, element_frame, st_fluid_post):
        self.phase_data_fluid = phase_data
        self.ext_data = ext_data
        self.temp = temperature
        self.new_fluid_V = new_fluid_V # REVIEW can be changed = self.phase_data_fluid['volume[ccm]']
        self.sys_H = sys_H # self.element_frame.loc['H', 'total:']
        self.element_frame = element_frame.copy()
        self.st_fluid_post = st_fluid_post

    def hydrogen_ext_all(self):
        """
        recalculates hydrogen content due to fluid extraction
        """
        new_extraction = self.phase_data_fluid
        self.ext_data = pd.concat([self.ext_data, new_extraction], axis=1)
        self.ext_data = self.ext_data.rename(
            columns={"water.fluid": self.temp})

        # diminish double total: in dataframe
        if self.element_frame.columns[-1] == self.element_frame.columns[-2]:
            self.element_frame = self.element_frame.iloc[:, :-1]

        # settings some important values from element data frame
        extracted_fluid_vol = self.phase_data_fluid['volume[ccm]']
        fluid_H = self.element_frame.loc['H', 'water.fluid']
        sys_H = self.element_frame.loc['H', 'total:']
        fluid_O = self.element_frame.loc['O', 'water.fluid']
        sys_O = self.element_frame.loc['O','total:']

        # remaining H and O in bulk
        total_hydrogen = sys_H - fluid_H
        total_oxygen = sys_O - fluid_O

        # IMPORTANT step - substracting H and O of fluid from total system - total system will be new bulk
        self.element_frame.loc['H', 'total:'] = total_hydrogen
        self.element_frame.loc['O', 'total:'] = total_oxygen

        # set fluid O and H to zero after this step
        self.element_frame.loc['O', 'water.fluid'] = 0.0
        self.element_frame.loc['H', 'water.fluid'] = 0.0
        self.st_fluid_post[-1] = self.st_fluid_post[-1] - extracted_fluid_vol
        # print("=====================================")
        print("=======yummy yummy yummy, I drunk all your water========")
        # print("=====================================")

    def hydrogen_partial_ext(self):
        pass


class System_status:
    """
    Compiling the system conditions for the P-T-t path. Ready for plotting and investigations.
    """

    def __init__(self, df_var_dictionary, df_h2o_content_dic, df_element_total, st_elements):

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
            [df_ABC[phase][1], df_ABC[phase][2], df_ABC[phase][3]])

        for phase in list(flatten(database_names)):
            a = df_ABC[phase][1] - norm_ABC[0]
            b = df_ABC[phase][2] - norm_ABC[1]
            c = df_ABC[phase][3] - norm_ABC[2]
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
        oxygen_moles = copy.deepcopy(np.array(self.element_oxy))
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
        moles_frac = list(np.nan_to_num(moles_frac, 0))
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
        print("Oxygen fractionation calculation done.")


class Garnet_recalc():

    def __init__(self, theriak, dataclass, temperature, pressure):
        self.mineral = dataclass
        self.recalc_volume = 0
        self.recalc_weight = 0
        self.temperature = temperature
        self.pressure = pressure
        self.theriak_path = theriak

    def recalculation_of_garnets(self):
        for garnet in self.mineral:
            vals = garnet.elements[0]
            index = garnet.elements[1]
            relements = pd.DataFrame(vals, index=index)

            # create the bulk from element entry - normalized to 1 mol for theriak
            # forward - volume needs to be back-converted to the moles of the shell (garnet.moles)
            bulk = garnet_bulk_from_dataframe(relements, garnet.moles, garnet.volPmole, garnet.volume)

            # FIXME static database tc55 for garnet recalculation
            db = "tc55.txt    GARNET"
            Data_dic, g_sys, pot_frame = read_theriak(
                self.theriak_path,
                database=db, temperature=self.temperature,
                pressure=self.pressure, whole_rock=bulk)
            grt = Data_dic['df_Vol_Dens'].columns[0]
            phase_list = list(Data_dic['df_Vol_Dens'].columns)

            # check for couble garnet
            double_grt_check = []
            selected_garnet_name = False
            for item in phase_list:
                if 'GARNET' in item:
                    selected_garnet_name = item
                    double_grt_check.append(item)

            if selected_garnet_name is False:
                print("WARNING: Garnet increment no longer stable. Press ESC to continue.")
                print(f"Bulk is: {bulk}")
                # keyboard.wait('esc')
                v = 0
            else:
                # backward - Volume of the garnet shell to be added (back-normlaized from 1 mole to x-mol (garnet.moles))
                rec_mole = garnet.volume/garnet.volPmole
                v = Data_dic['df_Vol_Dens'][grt]['volume/mol']*rec_mole
                g = Data_dic['df_Vol_Dens'][grt]['wt/mol']*rec_mole

            self.recalc_volume += v
            self.recalc_weight += g
