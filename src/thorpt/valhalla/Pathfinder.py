
"""
Written by
Thorsten Markmann
thorsten.markmann@unibe.ch
status: 16.07.2024
"""

import pandas as pd
import numpy as np
from pathlib import Path
import time
import scipy.special

import platform
import json
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from tkinter import Tk, filedialog, messagebox, simpledialog
from scipy.interpolate import splev, splrep
from scipy.interpolate import Rbf, InterpolatedUnivariateSpline


def getReferenceLength(index):
    '''
    Get the reference length in the requested direction

    Args:
        index (int): Index value representing the direction. Use 0 for x-direction or 1 for y-direction.

    Returns:
        tuple: A tuple containing the scaling factor, the coordinate value of the starting point, and the origin value.

    Raises:
        None

    Usage:
        factor, start_coord, origin_value = getReferenceLength(index)
    '''

    # define a 'direction' string
    direction = 'x' if index == 0 else 'y'
    # get the reference length
    length_selected = False
    while not length_selected:
        messagebox.showinfo(
            "Select reference length",
            "Use the mouse to select the reference length in {:s} direction.".format(direction) +
            "Click the start and the end of the reference length."
        )
        coord = plt.ginput(
            2,
            timeout=0,
            show_clicks=True
        )  # capture only two points
        # ask for a valid length
        valid_length = False
        while not valid_length:
            reflength = simpledialog.askfloat(
                "Enter reference length",
                "Enter the reference length in {:s} direction".format(direction))
            if isinstance(reflength, float):
                valid_length = True
            else:
                messagebox.showerror("Error", "Please provide a valid length.")
        # calculate scaling factor
        deltaref = coord[1][index]-coord[0][index]
        factor = reflength/deltaref
        # ask for origin values of plot
        valid_origin_value = False
        while not valid_origin_value:
            origin_value = simpledialog.askfloat(
                "Enter origin value",
                "Enter origin value {:s} direction".format(direction))
            if isinstance(origin_value, float):
                valid_origin_value = True
            else:
                messagebox.showerror(
                    "Error", "Please provide a valid origin value.")
        length_selected = messagebox.askyesno(
            "Length confirmation",
            "You selected {:4.0f} pixels in {:s} direction"
            "corresponding to {:4.4f} units. Is this correct?".format(
                deltaref, direction, reflength)
        )
    origin_value = float(origin_value)
    return factor, coord[0][index], origin_value

def read_temperature_pressure_txt():
    """
    Reads temperature and pressure information from a txt file.

    Returns:
        temperatures (numpy.ndarray): Array of temperature values.
        pressures (numpy.ndarray): Array of pressure values.
    """

    # GUI to select the txt file of P-T information
    filein = filedialog.askopenfilename(
            title="Select a digitized path file",
            filetypes=[("other", ".txt")]
            )

    # reading information to lines
    with open(filein) as f:
        lines = f.readlines()

    # read temperature from first entry and split the string to the separate values
    temperatures = lines[0]
    temperatures = temperatures.split()
    for i, item in enumerate(temperatures):
        temperatures[i] = np.float32(item)

    # read pressure from first entry and split the string to the separate values
    pressures = lines[1]
    pressures = pressures.split()
    for i, item in enumerate(pressures):
        pressures[i] = np.float32(item)

    # covert to array
    temperatures = np.array(temperatures)
    pressures = np.array(pressures)

    return temperatures, pressures

def layered_model_PTpatch(temperatures, pressures, layers, temperature_increase_to_bottom=100):
    """
    Creates a layered model for a P-T patch.

    Args:
        temperatures (list): The temperature array from the digitized P-T path.
        pressures (list): The pressure array from the digitized P-T path.
        layers_array (list): The array of layer thicknesses.
        temperature_increase_to_bottom (int): The temperature increase to the bottom.

    Returns:
        temperature_matrix (list): The temperature matrix.
        pressure_matrix (list): The pressure matrix.
    """

    layers_array = np.zeros(len(layers))
    # convert each entry in layers_array to float
    for i, entry in enumerate(layers):
        # print(entry)
        layers_array[i] = float(entry[0])
    slab_thickness = np.sum(layers_array)
    layers_array = layers_array
    layers_array = np.append(layers_array, 0)[::-1]
    positional_layer = np.cumsum(layers_array)[0:-1][::-1]

    # get temperature matrix, iterating slice-wise with increasing depth
    temperature_matrix = []
    for temperature_top in temperatures:
        temperature_bottom = temperature_top + temperature_increase_to_bottom
        temperature_array = temperature_top + (temperature_bottom - temperature_top) * \
                    scipy.special.erf(positional_layer/(slab_thickness/2))
        temperature_matrix.append(temperature_array)

    # get pressure matrix, iterating slice-wise with increasing depth
    pressure_matrix = []
    density = 3300
    for pressure_top in pressures:
        overburden = np.cumsum(layers_array)[:-1]
        pressure_array = pressure_top + density * 9.81 * overburden * 1e-5
        pressure_matrix.append(pressure_array)

    return temperature_matrix, pressure_matrix

def crust2layer_model(pressure_array, time, speed, angle, dt=10000):

    # read layermodel.txt
    lines = []
    try:
        with open('layermodel.txt') as f:
            lines = f.readlines()

            for item in lines:
                if 'layers' in item:
                    # split item by :
                    layers = item.split(':')
                    layers = layers[1].split(',')
                    number_of_layers = len(layers)
                if 'rho' in item:
                    rho = item.split(':')
                    rho = rho[1].split(',')
                    # add number to rho_list for each number of layers
                    density_list = np.zeros(number_of_layers)
                    for i in range(number_of_layers):
                        density_list[i] = float(rho[i])
                if 'thickness' in item:
                    thickness = item.split(':')
                    thickness = thickness[1].split(',')
                    layer_thickness_list = []
                    for val in thickness:
                        if val == 'increasing':
                            pass
                        else:
                            layer_thickness_list.append(float(val))

    except FileNotFoundError:
        print("Error: 'layermodel.txt' not found.")
        return

    c_p_zero = 0
    crust_depth = 0
    depth = []
    c_p_list = []
    time = []

    # calculated pressure with the layer model
    for i in range(len(layer_thickness_list)):
        c_p_zero += density_list[i] * layer_thickness_list[i] * 9.81 / 10**5

    c_p = c_p_zero
    d_step = speed * dt * abs(np.tan(angle/180*np.pi))
    if c_p_zero < pressure_array[0]:
        start_depth = pressure_array[0] * 10**5 / (density_list[0] * 9.81)
        depth.append(start_depth + sum(layer_thickness_list))
        start_time =  depth[-1] / speed / abs(np.sin(angle/180*np.pi))
        time.append(start_time)
        c_p = pressure_array[0]
        c_p_list.append(c_p)
        
    # c_p = c_p_zero + (density_list[-1]*(crust_depth) * 9.81 / 10**5)
    
    for c_p in pressure_array[1:]:
        calc_depth = c_p * 10**5 / (density_list[0] * 9.81)
        depth.append(calc_depth + sum(layer_thickness_list))
        calc_time =  depth[-1] / speed / abs(np.sin(angle/180*np.pi))
        time.append(calc_time)
        c_p_list.append(c_p)
        # print(calc_depth)


    if len(depth) == 0:
        print("Error: No depth calculated.")
    else:
        print('End-depth is: {} km'.format(depth[-1]))

    return c_p_list, time, depth

class Pathfinder_Theoule:
    """
    A class that represents a pathfinder for Theoule.

    Attributes:
    - temperatures (list): The temperature array from the digitized P-T path.
    - pressures (list): The pressure array from the digitized P-T path.
    - path_increment (list): The path increment values for pressure and temperature.
    - sub_angle (float): The subduction angle in degrees.
    - plate_speed (float): The plate speed in m/year.
    - dt (int): The time increment in years.
    - rho (list): The density values for the crust and mantle.
    - time (list): The time array for the path.
    - depth (list): The depth array for the path.
    - lower_t_bound (float): The lower temperature bound for filtering.

    Methods:
    - prograde(): Performs the prograde calculation for the path.
    - loop(): Performs the loop calculation for the path.
    """
    def __init__(self, temperatures, pressures, path_increment, sub_angle=False, plate_speed=False, dt=10000):
        """
        Initializes a Pathfinder_Theoule object.

        Parameters:
        - temperatures (list): The temperature array from the digitized P-T path.
        - pressures (list): The pressure array from the digitized P-T path.
        - path_increment (list): The path increment values for pressure and temperature.
        - sub_angle (float): The subduction angle in degrees. (default: False)
        - plate_speed (float): The plate speed in m/year. (default: False)
        - dt (int): The time increment in years. (default: 10000)
        """

        if plate_speed is False:
            self.speed = float(input("Give me a speed in m/year:\n"))
        else:
            self.speed = float(plate_speed)
        if sub_angle is False:
            self.angle = float(input("Give me a subduction angle in degree:\n"))
        else:
            self.angle = float(sub_angle)
        self.temp = temperatures
        self.pressure = pressures
        self.dt = dt
        self.rho = [2800, 3300]
        self.time = [0]
        self.depth = []
        self.p_increment = np.float64(path_increment[0])
        self.t_increment = np.float64(path_increment[1])
        self.lower_t_bound = np.float64(path_increment[2])

    def prograde(self):
        """
        Perform prograde calculation to determine the depth and temperature profile
        of a planetary body based on pressure and temperature data.

        Returns:
            None
        """
        # spl = self.prepare_spline()
        c_p_list = self.construct_layer_model()
        self.depth = c_p_list[2]
        self.depth_store = c_p_list[2]
        self.time_store = c_p_list[1]
        new_data = self.fit_model_to_path(c_p_list)
        new_data = self.select_steps(new_data)
        self.filter_steps(self.temp, [self.pressure, self.depth, self.time])

    def prepare_spline(self):
        """
        Prepare spline from input P-T data.

        Returns:
            spl (tuple): Spline representation of the input P-T data.
        """

        # when self.pressure and self.temp are len of 2, do a linear interpolation with 10 nodes
        if len(self.pressure) == 2:
            spl = (np.linspace(self.pressure[0], self.pressure[1], 10), np.linspace(self.temp[0], self.temp[1], 10), 1)
        else:
            spl = splrep(self.pressure, self.temp)

        return spl

    def construct_layer_model(self):
        """
        Construct layer model for pressure, depth, and time.

        Returns:
            c_p_list (list): List of calculated pressures.
        """
        return crust2layer_model(self.pressure, self.time, self.speed, self.angle, self.dt)

    def fit_model_to_path(self, c_p_list):
        """
        Fit the model to the P-T path.

        Args:
            c_p_list (list): List of calculated pressures.
            spl (tuple): Spline representation of the input P-T data.

        Returns:
            yinterp (numpy.ndarray): Interpolated temperature values.
        """
        from scipy.interpolate import splprep, splev
        data = [
            self.temp[:len(c_p_list[0])], 
            np.array(c_p_list[0]), 
            np.array(c_p_list[2]), 
            np.array(c_p_list[1])
            ]
        
        # Create a spline representation of the parametric curve
        tck, u = splprep(data, s=0)

        # Define new parameter values for evaluation
        new_u = np.linspace(0, 1, 1000)

        # Evaluate the spline at the new parameter values
        new_data = splev(new_u, tck)

        return new_data 

    def filter_steps(self, yinterp, c_p_list):
        """
        Filter the steps based on pressure and temperature increments.

        Args:
            yinterp (numpy.ndarray): Interpolated temperature values.
            c_p_list (list): List of calculated pressures.

        Returns:
            yinterp (numpy.ndarray): Filtered temperature values.
            c_p_list (numpy.ndarray): Filtered pressure values.
        """
        new_x = [yinterp[0]]
        new_y = [c_p_list[0][0]]
        if len(self.depth) == 0:
            self.depth = np.ones(len(c_p_list[0]))
        new_d = [self.depth[0]]
        if len(self.time) == 0 or len(self.time) == 1:
            self.time = np.ones(len(c_p_list[0]))
        new_t = [self.time[0]]

        
        for i, val in enumerate(c_p_list[0]):
            step_p = val - new_y[-1]
            step_t = yinterp[i] - new_x[-1]
            if step_p >= self.p_increment:
                if step_t >= self.t_increment:
                    new_y.append(val)
                    new_x.append(yinterp[i])
                    new_d.append(c_p_list[1][i])
                    new_t.append(c_p_list[2][i])
            elif step_t >= self.t_increment:
                if step_p >= self.p_increment:
                    new_y.append(val)
                    new_x.append(yinterp[i])
                    new_d.append(c_p_list[1][i])
                    new_t.append(c_p_list[2][i])
        yinterp = np.array(new_x)
        c_p_list = np.array(new_y)

        self.temp = yinterp
        self.pressure = c_p_list
        self.time = new_t
        self.depth = new_d


    def select_steps(self, data):
        """
        Select only steps with temperature >= lower temperature bound.

        Args:
            yinterp (numpy.ndarray): Filtered temperature values.
            c_p_list (numpy.ndarray): Filtered pressure values.

        Returns:
            None
        """

        frame = pd.DataFrame(data)
        cut_T = self.lower_t_bound
        yinterp = np.array(frame.iloc[0][frame.iloc[0] >= cut_T])
        c_p_list = np.array(frame.iloc[1][frame.iloc[0] >= cut_T])
        self.time = np.array(frame.iloc[3][frame.iloc[0] >= cut_T])
        self.depth = np.array(frame.iloc[2][frame.iloc[0] >= cut_T])
        self.temp = yinterp
        self.pressure = c_p_list


    def loop(self):
        """
        Performs the burial and exhumation process to generate a P-T path.
        """
        depth, crust_d, c_p_list, c_p, d_step = self._initialize_loop_variables()
        self._perform_burial(depth, crust_d, c_p_list, c_p, d_step)
        self._perform_exhumation(depth, crust_d, c_p_list, c_p, d_step)
        self._filter_and_update_results(c_p_list)

    def _initialize_loop_variables(self):
        depth = 1000  # in meter
        crust_depth = 1000  # in meter
        c_p_list = []
        c_p = (self.rho[0]*crust_depth + self.rho[1] * (depth-crust_depth)) * 9.81 / 10**5  # in Bar
        d_step = self.speed/100 * self.dt * abs(np.sin(self.angle/180*np.pi))  # in meter
        print("Start resampling into \x1B[3mP-T-t\x1B[0m path")
        return depth, crust_depth, c_p_list, c_p, d_step

    def _perform_burial(self, depth, crust_depth, c_p_list, c_p, d_step):
        while c_p <= max(self.pressure):
            c_p = (self.rho[0]*crust_depth + self.rho[1] * (depth-1)) * 9.81 / 10**5  # in Bar
            if c_p < self.pressure[0]:
                depth += d_step
            else:
                depth += d_step
                self.time.append(self.time[-1]+self.dt)
                self.depth.append(depth)
                c_p_list.append(c_p)
        print("Burial finished...")

    def _perform_exhumation(self, depth, crust_depth, c_p_list, c_p, d_step):
        if c_p > self.pressure[-1]:
            d_step = self.speed/100*5 * self.dt * abs(np.sin(self.angle/180*np.pi))  # in meter
            while c_p > self.pressure[-1]:
                c_p = (self.rho[0]*crust_depth + self.rho[1] * (depth-1)) * 9.81 / 10**5  # in Bar
                if c_p < self.pressure[-1]:
                    depth -= d_step
                else:
                    depth -= d_step
                    self.time.append(self.time[-1]+self.dt)
                    self.depth.append(depth)
                    c_p_list.append(c_p)
            print("Exhumation finished...")
        print('Final depth is: {} m'.format(depth))

    def _filter_and_update_results(self, c_p_list):
        self.time = self.time[1:]
        yinterp = self.temp
        c_p_list = self.pressure
        new_x, new_y, new_d, new_t = self._filter_steps(yinterp, c_p_list)
        yinterp = np.array(new_x)
        c_p_list = np.array(new_y)
        self.time = new_t
        self.depth = new_d
        self._select_steps_with_min_temp(yinterp, c_p_list)

    def _select_steps_with_min_temp(self, yinterp, c_p_list):
        frame = pd.DataFrame([yinterp, c_p_list, self.time, self.depth])
        cut_T = self.lower_t_bound
        yinterp = np.array(frame.iloc[0][frame.iloc[0] >= cut_T])
        c_p_list = np.array(frame.iloc[1][frame.iloc[0] >= cut_T])
        self.time = np.array(frame.iloc[2][frame.iloc[0] >= cut_T])
        self.depth = np.array(frame.iloc[3][frame.iloc[0] >= cut_T])
        self.temp = yinterp
        self.pressure = c_p_list
        self.rock_rho = [2800, 3300]  # kg/m3
        self.X_val = [0]
        self.Y_val = []

    def calc_time_model(self,
                        timestep=1000, t_end=33e6, start_depth=20,
                        end_depth=80_000, t_start=144e6, rate=1.5, angle=15,
                        start_T=350, end_T=600, dT=10, start_p=5000, end_p=20000, mode='time'):
        """
        iterating over time and creating P-T-t path or generating a line path
        """
        if mode == 'time':
            self.timestep = timestep  # in years
            self.t_end = t_end
            self.t_start = [t_start]  # default is 144_000_000 years (144 Ma)

            start_depth = start_depth
            crust_thickness = 1000  # default layer on top
            self.rate = rate/100  # m/year
            self.T = [start_depth*self.geotherm+200]  # Temperature in °C
            self.P = [self.rock_rho[1] * start_depth * 9.81]  # N/m2
            self.P = [(self.rock_rho[0]*crust_thickness + self.rock_rho[1]
                       * (start_depth-crust_thickness)) * 9.81]  # N/m2

            self.Y_val = [start_depth]
            self.end_depth = end_depth

            # self.angle = 15  # degree
            self.angle = angle

            nt = (self.t_start[-1] - self.t_end) / self.timestep
            print("Start path calculation. Please wait...")
            while self.t_start[-1] > self.t_end:

                # print(f"The time is: {self.t_start[-1]/1e6}")
                Y1 = self.Y_val[-1]
                x_step = self.rate * self.timestep
                y_step = self.rate * self.timestep * \
                    abs(np.sin(self.angle/180*np.pi))
                x = self.X_val[-1] + x_step
                y = self.Y_val[-1] + y_step

                temp_step = self.geotherm * (y-Y1)
                press_step = self.rock_rho[1] * (y-Y1) * 9.81
                T = self.T[-1] + temp_step
                P = self.P[-1] + press_step

                self.X_val.append(x)
                self.Y_val.append(y)
                self.T.append(T)
                self.P.append(P)

                self.t_start.append(self.t_start[-1] - self.timestep)

                if self.Y_val[-1] > self.end_depth:
                    print("Final depth is reached abort mission")
                    break

            print(f"Depth is {y} Meter")
            print(f"Pressure is {P/1e9} GPa")
            print(f"Temperature is {T} °C")
        elif mode == 'line':
            self.T = np.arange(start_T, end_T, dT)
            self.P = np.linspace(start_p, end_p, len(self.T))


class Pub_pathfinder:
    """
    Modul to read published P-T-paths from Penniston-Dorland (2015)
    """

    def __init__(self, name_tag="Colombia_Ecuador.txt"):

        self.file_name = name_tag
        self.temperatures = 0
        self.pressures = 0

    def published_path(self):
        """
        - reading file for P-T-path from Penniston-Dorland (2015) publication
        - files are in paper_path/D80/
        - data is passed to frame and only temperatures between 290 and 700 °C are regarded
        - P and T steps are extracted for theriak minimization
        """

        # selecting path and file depending on OS
        main_folder = Path(__file__).parent.absolute()
        data_folder = main_folder / "paper_paths/D80/"

        if platform.system() == 'Windows':
            file_to_open = data_folder / self.file_name
        else:
            file_to_open = data_folder / self.file_name[:-4]

        # reading file and save to DataFrame
        frame = pd.read_csv(file_to_open, sep=" ", header=None)
        # adjusting and cutting data
        frame = frame[frame[0] == 0]
        frame = frame[frame[3] < 700]
        frame = frame[frame[3] > 290]
        # storing important P-T values
        self.temperatures = frame[3]
        self.pressures = (frame[2]+frame[1])/100*10000


class PTPathDigitizer:
    """
    Modul to digitize a P-T path from plots and extract P and T values
    or use existing txt file saves from previous paths
    """

    def __init__(self):
        self.temperatures = 0
        self.pressures = 0

    def run(self):
        '''
        Main function of curve digitizer
        '''

        # open the dialog box
        # first hide the root window
        root = Tk()
        root.withdraw()
        # open the dialog
        filein = filedialog.askopenfilename(
            title="Select image to digitize",
            filetypes=(
                ("jpeg files", "*.jpg"),
                ("png files", "*.png"))
        )
        root.update()
        if len(filein) == 0:
            # nothing selected, return
            return

        # show the image
        img = mpimg.imread(filein)
        _, ax = plt.subplots()

        ax.imshow(img)
        ax.axis('off')  # clear x-axis and y-axis

        # get reference length in x direction
        xfactor, xorigin, base_val_x = getReferenceLength(0)

        # get the reference length in y direction
        yfactor, yorigin, base_val_y = getReferenceLength(1)

        # digitize curves until stoped by the user
        reply = True
        while reply:

            messagebox.showinfo(
                "Digitize curve",
                "Please digitize the curve. The first point is the origin." +
                "Left click: select point; Right click: undo; Middle click: finish"
            )

            # get the curve points
            x = plt.ginput(
                -1,
                timeout=0,
                show_clicks=True
            )
            x = np.array(x)

            ax.plot(x[:, 0], x[:, 1], 'g', 'linewidth', 1.5)
            plt.draw()

            # convert the curve points from pixels to coordinates
            coords = np.array([x[:, 0], x[:, 1]])
            coords[0] = base_val_x + (x[:, 0] - xorigin) * xfactor
            coords[1] = base_val_y + (x[:, 1] - yorigin) * yfactor

            # write the data to a file
            # first get the filename
            validFile = False

            while not validFile:
                fileout = filedialog.asksaveasfilename(
                    title="Select file to save the data",
                    filetypes=[("Simple text files (.txt)", "*.txt")],
                    defaultextension='txt'
                )
                if len(fileout) == 0:
                    # nothing selected, pop up message to retry
                    messagebox.showinfo(
                        "Filename error", "Please select a filename to save the data.")
                else:
                    validFile = True

            # write the data to file
            np.savetxt(fileout, coords)
            # with open(fileout, 'w', encoding='utf-8') as f:
            #     f.write(np.array2string(coords[0], precision=4, separator='\t'))
            #     f.write("\n")
            #     f.write(np.array2string(coords[1], precision=4, separator='\t'))
            # coords.tofile(fileout, sep='\t', format='%s')

            self.temperatures = coords[0]
            self.pressures = coords[1]

            reply = messagebox.askyesno(
                "Finished?",
                "Digitize another curve?"
            )

        # clear the figure
        plt.clf()
        plt.close()

    def stored_digitization(self):
        """
        Main method to handle the stored digitization process.
        """
        filein = self._get_file_path()
        lines = self._read_file(filein)
        self.temperatures, self.pressures = self._parse_lines(lines)

    def _get_file_path(self):
        """
        Get the file path for the digitized P-T path file.
        """
        return filedialog.askopenfilename(
            title="Select a digitized path file",
            filetypes=[("other", ".txt")]
        )

    def _read_file(self, filein):
        """
        Read the contents of the file.
        """
        with open(filein) as f:
            return f.readlines()

    def _parse_lines(self, lines):
        """
        Parse the lines from the file to extract temperatures and pressures.
        """
        temperatures = lines[0].split()
        pressures = lines[1].split()
        temperatures = [np.float32(item) for item in temperatures]
        pressures = [np.float32(item) for item in pressures]
        return np.array(temperatures), np.array(pressures)



        """
        Executes the second modified digitization process for a P-T path.

        Args:
            path_arguments (bool or list, optional): The path arguments (default is False).
            path_increment (bool, optional): The path increment value (default is False).
        """
        # Code implementation...


class core_Pathfinder:
    """
    A class to handle the Pathfinder module for digitizing or using stored P-T paths.
    Methods:
        __init__(temp=0, pressure=0, dt=1000):
            Initializes a Pathfinder object with default temperature, pressure, and time step values.
        execute_digi_prograde():
        execute_digi():
        execute_digi_mod():
        digitilization_module(path_arguments=False, path_increment=False):
        loop_digi(path_arguments=False, path_increment=False):
        gridding(path_arguments, path_increment):
    Private Methods:
        _select_digitization_method(path_arguments):
            Selects the digitization method based on user input or provided arguments.
        _store_image_pt_path(path_arguments):
            Stores the image P-T path in arrays.
        _convert_pressure_units(pressures, path_arguments):
            Converts the pressure units based on user input or provided arguments.
        _create_f_path_dataframe(temperatures, pressures):
            Creates a DataFrame for the P-T path.
        _test_prograde_peak_retrograde(f_path):
            Tests for prograde, peak, and retrograde paths in the P-T path.
        _apply_interpolation_loop(f_path, pressures):
            Applies interpolation to the P-T path.
        _create_pathfinder(temperatures, pressures, path_arguments, path_increment):
            Creates a Pathfinder_Theoule object and updates the instance variables.
        _choose_digitization_method(path_arguments):
            Chooses the digitization method based on user input or provided arguments.
        _extract_prograde_retrograde_paths(f_path):
            Extracts prograde and retrograde paths from the P-T path.
        _interpolate_paths(t_pro, p_pro, t_ret, p_ret, pressures):
            Interpolates the prograde and retrograde paths.
        _create_pt_path(temperatures, pressures, path_arguments, path_increment):
            Creates a P-T path and updates the instance variables.

        # External call to use Pathfinder module
        # Decide by input if you want to digitise a P-T path from image
        # or
        # use a stored P-T path from txt file
        """

    def __init__(self, temp=0, pressure=0, dt=1000):
        """
        Initializes a Pathfinder object.

        Args:
            temp (int): The temperature value (default is 0).
            pressure (int): The pressure value (default is 0).
            dt (int): The time step value (default is 1000).
        """
        self.temp = temp
        self.pressure = pressure
        self.time_var = 0
        self.depth = 0
        self.dt = dt
        # Calling digitizing module
        self.new_or_read = PTPathDigitizer()


    def digitilization_module(self, path_arguments=False, path_increment=False):
        """
        Executes the second modified digitization process for a P-T path.

        Args:
            path_arguments (bool or list, optional): The path arguments (default is False).
            path_increment (bool, optional): The path increment value (default is False).
        """
        # TODO writing argument for multiple peak and retrograde paths - complex loop
        # idea get pieces of prograde and retrograde path and the number of slices then iterate
        # is the BACK-UP that 'vho path' is using

        self._select_digitization_method(path_arguments)
        temperatures, pressures = self._store_image_pt_path(path_arguments)
        pressures = self._convert_pressure_units(pressures, path_arguments)
        f_path = self._create_f_path_dataframe(temperatures, pressures)

        # test in P-T path for prograde, peak and retrograde paths (temperature is f_path[0] and pressure is f_path[1])
        t_peak, p_peak, t_ret, p_ret = self._test_prograde_peak_retrograde(f_path)

        if len(t_peak) > 0 and len(t_ret) > 0:
            # interpolate the prograde and retrograde paths
            temperatures, pressures = self._apply_interpolation_loop(f_path, pressures)

            if path_arguments is False:
                nasa = Pathfinder_Theoule(
                    temperatures, pressures,
                    path_increment=[500, 15, 350],
                    dt=self.dt
                )
            else:
                nasa = Pathfinder_Theoule(
                    temperatures, pressures,
                    plate_speed=path_arguments[3],
                    sub_angle=path_arguments[4],
                    dt=self.dt,
                    path_increment=path_increment
                )

            nasa.loop()

            self.temp = nasa.temp
            self.pressure = nasa.pressure
            self.time_var = nasa.time
            self.depth = nasa.depth
            self.sub_angle = nasa.angle
            self.plate_v = nasa.speed
            loop = True

        elif len(t_ret) > 0 and len(t_peak) == 0:
            # Using a retrograde path

            if np.shape(f_path)[0] == 2:
                temperatures = np.linspace(min(temperatures), max(temperatures), 100)[::-1]
                pressures = np.linspace(min(pressures), max(pressures), 100)[::-1]
            else:
                temperatures, pressures = self._apply_interpolation_retrograde(f_path, pressures)

            if path_arguments is False:
                nasa = Pathfinder_Theoule(
                    temperatures, pressures,
                    path_increment=[500, 15, 350],
                    dt=self.dt
                )
            else:
                nasa = Pathfinder_Theoule(
                    temperatures, pressures,
                    plate_speed=path_arguments[3],
                    sub_angle=path_arguments[4],
                    dt=self.dt,
                    path_increment=path_increment
                )

            nasa.loop()

            self.temp = nasa.temp
            self.pressure = nasa.pressure
            self.time_var = nasa.time
            self.depth = nasa.depth
            self.sub_angle = nasa.angle
            self.plate_v = nasa.speed
            loop = True

        else:
            # Using a prograde path

            if np.shape(f_path)[0] == 2:
                temperatures = np.linspace(min(temperatures), max(temperatures), 100)
                pressures = np.linspace(min(pressures), max(pressures), 100)

            if path_arguments is False:
                nasa = Pathfinder_Theoule(
                    temperatures, pressures,
                    path_increment=[500, 15, 350],
                    dt=self.dt
                )
            else:
                nasa = Pathfinder_Theoule(
                    temperatures, pressures,
                    plate_speed=path_arguments[3],
                    sub_angle=path_arguments[4],
                    dt=self.dt,
                    path_increment=path_increment
                )

            nasa.prograde()

            self.temp = nasa.temp
            self.pressure = nasa.pressure
            self.time_var = nasa.time
            self.depth = nasa.depth
            self.sub_angle = nasa.angle
            self.plate_v = nasa.speed

    def _select_digitization_method(self, path_arguments):
        answers = ["new", "stored"]
        if path_arguments is False:
            for val in answers:
                print(val)
            answer = input("Pathfinder function - new or stored path? Select answer\n")
        else:
            answer = path_arguments[1]

        if answer == answers[0]:
            self.new_or_read.run()
        elif answer == answers[1]:
            self.new_or_read.stored_digitization()
        else:
            print("Unexpected end - no P-T file input")
            time.sleep(10)
            exit()

    def _store_image_pt_path(self, path_arguments):
        temperatures = self.new_or_read.temperatures
        pressures = self.new_or_read.pressures
        return temperatures, pressures

    def _convert_pressure_units(self, pressures, path_arguments):
        if path_arguments is False:
            units = input("What was the unit of pressure input? Bar/kbar/GPa/MPa?\n")
        else:
            units = path_arguments[2]

        if units == 'GPa':
            pressures = pressures * 10000
        if units == 'MPa':
            pressures = pressures * 10
        if units == 'kbar':
            pressures = pressures * 1000

        return pressures

    def _create_f_path_dataframe(self, temperatures, pressures):
        if np.diff(temperatures)[-1] > 0 and np.diff(pressures)[-1] < 0:
            f_path = pd.DataFrame(
                [temperatures, pressures,
                 np.diff(temperatures, append=temperatures[-1]+1),
                 np.diff(pressures, append=pressures[-1]-1)]
            ).T
        elif np.diff(temperatures)[-1] <= 0 and np.diff(pressures)[-1] < 0:
            f_path = pd.DataFrame(
                [temperatures, pressures,
                 np.diff(temperatures, append=temperatures[-1]-1),
                 np.diff(pressures, append=pressures[-1]-1)]
            ).T
        else:
            f_path = pd.DataFrame(
                [temperatures, pressures,
                 np.diff(temperatures), np.diff(pressures)]
            ).T
        return f_path

    def _test_prograde_peak_retrograde(self, f_path):
        
        # Test for prograde, peak, and retrograde paths in the P-T path
        t_peak = f_path[0][(f_path[3] <= 0) & (f_path[2] > 0)]
        p_peak = f_path[1][(f_path[3] <= 0) & (f_path[2] > 0)]
        t_ret = f_path[0][(f_path[3] < 0) & (f_path[2] <= 0)]
        p_ret = f_path[1][(f_path[3] < 0) & (f_path[2] <= 0)]

        if not t_peak.empty:
            t_peak = np.array(t_peak)
        if not p_peak.empty:
            p_peak = np.array(p_peak)
        t_ret = np.array(t_ret)
        p_ret = np.array(p_ret)

        return t_peak, p_peak, t_ret, p_ret

    def _apply_interpolation_loop(self, f_path, pressures):
        peak_index = f_path.index[f_path[1] == f_path[1].max()]
        peak_index = peak_index[0]

        t_pro = list(f_path[0][0:peak_index+1])
        p_pro = list(f_path[1][0:peak_index+1])

        t_ret = list(f_path[0][peak_index:])
        p_ret = list(f_path[1][peak_index:])

        ius = InterpolatedUnivariateSpline(p_pro, t_pro)
        rbf = Rbf(p_ret, t_ret)

        p_line = np.linspace(min(pressures), max(pressures), 30)
        p_line2 = np.linspace(max(pressures), pressures[-1], 30)

        yi = ius(p_line)
        fi = rbf(p_line2)

        temperatures = list(np.around(yi, 2)) + list(np.around(fi, 2))
        pressures = list(np.around(p_line, 2)) + list(np.around(p_line2, 2))
        return temperatures, pressures

    def _apply_interpolation_retrograde(self, f_path, pressures):
        peak_index = f_path.index[f_path[1] == f_path[1].max()]
        peak_index = peak_index[0]

        t_ret = list(f_path[0][peak_index:])
        p_ret = list(f_path[1][peak_index:])

        rbf = Rbf(p_ret, t_ret)

        p_line = np.linspace(min(pressures), max(pressures), 30)
        p_line2 = np.linspace(max(pressures), pressures[-1], 30)

        fi = rbf(p_line2)

        temperatures = list(np.around(fi, 2))
        pressures = list(np.around(p_line2, 2))
        return temperatures, pressures

    def _choose_digitization_method(self, path_arguments):
        answers = ["new", "stored"]
        if path_arguments is False:
            for val in answers:
                print(val)
            answer = input("Pathfinder function - new or stored path? Select answer\n")
        else:
            answer = path_arguments[1]

        if answer == answers[0]:
            self.new_or_read.run()
        elif answer == answers[1]:
            self.new_or_read.stored_digitization()
        else:
            print("Unexpected end - no P-T file input")
            time.sleep(10)
            exit()

    def _store_image_pt_path(self, path_arguments):
        temperatures = self.new_or_read.temperatures
        pressures = self.new_or_read.pressures
        return temperatures, pressures

    def _convert_pressure_units(self, pressures, path_arguments):
        if path_arguments is False:
            units = input("What was the unit of pressure input? Bar/kbar/GPa/MPa?\n")
        else:
            units = path_arguments[2]

        if units == 'GPa':
            pressures = pressures * 10000
        if units == 'MPa':
            pressures = pressures * 10
        if units == 'kbar':
            pressures = pressures * 1000
        return pressures

    def _create_f_path_dataframe(self, temperatures, pressures):
        if np.diff(temperatures)[-1] > 0 and np.diff(pressures)[-1] < 0:
            f_path = pd.DataFrame(
                [temperatures, pressures,
                 np.diff(temperatures, append=temperatures[-1]+1),
                 np.diff(pressures, append=pressures[-1]-1)]
            ).T
        elif np.diff(temperatures)[-1] <= 0 and np.diff(pressures)[-1] < 0:
            f_path = pd.DataFrame(
                [temperatures, pressures,
                 np.diff(temperatures, append=temperatures[-1]-1),
                 np.diff(pressures, append=pressures[-1]-1)]
            ).T
        else:
            f_path = pd.DataFrame(
                [temperatures, pressures,
                 np.diff(temperatures), np.diff(pressures)]
            ).T
        return f_path

    def _extract_prograde_retrograde_paths(self, f_path):
        t_peak = f_path[0][f_path[3] <= 0][f_path[2] > 0]
        if len(t_peak.index) > 0:
            t_peak = f_path[0].loc[t_peak.index[0]:t_peak.index[-1]+1]
            t_peak = np.array(t_peak)
        p_peak = f_path[1][f_path[3] <= 0][f_path[2] > 0]
        if len(p_peak.index) > 0:
            p_peak = f_path[1].loc[p_peak.index[0]:p_peak.index[-1]+1]
            p_peak = np.array(p_peak)
        t_ret = f_path[0][f_path[3] < 0][f_path[2] <= 0]
        t_ret = np.array(t_ret)
        p_ret = f_path[1][f_path[3] < 0][f_path[2] <= 0]
        p_ret = np.array(p_ret)

        peak_index = f_path.index[f_path[1] == f_path[1].max()]
        peak_index = peak_index[0]

        t_pro = list(f_path[0][0:peak_index+1])
        p_pro = list(f_path[1][0:peak_index+1])

        t_ret = list(f_path[0][peak_index:])
        p_ret = list(f_path[1][peak_index:])
        return t_pro, p_pro, t_ret, p_ret

    def gridding(self, path_arguments, path_increment):
        """
        Grids the temperature and pressure arrays based on the given path arguments and increments.

        Args:
            path_arguments (list): List of path arguments, where the third element represents the pressure unit.
            path_increment (list): List of path increments, where the first element represents the pressure increment
                                   and the second element represents the temperature increment.

        Returns:
            None
        """

        self.new_or_read.stored_digitization()
        temperatures = self.new_or_read.temperatures
        pressures = self.new_or_read.pressures
        # transform pressure based on path_arguments to bar
        if path_arguments[2] == 'GPa':
            pressures = pressures * 10000
        if path_arguments[2] == 'kbar':
            pressures = pressures * 1000
        # round the pressure array based on the increment
        if np.float32(path_increment[0]) >=10:
            pressures = np.round(pressures,-1)
        # round the temperature array based on the increment
        if np.float32(path_increment[1]) >=10:
            temperatures = np.round(temperatures,-1)
        # Creating arrays for temperature and pressure based on input increments
        x = np.arange(min(temperatures), max(temperatures), np.int32(path_increment[1]))
        y = np.arange(min(pressures), max(pressures), np.int32(path_increment[0]))
        # generating mesh array
        xv, yv = np.meshgrid(x, y)
        # flatten the mesh array for node input
        temperatures = xv.flatten()
        pressures = yv.flatten()
        # Write infromation to function variables
        self.temp = temperatures
        self.pressure = pressures
        self.time_var = np.full(len(temperatures),np.nan)
        self.depth = np.full(len(temperatures),np.nan)
        self.sub_angle = "Undefined"
        self.plate_v = "Undefined"


class Pathfinder_calc:
    def __init__(self):
        self.rock_rho = [2800, 3300]
        self.geotherm = 25  # degree C/km
        self.timestep = 1000  # in years
        self.rate = 0.02  # m/year
        self.angle = 15  # degree
        self.T = []
        self.P = []
        self.Y_val = []
        self.t_start = []

    def calc_time_model(self, timestep=1000, t_end=33e6, start_depth=20, end_depth=80000, rate=1.5, angle=15):
        self.timestep = timestep
        self.t_end = t_end
        self.t_start = [144e6]  # default is 144_000_000 years (144 Ma)
        self.rate = rate / 100  # m/year
        self.angle = angle
        self.Y_val = [start_depth]
        self.T = [start_depth * self.geotherm + 200]  # Temperature in °C
        self.P = [(self.rock_rho[0] * 1000 + self.rock_rho[1] * (start_depth - 1000)) * 9.81]  # N/m2

        while self.t_start[-1] > self.t_end:
            Y1 = self.Y_val[-1]
            x_step = self.rate * self.timestep
            y_step = self.rate * self.timestep * abs(np.sin(self.angle / 180 * np.pi))
            x = self.Y_val[-1] + x_step
            y = self.Y_val[-1] + y_step

            temp_step = self.geotherm * (y - Y1)
            press_step = self.rock_rho[1] * (y - Y1) * 9.81
            T = self.T[-1] + temp_step
            P = self.P[-1] + press_step

            self.Y_val.append(y)
            self.T.append(T)
            self.P.append(P)
            self.t_start.append(self.t_start[-1] - self.timestep)

            if self.Y_val[-1] > end_depth:
                break

class Pathfinder:
    """
    Class representing a pathfinder for calculating temperature and pressure values along a path.

    Attributes:
    - temperature: The temperature values along the path.
    - pressure: The pressure values along the path.
    - time: The time values along the path.
    - depth: The depth values along the path.
    - metadata: Additional metadata associated with the path.
    - dt: The time step used for calculations.

    Methods:
    - connect_extern: Connects to external modules and performs calculations based on the selected mode.
    """

    def __init__(self):
        self.temperature = 0
        self.pressure = 0
        self.time = 0
        self.depth = 0
        self.metadata = {}
        self.dt = 0

        self.mod2 = PTPathDigitizer()
        self.mod3 = Pub_pathfinder()
        self.theoule = core_Pathfinder()

    def connect_extern(self, path_arguments=False, path_increment=False):
        """
        Connects to external modules and performs calculations based on the selected mode.

        Parameters:
        - path_arguments: Optional. List of path arguments.
        - path_increment: Optional. List of path increments.

        Returns:
        None
        """
        main_folder = Path(__file__).parent.absolute()
        file_to_open = main_folder / "output_Pathfinder.txt"

        if path_arguments is False:
            answer = input(
                "What mode do you want to use?\n[Mod1, Mod2, Mod3, Mod4, Mod5]")
        else:
            # Take stated answer from init file
            answer = path_arguments[0]

        # default setting if no path increments are given as input
        if path_increment is False:
            path_increment = [500, 15, 350]

        if answer == 'Mod1':
            self._handle_mod1()
        elif answer == 'Mod2':
            self._handle_mod2(path_arguments, path_increment)
        elif answer == 'Mod3':
            self._handle_mod3()
        elif answer == 'Mod4':
            self._handle_mod4(path_arguments, path_increment)
        elif answer == 'Mod5':
            self._handle_mod5(path_arguments, path_increment)
        elif answer == 'Mod6':
            self._handle_mod6(path_arguments, path_increment)

        # Path and Metadata
        # Create the data variable and generate output
        df = pd.DataFrame(
            [self.temperature, self.pressure, self.time, self.depth]).T
        df.columns = ['Temperature', 'Pressure', 'Time', 'Depth']

        with open('meta.json', 'w') as f:
            json.dump(self.metadata, f, indent=4, sort_keys=True)

        df.to_csv(file_to_open, sep=',', header=True, index=False)

    def _handle_mod1(self):
        # Subduction path
        dt = 1000
        time_end = 33e6
        depth_0 = int(input("Please enter a starting depth in meter..."))
        depth_1 = int(input(
            "Please enter the maximum depth for your model in meter..."))
        rate = float(input(
            "Please enter a convergence rate in m/year (e.g., 0.02)..."))
        angle = int(input("Please enter a subduction angle in degree..."))
        calc_path = Pathfinder_calc()
        calc_path.calc_time_model(
            timestep=dt, t_end=time_end, start_depth=depth_0, end_depth=depth_1, rate=rate, angle=angle)

        # store P and T values
        self.temperature = calc_path.T
        self.pressure = calc_path.P
        self.time = calc_path.t_start
        self.depth = calc_path.Y_val
        self.dt = dt

        # Store metadata
        self.metadata['Rock density [kg/m3]'] = calc_path.rock_rho
        self.metadata['Geotherm [degree C/km]'] = calc_path.geotherm
        self.metadata['Time step [year]'] = calc_path.timestep
        self.metadata['Sub./burial rate [km/year]'] = calc_path.rate
        self.metadata['Burial angle [Degree]'] = calc_path.angle
        self.metadata['Temperature unit'] = 'Degree C'
        self.metadata['Pressure unit'] = 'Pascal'

    def _handle_mod2(self, path_arguments, path_increment):
        """ Line path digitizing module"""

        # do linear path from digitized path and filter it by path increments
        if path_arguments is False:
            self.mod2.run()
        else:
            self.mod2.stored_digitization()

        
        _t_arr = np.linspace(min(self.mod2.temperatures), max(self.mod2.temperatures), 100)
        _p_arr = np.linspace(min(self.mod2.pressures), max(self.mod2.pressures), 100)

        # filter the path by increments using Pathfinder_Theoule module
        self.theoule = Pathfinder_Theoule(
                _t_arr, _p_arr, plate_speed=path_arguments[3],
                sub_angle=path_arguments[4], dt=np.ones(len(_p_arr)), path_increment=path_increment)
        self.theoule.filter_steps(_t_arr, [_p_arr, np.ones(len(_p_arr)), np.ones(len(_t_arr))])


        # store P and T values
        self.temperature = self.theoule.temp[::-1]
        self.pressure = self.theoule.pressure[::-1]
        self.depth = self.theoule.depth[::-1]
        self.dt = self.theoule.dt

        # Store metadata
        self.metadata['Convergence rate [cm/year]'] = path_arguments[3]
        self.metadata['Burial angle [Degree]'] = path_arguments[4]
        self.metadata['Temperature unit'] = 'Degree C'
        self.metadata['Pressure unit'] = 'Bar'
        self.metadata['Time step [years]'] = self.dt

    def _handle_mod3(self):
        # Plain digitizing module
        digitizer = PTPathDigitizer()
        digitizer.run()

        # store P and T values
        self.temperature = digitizer.temperatures
        self.pressure = digitizer.pressures

        # Store metadata
        self.metadata['Temperature unit'] = 'Degree C'
        self.metadata['Pressure unit'] = 'Bar'

    def _handle_mod4(self, path_arguments, path_increment):
        # Theoule mod for fitting a subduction rate to a digitized P-T path - only prograde
        if path_arguments is False:
            self.theoule.digitilization_module(path_increment=path_increment)
        else:
            self.theoule.digitilization_module(path_arguments, path_increment)

        # store P and T values
        self.temperature = self.theoule.temp
        self.pressure = self.theoule.pressure
        self.time = self.theoule.time_var
        self.depth = self.theoule.depth
        self.dt = self.theoule.dt

        # Store metadata
        self.metadata['Convergence rate [cm/year]'] = self.theoule.plate_v
        self.metadata['Burial angle [Degree]'] = self.theoule.sub_angle
        self.metadata['Temperature unit'] = 'Degree C'
        self.metadata['Pressure unit'] = 'Bar'
        self.metadata['Time step [years]'] = self.dt

    def _handle_mod5(self, path_arguments, path_increment):
        # Theoule mod for fitting a subduction rate to a digitized P-T path - loop for pro and retrograde
        self.theoule.loop_digi(path_arguments, path_increment)

        # store P and T values
        self.temperature = self.theoule.temp
        self.pressure = self.theoule.pressure
        self.time = self.theoule.time_var
        self.depth = self.theoule.depth
        self.dt = self.theoule.dt

        # Store metadata
        self.metadata['Convergence rate [cm/year]'] = self.theoule.plate_v
        self.metadata['Burial angle [Degree]'] = self.theoule.sub_angle
        self.metadata['Temperature unit'] = 'Degree C'
        self.metadata['Pressure unit'] = 'Bar'
        self.metadata['Time step [years]'] = self.dt

    def _handle_mod6(self, path_arguments, path_increment):
        self.theoule.gridding(path_arguments, path_increment)
        # store P and T values
        self.temperature = self.theoule.temp
        self.pressure = self.theoule.pressure
        self.time = self.theoule.time_var
        self.depth = self.theoule.depth
        self.dt = self.theoule.dt

        # Store metadata
        self.metadata['Convergence rate [cm/year]'] = self.theoule.plate_v
        self.metadata['Burial angle [Degree]'] = self.theoule.sub_angle
        self.metadata['Temperature unit'] = 'Degree C'
        self.metadata['Pressure unit'] = 'Bar'
        self.metadata['Time step [years]'] = self.dt



if __name__ == '__main__':

    nasa = Pathfinder()
    nasa.connect_extern()



