
"""
Written by
Thorsten Markmann
thorsten.markmann@geo.unibe.ch
status: 25.01.2024
"""

import pandas as pd
import numpy as np
from pathlib import Path
import time

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
    reply = False
    while not reply:
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
        validLength = False
        while not validLength:
            reflength = simpledialog.askfloat(
                "Enter reference length",
                "Enter the reference length in {:s} direction".format(direction))
            if isinstance(reflength, float):
                validLength = True
            else:
                messagebox.showerror("Error", "Please provide a valid length.")
        # calculate scaling factor
        deltaref = coord[1][index]-coord[0][index]
        factor = reflength/deltaref
        # ask for origin values of plot
        validOrigin = False
        while not validOrigin:
            base = simpledialog.askfloat(
                "Enter origin value",
                "Enter origin value {:s} direction".format(direction))
            if isinstance(base, float):
                validOrigin = True
            else:
                messagebox.showerror(
                    "Error", "Please provide a valid origin value.")
        reply = messagebox.askyesno(
            "Length confirmation",
            "You selected {:4.0f} pixels in {:s} direction"
            "corresponding to {:4.4f} units. Is this correct?".format(
                deltaref, direction, reflength)
        )
    base = float(base)
    return factor, coord[0][index], base

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

        # prepare spline from input P-T
        spl = splrep(self.pressure, self.temp)

        depth = 1000
        crust_d = 1000
        c_p_list = []
        # calculated pressure with 1km of starting cont. crust of density 2800
        c_p = (self.rho[0]*crust_d + self.rho[1]
               * (depth-crust_d)) * 9.81 / 10**5
        # c_p = self.rho[1] * depth * 1000 * 9.81 / 10**5
        d_step = self.speed * self.dt * abs(np.sin(self.angle/180*np.pi))

        while c_p < self.pressure[-1]:
            # calc pressure in bar
            c_p = (self.rho[0]*crust_d + self.rho[1]*(depth-1)) * 9.81 / 10**5
            # c_p = self.rho[1] * depth * 1000 * 9.81 / 10**5

            if c_p < self.pressure[0]:
                depth = depth + d_step
            else:
                depth = depth + d_step
                self.time.append(self.time[-1]+self.dt)
                self.depth.append(depth)
                c_p_list.append(c_p)

        print('End-depth is: {} km'.format(depth))

        yinterp = splev(c_p_list, spl)
        # plt.plot(c_p_list, yinterp, 'x--', markersize = 5)
        # plt.plot(self.pressure, self.temp, '-.')

        f_option = 3

        # filter option #1 - rough 300 bar filter, no caution for temperature
        if f_option == 1:
            while np.diff(c_p_list)[0] < self.p_increment:
                num = round(300/np.diff(c_p_list)[0])
                del c_p_list[::num]
                del self.time[::num]
                del self.depth[::num]
                yinterp = splev(c_p_list, spl)
            yinterp = splev(c_p_list, spl)

        if f_option == 2:
            # filter option #2 - min 300 bar and 10 DC steps
            newlist = [c_p_list[0]]
            for val in c_p_list:
                step = val - newlist[-1]
                if step >= self.p_increment:
                    newlist.append(val)
            c_p_list = np.array(newlist)
            yinterp = splev(c_p_list, spl)

            yinterp = list(yinterp)
            new_x = [yinterp[0]]
            new_y = [c_p_list[0]]
            for val in yinterp:
                step = val - new_x[-1]
                if step >= self.t_increment:
                    new_x.append(val)
                    index = yinterp.index(val)
                    new_y.append(newlist[index])

            yinterp = np.array(new_x)
            c_p_list = np.array(new_y)

        # option 3 - most convenient approach 13.02.2022
        # TODO: energy potential argument
        if f_option == 3:
            new_x = [yinterp[0]]
            new_y = [c_p_list[0]]
            new_d = [self.depth[0]]
            new_t = [self.time[0]]
            for i, val in enumerate(c_p_list):
                step_p = val - new_y[-1]
                step_t = yinterp[i] - new_x[-1]
                # define minimum pressure difference for step
                # NOTE "pro pressure steps?"
                if step_p >= self.p_increment:
                # if step_p >= 100 and step_t >= 1:
                    new_y.append(val)
                    new_x.append(yinterp[i])
                    new_d.append(self.depth[i])
                    new_t.append(self.time[i])
                # define minimum temperature difference for step
                # NOTE "pro temperature steps?"
                elif step_t >= self.t_increment:
                # elif step_t >= 1:
                    new_y.append(val)
                    new_x.append(yinterp[i])
                    new_d.append(self.depth[i])
                    new_t.append(self.time[i])

            yinterp = np.array(new_x)
            c_p_list = np.array(new_y)
            self.time = new_t
            self.depth = new_d

        # Selecting only steps with temperature >= 350 °C
        frame = pd.DataFrame([yinterp, c_p_list, self.time, self.depth])
        cut_T = self.lower_t_bound
        yinterp = np.array(frame.iloc[0][frame.iloc[0] >= cut_T])
        c_p_list = np.array(frame.iloc[1][frame.iloc[0] >= cut_T])
        self.time = np.array(frame.iloc[2][frame.iloc[0] >= cut_T])
        self.depth = np.array(frame.iloc[3][frame.iloc[0] >= cut_T])

        # test plot
        # plt.plot(c_p_list, yinterp, 'd', markersize=10)
        # plt.legend(['Spline points', 'Original', 'Filtered spline'])
        # plt.show()

        self.temp = yinterp
        self.pressure = c_p_list

    def loop(self):
        """
        Performs the burial and exhumation process to generate a P-T path.

        This method calculates the pressure-temperature (P-T) path by performing burial and exhumation steps.
        The burial step starts from the initial depth and increases the depth until the calculated pressure
        reaches the maximum pressure in the provided pressure array. The exhumation step starts from the
        maximum depth reached during burial and decreases the depth until the calculated pressure is greater
        than the minimum pressure in the provided pressure array.

        Returns:
            None
        """

        # NOTE "doing here different function for prograde, retrograde, loop and more?"
        # prepare spline from input P-T

        depth = 1000  # in meter
        crust_d = 1000  # in meter
        c_p_list = []
        # calculated pressure with 1km of starting cont. crust of density 2800
        c_p = (self.rho[0]*crust_d + self.rho[1] *
               (depth-crust_d)) * 9.81 / 10**5  # in Bar
        # c_p = self.rho[1] * depth * 1000 * 9.81 / 10**5
        d_step = self.speed/100 * self.dt * \
            abs(np.sin(self.angle/180*np.pi))  # in meter
        print("Start resampling into \x1B[3mP-T-t\x1B[0m path")
        # burial
        while c_p <= max(self.pressure):
            # calc pressure in bar
            c_p = (self.rho[0]*crust_d + self.rho[1] *
                   (depth-1)) * 9.81 / 10**5  # in Bar
            # c_p = self.rho[1] * depth * 1000 * 9.81 / 10**5

            if c_p < self.pressure[0]:
                depth = depth + d_step
            else:
                depth = depth + d_step
                self.time.append(self.time[-1]+self.dt)
                self.depth.append(depth)
                c_p_list.append(c_p)
        print("Burial finished...")
        # exhumation
        if c_p > self.pressure[-1]:
            # d_step = self.speed*3 * self.dt * abs(np.sin(80/180*np.pi))
            d_step = self.speed/100*5 * self.dt * \
                abs(np.sin(self.angle/180*np.pi))  # in meter
            while c_p > self.pressure[-1]:
                # calc pressure in bar
                c_p = (self.rho[0]*crust_d + self.rho[1] *
                       (depth-1)) * 9.81 / 10**5  # in Bar
                # c_p = self.rho[1] * depth * 1000 * 9.81 / 10**5

                if c_p < self.pressure[-1]:
                    depth = depth - d_step
                else:
                    depth = depth - d_step
                    self.time.append(self.time[-1]+self.dt)
                    self.depth.append(depth)
                    c_p_list.append(c_p)
            print("Exhumation finished...")
        print('Final depth is: {} m'.format(depth))

        # plt.plot(c_p_list, yinterp, 'x--', markersize = 5)
        # plt.plot(self.pressure, self.temp, '-.')
        # REVIEW slicing of time and depth probably crashes for other cases
        # It is adjusted for the Vho-extended-loop test (60 nodes from P and T)
        self.time = self.time[1:]
        self.depth

        # TODO loop function gives unfiltered data - no overstepping yet as for prograde version

        f_option = 3
        yinterp = self.temp
        c_p_list = self.pressure
        # option 3 - most convenient approach 13.02.2022
        if f_option == 3:
            new_x = [yinterp[0]]
            new_y = [c_p_list[0]]
            new_d = [self.depth[0]]
            new_t = [self.time[0]]
            for i, val in enumerate(c_p_list):
                step_p = val - new_y[-1]
                step_t = yinterp[i] - new_x[-1]
                # define minimum pressure difference for step
                # REVIEW "pressure steps?"
                if step_p >= 1000 or -step_p >= 1000:
                    new_y.append(val)
                    new_x.append(yinterp[i])
                    new_d.append(self.depth[i])
                    new_t.append(self.time[i])
                # define minimum temperature difference for step
                # REVIEW "temperature steps?"
                elif step_t >= 15 or -step_t >= 15:
                    new_y.append(val)
                    new_x.append(yinterp[i])
                    new_d.append(self.depth[i])
                    new_t.append(self.time[i])

            yinterp = np.array(new_x)
            c_p_list = np.array(new_y)
            self.time = new_t
            self.depth = new_d

        # Selecting only steps with temperature >= 350 °C
        frame = pd.DataFrame([yinterp, c_p_list, self.time, self.depth])
        cut_T = self.lower_t_bound
        yinterp = np.array(frame.iloc[0][frame.iloc[0] >= cut_T])
        c_p_list = np.array(frame.iloc[1][frame.iloc[0] >= cut_T])
        self.time = np.array(frame.iloc[2][frame.iloc[0] >= cut_T])
        self.depth = np.array(frame.iloc[3][frame.iloc[0] >= cut_T])

        # test plot
        # plt.plot(c_p_list, yinterp, 'd', markersize=10)
        # plt.legend(['Spline points', 'Original', 'Filtered spline'])
        # plt.show()

        self.temp = yinterp
        self.pressure = c_p_list


class Pathfinder_calc:
    """
    Creating very simple P-T-t path over subduction
    """

    def __init__(self):
        """
        Initialize timestep for subduction mechanism

        Args:
            timestep ([int]): [timestep for subduction eg. between 144 Ma and 44 Ma years, every 10000 years]
        """
        self.timestep = 0
        self.T = 0
        self.P = 0

        self.t_start = 0
        self.t_end = 0
        self.end_depth = 0
        self.rate = 0  # m/year
        self.angle = 0
        self.geotherm = 5/1000  # degree cel per m
        self.rock_rho = [2800, 3300]  # kg/m3
        self.X_val = [0]
        self.Y_val = []

    def calc_time_model(self,
                        timestep=1000, t_end=33e6, start_depth=20,
                        end_depth=80_000, t_start=144e6, rate=1.5, angle=15):
        """
        iterating over time and creating P-T-t path
        """
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

    def line_path(self, start_T=350, end_T=600, dT=10, start_p=5000, end_p=20000):
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


class Create_new_or_read_txt_pt_path:
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

        # give path to txt file for P-T path to open
        # FIXME - deactivated open file in stored digitization
        filein = filedialog.askopenfilename(
            title="Select a digitized path file",
            filetypes=[("other", ".txt")]
            )
        # FIXME - static path
        # filein = r"C:\Users\Markmann\PhD\Projects\cpag\Thorsten\DataFiles\CascadiaC_Condit2020.txt"
        # filein = r"C:\Users\Markmann\PhD\Projects\cpag\Thorsten\DataFiles\Syros_retro.txt"

        # for default file un-comment this:
        # filein = r"C:\Users\Markmann\PhD\Projects\cpag_code\Thorsten\DataFiles\Syros2.txt"

        with open(filein) as f:
            lines = f.readlines()

        temperatures = lines[0]
        temperatures = temperatures.split()
        for i, item in enumerate(temperatures):
            temperatures[i] = np.float32(item)

        pressures = lines[1]
        pressures = pressures.split()
        for i, item in enumerate(pressures):
            pressures[i] = np.float32(item)

        self.temperatures = np.array(temperatures)
        self.pressures = np.array(pressures)


class call_Pathfinder:
    """
    A class representing a Pathfinder object.

    Attributes:
        temp (int): The temperature value.
        pressure (int): The pressure value.
        time_var (int): The time variable.
        depth (int): The depth value.
        dt (int): The time step value.
        new_or_read (Create_new_or_read_txt_pt_path): An instance of Create_new_or_read_txt_pt_path class.

    Methods:
        __init__(self, temp=0, pressure=0, dt=1000): Initializes a Pathfinder object.
        execute_digi_prograde(self): Executes the digitization process for a prograde P-T path.
        execute_digi(self): Executes the digitization process for a P-T path.
        execute_digi_mod(self): Executes the modified digitization process for a P-T path.
        execute_digi_mod2(self, path_arguments=False, path_increment=False): Executes the second modified digitization process for a P-T path.
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
        self.new_or_read = Create_new_or_read_txt_pt_path()

    def execute_digi_prograde(self):
        """
        Executes the digitization process for a prograde P-T path.
        """
        # Code implementation...

    def execute_digi(self):
        """
        Executes the digitization process for a P-T path.
        """
        # Code implementation...

    def execute_digi_mod(self):
        """
        Executes the modified digitization process for a P-T path.
        """
        # Code implementation...

    def execute_digi_mod2(self, path_arguments=False, path_increment=False):
        """
        Executes the second modified digitization process for a P-T path.

        Args:
            path_arguments (bool or list, optional): The path arguments (default is False).
            path_increment (bool, optional): The path increment value (default is False).
        """
        # Code implementation...
class call_Pathfinder:
    """
    A class representing a Pathfinder object.

    Attributes:
        temp (int): The temperature value.
        pressure (int): The pressure value.
        time_var (int): The time variable.
        depth (int): The depth value.
        dt (int): The time step value.
        new_or_read (Create_new_or_read_txt_pt_path): An instance of Create_new_or_read_txt_pt_path class.
    """

    # External call to use Pathfinder module
    # Decide by input if you want to digitise a P-T path from image
    # or
    # use a stored P-T path from txt file

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
        self.new_or_read = Create_new_or_read_txt_pt_path()

    def execute_digi_prograde(self):
        """
        Executes the digitization process for a prograde P-T path.
        """
        answers = ["new", "stored"]

        # Choose new digitization or stored
        for val in answers:
            print(val)
        answer = input(
            "Pathfinder function - new or stored path? Select answer\n")
        if answer == answers[0]:
            self.PathfinderV2.run()
        elif answer == answers[1]:
            self.PathfinderV2.stored_digitization()
        else:
            exit()

        # Store image P-T path in array
        temperatures = self.PathfinderV2.temperatures
        pressures = self.PathfinderV2.pressures

        units = input(
            "What was the unit of pressure input? Bar/kbar/GPa/MPa?\n")
        if units == 'GPa':
            pressures = pressures * 10000
        if units == 'MPa':
            pressures = pressures * 10
        if units == 'kbar':
            pressures = pressures * 1000

        f_path = pd.DataFrame(
            [temperatures, pressures, np.diff(temperatures), np.diff(pressures)]
            ).T

        # testing for prograde, peak and retrograde
        t_peak = f_path[0][f_path[3] <= 0][f_path[2] > 0]
        p_peak = f_path[1][f_path[3] <= 0][f_path[2] > 0]
        t_ret = f_path[0][f_path[3] < 0][f_path[2] <= 0]
        p_ret = f_path[1][f_path[3] < 0][f_path[2] <= 0]

        # Theoule pathfinder creator
        nasa = Pathfinder_Theoule(
            temperatures, pressures,
            dt=self.dt
        )

        nasa.prograde()  # only for prograde P-T path

        # Update self variables
        self.temp = nasa.temp
        self.pressure = nasa.pressure
        self.time_var = nasa.time
        self.depth = nasa.depth

    def execute_digi(self):
        """
        Executes the digitization process for a P-T path.
        """
        answers = ["new", "stored"]

        # Choose new digitization or stored
        for val in answers:
            print(val)
        answer = input(
            "Pathfinder function - new or stored path? Select answer\n")
        if answer == answers[0]:
            self.PathfinderV2.run()
        elif answer == answers[1]:
            self.PathfinderV2.stored_digitization()
        else:
            exit()

        # Store image P-T path in array
        temperatures = self.PathfinderV2.temperatures
        pressures = self.PathfinderV2.pressures

        units = input(
            "What was the unit of pressure input? Bar/kbar/GPa/MPa?\n")
        if units == 'GPa':
            pressures = pressures * 10000
        if units == 'MPa':
            pressures = pressures * 10
        if units == 'kbar':
            pressures = pressures * 1000

        f_path = pd.DataFrame(
            [temperatures, pressures,
             np.diff(temperatures), np.diff(pressures)]
        ).T

        # testing for prograde, peak and retrograde
        t_peak = f_path[0][f_path[3] <= 0][f_path[2] > 0]
        p_peak = f_path[1][f_path[3] <= 0][f_path[2] > 0]
        t_ret = f_path[0][f_path[3] < 0][f_path[2] <= 0]
        p_ret = f_path[1][f_path[3] < 0][f_path[2] <= 0]

        # Theoule pathfinder creator
        nasa = Pathfinder_Theoule(
            temperatures, pressures,
            dt=self.dt
        )

        nasa.prograde()  # only for prograde P-T path

        # Update self variables
        self.temp = nasa.temp
        self.pressure = nasa.pressure
        self.time_var = nasa.time
        self.depth = nasa.depth

    # TODO writing argument for multiple peak and retrograde paths - complex loop
    # idea get pieces of prograde and retrograde path and the number of slices then iterate
    # is the BACK-UP that 'vho path' is using
    def execute_digi_mod(self):
        """
        Executes the modified digitization process for a P-T path.
        """
        pressures = self.pressure
        temperatures = self.temp

        # Consider prograde, peak and retrograde path units using differences
        # Create f_path dataframe
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

        # testing for prograde, peak and retrograde
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

        # REVIEW P-T-path loop flag!!!!
        loop = True

        if loop is True:
            peak_index = f_path.index[f_path[1] == f_path[1].max()]
            peak_index = peak_index[0]

            t_pro = list(f_path[0][0:peak_index+1])
            p_pro = list(f_path[1][0:peak_index+1])

            t_ret = list(f_path[0][peak_index:])
            p_ret = list(f_path[1][peak_index:])

            from scipy.interpolate import Rbf, InterpolatedUnivariateSpline

            # apply interpolation
            ius = InterpolatedUnivariateSpline(p_pro, t_pro)
            rbf = Rbf(p_ret, t_ret)

            # create array
            p_line = np.linspace(13000, 26000, 30)
            p_line2 = np.linspace(26000, 18000, 30)

            # Use interpolation
            yi = ius(p_line)
            fi = rbf(p_line2)

            # plt.plot(yi, p_line, 'xr-', fi, p_line2, 'xb-', temperatures, pressures, 'd')

            temperatures = list(np.around(yi, 2)) + list(np.around(fi, 2))
            pressures = list(np.around(p_line, 2)) + \
                list(np.around(p_line2, 2))

        # Theoule pathfinder creator
        nasa = Pathfinder_Theoule(
            temperatures, pressures,
            sub_angle=self.sub_angle,
            plate_speed=self.plate_v,
            dt=self.dt
        )
        if loop is True:
            # ANCHOR new loop routine for P-T path
            nasa.loop()  # first try to do a loop
        else:
            nasa.prograde()  # only for prograde P-T path

        # TODO double peak pressure or cut the first of retrograde path?
        # Update self variables
        self.temp = nasa.temp
        self.pressure = nasa.pressure
        self.time_var = nasa.time
        self.depth = nasa.depth

    # Deactivated loop - this is the up-to-date version for prograde digitized path modelling
    def execute_digi_mod2(self, path_arguments=False, path_increment=False):
        """
        Executes the second modified digitization process for a P-T path.

        Args:
            path_arguments (bool or list, optional): The path arguments (default is False).
            path_increment (bool, optional): The path increment value (default is False).
        """
        answers = ["new", "stored"]
        # Choose new digitization or stored
        # init inout or manual input
        if path_arguments is False:
            # manual
            for val in answers:
                print(val)
            answer = input(
                "Pathfinder function - new or stored path? Select answer\n")

            # Selected method from answer - 0 new digitisation - 1 stored txt
            if answer == answers[0]:
                self.new_or_read.run()
            elif answer == answers[1]:
                self.new_or_read.stored_digitization()
            else:
                print("Unexpected end - no P-T file input")
                time.sleep(10)
                exit()
        else:
            # init input
            answer = path_arguments[1]
            # Selected method from answer - 0 new digitisation - 1 stored txt
            if answer == answers[0]:
                self.new_or_read.run()
            elif answer == answers[1]:
                self.new_or_read.stored_digitization()
            else:
                print("Unexpected end - no P-T file input")
                time.sleep(10)
                exit()

        # Store image P-T path in array
        temperatures = self.new_or_read.temperatures
        pressures = self.new_or_read.pressures

        # convert units into Bar - THIS IS NOT SI units because of theriak input
        if path_arguments is False:
            # Manual
            units = input(
                "What was the unit of pressure input? Bar/kbar/GPa/MPa?\n")
        else:
            # init input
            units = path_arguments[2]

        if units == 'GPa':
            pressures = pressures * 10000
        if units == 'MPa':
            pressures = pressures * 10
        if units == 'kbar':
            pressures = pressures * 1000

        # Consider prograde, peak and retrograde path units using differences
        # Create f_path dataframe
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

        # testing for prograde, peak and retrograde
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

        # REVIEW P-T-path loop flag!!!!
        loop = False

        if loop is True:
            peak_index = f_path.index[f_path[1] == f_path[1].max()]
            peak_index = peak_index[0]

            t_pro = list(f_path[0][0:peak_index+1])
            p_pro = list(f_path[1][0:peak_index+1])

            t_ret = list(f_path[0][peak_index:])
            p_ret = list(f_path[1][peak_index:])

            # apply interpolation
            ius = InterpolatedUnivariateSpline(p_pro, t_pro)
            rbf = Rbf(p_ret, t_ret)

            # create array
            p_line = np.linspace(min(pressures), max(pressures), 30)
            p_line2 = np.linspace(max(pressures), pressures[-1], 30)

            # Use interpolation
            yi = ius(p_line)
            fi = rbf(p_line2)

            # plt.plot(yi, p_line, 'xr-', fi, p_line2, 'xb-', temperatures, pressures, 'd')

            temperatures = list(np.around(yi, 2)) + list(np.around(fi, 2))
            pressures = list(np.around(p_line, 2)) + \
                list(np.around(p_line2, 2))

        # Theoule pathfinder creator
        if path_arguments is False:
            nasa = Pathfinder_Theoule(
                temperatures, pressures,
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

        if loop is True:
            # ANCHOR new loop routine for P-T path
            nasa.loop()  # first try to do a loop
        else:
            nasa.prograde()  # only for prograde P-T path

        # TODO double peak pressure or cut the first of retrograde path?
        # Update self variables
        self.temp = nasa.temp
        self.pressure = nasa.pressure
        self.time_var = nasa.time
        self.depth = nasa.depth
        self.sub_angle = nasa.angle
        self.plate_v = nasa.speed

    # copy of execute_digi_mod2 - Active loop - solution for prograde plus retrograde modelling
    def loop_digi(self):
        """
        Performs the loop digitization process.

        This method prompts the user to choose between new digitization or stored path.
        It then stores the image P-T path in an array, converts the pressure units,
        considers prograde, peak, and retrograde path units, performs interpolation,
        and updates the self variables.

        Returns:
            None
        """
        answers = ["new", "stored"]
        # Choose new digitization or stored
        for val in answers:
            print(val)
        answer = input("Pathfinder function - new or stored path? Select answer\n")
        if answer == answers[0]:
            self.PathfinderV2.run()
        elif answer == answers[1]:
            self.PathfinderV2.stored_digitization()
        else:
            exit()
        # Store image P-T path in array
        temperatures = self.PathfinderV2.temperatures
        pressures = self.PathfinderV2.pressures

        # convert units into Bar - THIS IS NOT SI units because of theriak input
        units = input("What was the unit of pressure input? Bar/kbar/GPa/MPa?\n")
        if units == 'GPa':
            pressures = pressures * 10000
        if units == 'MPa':
            pressures = pressures * 10
        if units == 'kbar':
            pressures = pressures * 1000

        # Consider prograde, peak and retrograde path units using differences
        # Create f_path dataframe
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

        # testing for prograde, peak and retrograde
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

        # resembling loop
        peak_index = f_path.index[f_path[1] == f_path[1].max()]
        peak_index = peak_index[0]

        t_pro = list(f_path[0][0:peak_index+1])
        p_pro = list(f_path[1][0:peak_index+1])

        t_ret = list(f_path[0][peak_index:])
        p_ret = list(f_path[1][peak_index:])

        # apply interpolation
        ius = InterpolatedUnivariateSpline(p_pro, t_pro)
        rbf = Rbf(p_ret, t_ret)

        # create array
        p_line = np.linspace(min(pressures), max(pressures), 30)
        p_line2 = np.linspace(max(pressures), pressures[-1], 30)

        # Use interpolation
        yi = ius(p_line)
        fi = rbf(p_line2)

        # plt.plot(yi, p_line, 'xr-', fi, p_line2, 'xb-', temperatures, pressures, 'd')

        temperatures = list(np.around(yi, 2)) + list(np.around(fi, 2))
        pressures = list(np.around(p_line, 2)) + list(np.around(p_line2, 2))

        # Theoule pathfinder creator
        nasa = Pathfinder_Theoule(
            temperatures, pressures,
            sub_angle=self.sub_angle,
            plate_speed=self.plate_v,
            dt=self.dt
        )
        nasa.loop()  # first try to do a loop

        # TODO double peak pressure or cut the first of retrograde path?
        # Update self variables
        self.temp = nasa.temp
        self.pressure = nasa.pressure
        self.time_var = nasa.time
        self.depth = nasa.depth

    def simple_digi(self, path_arguments=False):
        """
        This method allows the user to choose between new digitization or stored path.
        If `path_arguments` is False, the user is prompted to select the method manually.
        If `path_arguments` is provided, the method uses the provided argument.
        The selected method is then executed using the `PathfinderV2` object.
        The image P-T path is stored in the `temp` and `pressure` arrays.
        """
        answers = ["new", "stored"]
        
        if path_arguments is False:
            # manual input
            for val in answers:
                print(val)
            answer = input("Pathfinder function - new or stored path? Select answer\n")

            if answer == answers[0]:
                self.PathfinderV2.run()
            elif answer == answers[1]:
                self.PathfinderV2.stored_digitization()
            else:
                print("Unexpected end - no P-T file input")
                time.sleep(10)
                exit()
        else:
            # init input
            answer = path_arguments[1]

            if answer == answers[0]:
                self.PathfinderV2.run()
            elif answer == answers[1]:
                self.PathfinderV2.stored_digitization()
            else:
                print("Unexpected end - no P-T file input")
                time.sleep(10)
                exit()

        self.temp = self.PathfinderV2.temperatures
        self.pressure = self.PathfinderV2.pressures

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

        self.mod2 = Create_new_or_read_txt_pt_path()
        self.mod3 = Pub_pathfinder()
        self.theoule = call_Pathfinder()

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
            # Subduction path
            dt = 1000
            # dt = input("Enter a value for the time steps in years (1000 is default)...")
            time_end = 33e6
            # depth_0 = 20_000  # in meter
            # depth_1 = 80_000  # in meter
            depth_0 = int(input("Please enter a starting depth in meter..."))
            depth_1 = int(input(
                "Please enter the maximum depth for you model in meter..."))
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

        if answer == 'Mod2':
            # Line path
            calc_path = Pathfinder_calc()
            calc_path.line_path()

            # store P and T values
            self.temperature = calc_path.T
            self.pressure = calc_path.P

            # Store metadata
            self.metadata['Temperature unit'] = 'Degree C'
            self.metadata['Pressure unit'] = 'Pascal'

        if answer == 'Mod3':
            # Plain digitizing module
            digitizer = Create_new_or_read_txt_pt_path()
            digitizer.run()

            # store P and T values
            self.temperature = digitizer.temperatures
            self.pressure = digitizer.pressures

            # Store metadata
            self.metadata['Temperature unit'] = 'Degree C'
            self.metadata['Pressure unit'] = 'Bar'

        if answer == 'Mod4':
            # Theoule mod for fitting a subduction rate to a digitized P-T path - only prograde
            if path_arguments is False:
                self.theoule.execute_digi_mod2(path_increment=path_increment)
            else:
                self.theoule.execute_digi_mod2(path_arguments, path_increment)

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

        if answer == 'Mod5':
            # Theoule mod for fitting a subduction rate to a digitized P-T path - loop for pro and retrograde
            self.theoule.loop_digi()

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

        if answer == 'Mod6':
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
            print()

        # Path and Metadata
        # Create the data variable and generate output
        df = pd.DataFrame(
            [self.temperature, self.pressure, self.time, self.depth]).T
        df.columns = ['Temperature', 'Pressure', 'Time', 'Depth']

        with open('meta.json', 'w') as f:
            json.dump(self.metadata, f, indent=4, sort_keys=True)

        df.to_csv(file_to_open, sep=',', header=True, index=False)


if __name__ == '__main__':

    nasa = Pathfinder()
    nasa.connect_extern()



