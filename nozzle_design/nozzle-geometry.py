'''
This file will contain the nozzle geometry.

'''

import numpy as np

def read_nozzle_geometry(filename):
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    x = data[:, 0]
    A = data[:, 1]
    y = A/25.4
    return x, A, y


def calculate_throat_area(x, A):
    """
    Calculate the throat area of the nozzle.
    
    Parameters:
    x (array): The x-coordinates of the nozzle geometry.
    A (array): The cross-sectional areas at the corresponding x-coordinates.
    
    Returns:
    float: The throat area of the nozzle.
    """
    return np.min(A)

def calculate_thetas(x, y):
    """
    Calculate the angles of the nozzle at each point.
    
    Parameters:
    x (array): The x-coordinates of the nozzle geometry.
    y (array): The y-coordinates of the nozzle geometry at each x position.
    
    Returns:
    array: The angles of the nozzle at each point.
    """
    dy_dx = np.gradient(y, x)
    theta = np.arctan(dy_dx)
    return theta

def calculate_throat_location(x, A):
    """
    Find the location of the throat in the nozzle geometry.
    
    Parameters:
    x (array): The x-coordinates of the nozzle geometry.
    A (array): The cross-sectional areas at the corresponding x-coordinates.
    
    Returns:
    float: The x-coordinate of the throat location.
    """
    throat_index = np.argmin(A)
    return x[throat_index]

def get_only_divergent_section(x, y, throat_location, A):
    """
    Extract the divergent section of the nozzle geometry.
    
    Parameters:
    x (array): The x-coordinates of the nozzle geometry.
    y (array): The y-coordinates of the nozzle geometry at each x position.
    A (array): The cross-sectional areas at the corresponding x-coordinates.
    
    Returns:
    tuple: x-coordinates, y-coordinates, and areas of the divergent section.
    """
    throat_index = np.where(x == throat_location)[0][0]
    return x[throat_index:], y[throat_index:], A[throat_index:]

def initialize_nozzle_geometry(filename):
    """
    Initialize the nozzle geometry from a file.
    
    Parameters:
    filename (str): The name of the file containing the nozzle geometry data.
    
    Returns:
    tuple: x-coordinates and cross-sectional areas of the nozzle.
    """
    x, A, y = read_nozzle_geometry(filename)
    throat_area = calculate_throat_area(x, A)
    throat_location = calculate_throat_location(x, A)
    x, y, A = get_only_divergent_section(x, y, throat_location, A)
    thetas = calculate_thetas(x, y)
    
    return x,y, throat_area, thetas, throat_location


