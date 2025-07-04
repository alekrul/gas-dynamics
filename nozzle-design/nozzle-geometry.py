'''
This file will contain the nozzle geometry.

'''

import numpy as np

def read_nozzle_geometry(filename):
    data = np.loadtxt(filename, delimiter=',', skiprows=1)
    x = data[:, 0]
    A = data[:, 1]
    return x, A

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

def calculate_thetas(x, A):
    """
    Calculate the angles of the nozzle at each point.
    
    Parameters:
    x (array): The x-coordinates of the nozzle geometry.
    A (array): The cross-sectional areas at the corresponding x-coordinates.
    
    Returns:
    array: The angles of the nozzle at each point.
    """
    dA_dx = np.gradient(A, x)
    theta = np.arctan(dA_dx / (2 * np.sqrt(A)))
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


def initialize_nozzle_geometry(filename):
    """
    Initialize the nozzle geometry from a file.
    
    Parameters:
    filename (str): The name of the file containing the nozzle geometry data.
    
    Returns:
    tuple: x-coordinates and cross-sectional areas of the nozzle.
    """
    x, A = read_nozzle_geometry(filename)
    throat_area = calculate_throat_area(x, A)
    thetas = calculate_thetas(x, A)
    throat_location = calculate_throat_location(x, A)
    
    return x, A, throat_area, thetas, throat_location