import numpy as np

def num_points(z):
    """
    Calculate the number of points in the MOC based on the defined initial z divisions.

    Parameters:
    z (array): The number of initial divisions.

    Returns:
    int: The number of points in the MOC.
    """
    points = 0
    for i in range(z+1):
        points += i+1
    return points
