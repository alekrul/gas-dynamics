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

def wall_indices(z):
    """
    Calculate the indices of the wall points in the MOC based on the number of initial divisions.

    Parameters:
    z (int): The number of initial divisions.

    Returns:
    list: A list of indices corresponding to the wall points.
    """
    indices = []
    indices.append(0)
    for i in range(z):
        ind = indices[i] + (z-i)+1
        indices.append(ind)
    
    return indices

def axis_indices(z):
    """
    Calculate the indices of the axis points in the MOC based on the number of initial divisions.

    Parameters:
    z (int): The number of initial divisions.

    Returns:
    list: A list of indices corresponding to the axis points.
    """
    indices = []
    indices.append(1)
    for i in range(z-1):
        ind = indices[i] + (z-i)+1
        indices.append(ind)
    
    return indices

def MOC(z, theta_geometry, M0, P0, T0):
    """
    Perform the Method of Characteristics (MOC) for a nozzle flow.

    Parameters:
    z (int): The number of initial divisions.
    x (array): The x-coordinates of the nozzle geometry.
    A (array): The cross-sectional areas at the corresponding x-coordinates.
    theta (array): The angles of the nozzle at each point.
    M0 (float): Initial Mach number.
    P0 (float): Initial pressure.
    T0 (float): Initial temperature.
    gamma (float): Ratio of specific heats.
    R (float): Specific gas constant.

    Returns:
    tuple: Arrays containing the Mach number, pressure, and temperature at each point in the MOC.
    """
    
    points = num_points(z)
    
    M = np.zeros(points)
    P = np.zeros(points)
    T = np.zeros(points)
    
    R = np.zeros(points)
    Q = np.zeros(points)
    nu = np.zeros(points)
    mi = np.zeros(points)
    theta = np.zeros(points)


    # Initialize first point
    M[0] = M0
    P[0] = P0
    T[0] = T0
    
    init_delta_theta = theta_geometry[0] / z
    
    wall = wall_indices(z)
    axis = axis_indices(z)
    # Calculate characteristics
    for i in range(1, points):        
        #initialize first points with Q from point 0
        if i <= z:
            Q[i] = 2*init_delta_theta*i

        if i in axis:
            theta[i] = 0 #theta is zero in the axis
            if i == 1:
                nu[i] = Q[i] #theta is zero so nu equal Q
                R[i] = theta[i] - nu[i] #same to say that R is minus nu
            else:
                axis_pos = axis.index(i) #this returns the position of the axis in the axis list
                Q[i] = Q[axis[axis_pos-1]+1] #the Q value is the same as the previous axis point plus one
                nu[i] = Q[i] #theta is zero so nu equal Q
                R[i] = theta[i] - nu[i] #same to say that R is minus nu
        elif i in wall:
            theta[i] = theta_geometry[i] #rever esse aqui, pois não necessariamente o angulo do indice da parece é o mesmo da localizacao que estou analizando
            R[i] = R[i-1]
            nu[i] = theta[i] - R[i]
            Q[i] = theta[i] + nu[i]
        else: #internal point
            if i > z:
                axis_idx = max(j for j, idx in enumerate(axis) if idx <= i)
                Q[i] = Q[i-z+(axis_idx-1)]    

            R[i] = R[i-1]
            theta[i] = (Q[i] + R[i])/2
            nu[i] = Q[i] - theta[i]

    return nu, R, theta, Q



z = 4
theta_geometry = [13.19,0,0,0,0,9.8925,0,0,0,6.595,0,0,3.2975,0,0]
nu,R,theta,Q = MOC(z, theta_geometry, 1.0, 1, 288)

for i in range(len(nu)):
    print(f"Point {i}: nu = {nu[i]:.4f}, R = {R[i]:.4f}, theta = {theta[i]:.4f}, Q = {Q[i]:.4f}")
