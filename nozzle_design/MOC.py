import numpy as np
import nozzle_geometry as geometry
import equations as eq
import constants as const

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

def calculate_x_y_coordinates(tan1, tan2, x1, y1, x2, y2):
    """
    Calculate the x and y coordinates of the point.

    Parameters:.
    

    Returns:
    x and y coordinates of the point.
    """
    # Solve for intersection of two lines: y - y1 = tan1*(x - x1), y - y2 = tan2*(x - x2)
    A = np.array([[-tan1, 1], [-tan2, 1]])
    b = np.array([y1 - tan1 * x1, y2 - tan2 * x2])
    x_p, y_p = np.linalg.solve(A, b)

    return x_p, y_p


def find_intersection_with_contour(x_contour, y_contour, x0, y0, slope):
    """
    Find the intersection point of a line (defined by point x0, y0 and slope)
    with a piecewise linear contour defined by x_contour, y_contour.

    Parameters:
    x_contour (array): x-coordinates of the nozzle wall contour.
    y_contour (array): y-coordinates of the nozzle wall contour.
    x0, y0 (float): Point through which the line passes.
    slope (float): Slope of the line (e.g., C+ characteristic slope).

    Returns:
    tuple: (x_int, y_int, wall_slope) intersection point and local slope of the wall at intersection
    """
    for i in range(len(x_contour) - 1):
        x1, x2 = x_contour[i], x_contour[i + 1]
        y1, y2 = y_contour[i], y_contour[i + 1]

        dx = x2 - x1
        dy = y2 - y1

        wall_slope = dy / dx

        if slope == wall_slope:
            continue  # Parallel lines, skip

        # Set the two equations equal and solve for x:
        # slope*(x - x0) + y0 = wall_slope*(x - x1) + y1
        # slope*x - slope*x0 + y0 = wall_slope*x - wall_slope*x1 + y1
        # (slope - wall_slope)*x = slope*x0 - wall_slope*x1 + y1 - y0
        x_int = (slope * x0 - wall_slope * x1 + y1 - y0) / (slope - wall_slope)
        y_int = slope * (x_int - x0) + y0

        # Check if intersection is within the segment [x1, x2]
        if (min(x1, x2) <= x_int <= max(x1, x2)):
            return x_int, y_int, wall_slope
        
    raise ValueError("No intersection found between the line and the wall contour.")



def MOC(z, x, y, theta_geometry, M0, P0, T0):
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
    x_p =  np.zeros(points)
    y_p =  np.zeros(points)
    C_minus = np.zeros(points)
    C_plus = np.zeros(points)

    # Initialize first point
    M[0] = M0
    P[0] = P0
    T[0] = T0
    x_p[0] = x[0]
    y_p[0] = y[0]
    
    nu[0] = eq.prandtl_meyer(const.GAMMA, M[0])
    mi[0] = eq.mach_angle(M[0])
    theta[0] = theta_geometry[0]
    
    init_delta_theta = np.max(theta_geometry) / z
    
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
                M[i] = eq.inverse_prandtl_meyer(const.GAMMA, nu[i], 'newton') #calculate the Mach number from the Prandtl-Meyer angle
                mi[i] = eq.mach_angle(M[i]) #calculate the Mach angle from the Mach number
                C_plus[i] = 0
                C_minus[i] = np.tan(((init_delta_theta*i)/2)-((mi[i]+eq.mach_angle(eq.inverse_prandtl_meyer(const.GAMMA,Q[i]/2,'newton')))/2)) 
                x_p[i], y_p[i] = calculate_x_y_coordinates(C_minus[i], C_plus[i], x[0], y[0], x[0], 0) #calculate the coordinates of the point
            else:
                axis_pos = axis.index(i) #this returns the position of the axis in the axis list
                Q[i] = Q[axis[axis_pos-1]+1] #the Q value is the same as the previous axis point plus one
                nu[i] = Q[i] #theta is zero so nu equal Q
                R[i] = theta[i] - nu[i] #same to say that R is minus nu
                M[i] = eq.inverse_prandtl_meyer(const.GAMMA, nu[i], 'newton') #calculate the Mach number from the Prandtl-Meyer angle
                mi[i] = eq.mach_angle(M[i]) #calculate the Mach angle from the Mach number
                C_plus[i] = 0
                C_minus[i] = np.tan((theta[axis[axis_pos-1]+1]/2)-((mi[i]+mi[axis[axis_pos-1]+1])/2))
                axis_pos = axis.index(i)
                prev_axis_idx = axis[axis_pos - 1] if axis_pos > 0 else None
                x_p[i], y_p[i] = calculate_x_y_coordinates(C_minus[i], C_plus[i], x_p[axis[axis_pos-1]+1], y_p[axis[axis_pos-1]+1], x_p[prev_axis_idx], y_p[prev_axis_idx]) #calculate the coordinates of the point
        elif i in wall:
            #theta[i] = theta_geometry[i] #rever esse aqui, pois não necessariamente o angulo do indice da parece é o mesmo da localizacao que estou analizando
            c_plus = np.tan(theta[i-1] + mi[i-1]) 
            x_p[i], y_p[i], theta[i] = find_intersection_with_contour(x, y, x_p[i-1], y_p[i-1], c_plus)
            R[i] = R[i-1]
            nu[i] = theta[i] - R[i]
            Q[i] = theta[i] + nu[i]
            M[i] = eq.inverse_prandtl_meyer(const.GAMMA, nu[i], 'newton')
            mi[i] = eq.mach_angle(M[i])
            C_minus[i] = 0
            #wall_pos = wall.index(i)
            #if wall_pos == 0: #first wall point
            #    C_plus[i] = np.tan((theta[i]+theta[i-1])/2 + ((mi[i]+mi[i-1])/2))
            #    x_p[i], y_p[i] = calculate_x_y_coordinates((theta[i]+theta[0])/2, C_plus[i], x[0], y[0], x_p[i-1], y_p[i-1])
            #else:
            #    prev_wall_idx = wall[wall_pos - 1] if wall_pos > 0 else None
            #    C_plus[i] = np.tan((theta[i]+theta[i-1])/2 + ((mi[i]+mi[i-1])/2))
            #    x_p[i], y_p[i] = calculate_x_y_coordinates((theta[i]+theta[prev_wall_idx])/2, C_plus[i], x_p[prev_wall_idx], y_p[prev_wall_idx], x_p[i-1], y_p[i-1])
            
        else: #internal point
            if i > z:
                axis_idx = max(j for j, idx in enumerate(axis) if idx <= i)
                Q[i] = Q[i-z+(axis_idx-1)]    

            R[i] = R[i-1]
            theta[i] = (Q[i] + R[i])/2
            nu[i] = Q[i] - theta[i]
            M[i] = eq.inverse_prandtl_meyer(const.GAMMA, nu[i], 'newton')
            mi[i] = eq.mach_angle(M[i])
            C_plus[i] = np.tan((theta[i]+theta[i-1])/2 + ((mi[i]+mi[i-1])/2))
            if i <=z:
                C_minus[i] = np.tan((init_delta_theta*i)/2 - ((mi[i]+eq.mach_angle(eq.inverse_prandtl_meyer(const.GAMMA,Q[i]/2,'newton')))/2))
                x_p[i], y_p[i] = calculate_x_y_coordinates(C_minus[i], C_plus[i], x[0], y[0],x_p[i-1], y_p[i-1])
            else:
                C_minus[i] = np.tan((theta[i]+theta[i-z+(axis_idx-1)])/2 - ((mi[i]+mi[i-z+(axis_idx-1)])/2))
                x_p[i], y_p[i] = calculate_x_y_coordinates(C_minus[i], C_plus[i], x_p[i-z+(axis_idx-1)], y_p[i-z+(axis_idx-1)],x_p[i-1], y_p[i-1])
            
    return nu, R, theta, Q, M, mi, x_p, y_p, C_minus, C_plus



#preciso importar a geometria do nozzle, definir o contorno, calcular as coordenadas dos pontos, a reta de cada linha característica, e aí o ponto na parede vai ser a interseção
# da reta do ponto anterior com o contorno do nozzle, com isso eu pego o theta e sigo o cálculo já implementado
#a inicializacao deve pegar o theta do primeiro angulo na garganta