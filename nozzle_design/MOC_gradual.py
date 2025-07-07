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

def wall_indices(z,points):
    """
    Calculate the indices of the wall points in the MOC based on the number of initial divisions.

    Parameters:
    z (int): The number of initial divisions.

    Returns:
    list: A list of indices corresponding to the wall points.
    """
    indices = []
    indices.append(0)
    for i in range(points):
        ind = indices[i] + 2*z - 1
        indices.append(ind)
    
    return indices

def axis_indices(z,points):
    """
    Calculate the indices of the axis points in the MOC based on the number of initial divisions.

    Parameters:
    z (int): The number of initial divisions.

    Returns:
    list: A list of indices corresponding to the axis points.
    """
    indices = []
    indices.append(z-1)
    for i in range(points):
        ind = indices[i] + 2*z -1
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

def calculate_x_axis_coordinates(tan, x0, y0):
    """
    Calculate the x-coordinate of the point on the axis given a tangent and a point.

    Parameters:
    tan (float): The tangent of the angle at the point.
    x0 (float): The x-coordinate of the point.
    y0 (float): The y-coordinate of the point.

    Returns:
    float: The x-coordinate of the point on the axis.
    """
    # Solve for intersection with the x-axis: y = 0
    x_axis = (-y0) / tan + x0

    return x_axis, 0
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



def MOC(z, x, y, A, theta_geometry):
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
    points = len(x)*(z+1)
    
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

    # Initialize first points
    initial_A_ratio = A[1] / A[0]
    M0 = eq.solve_mach(initial_A_ratio, const.GAMMA)
    y_div = y[1]/(z-1)
    delta_theta = theta_geometry[1]/(z-1) #<<<<<<
    for i in range(z):
        M[i] = M0
        
        x_p[i] = x[1]
        y_p[i] = y[1] - (i * y_div)
        nu[i] = eq.prandtl_meyer(const.GAMMA, M[i])
        mi[i] = eq.mach_angle(M[i])
        theta[i] = theta_geometry[1] - (delta_theta * i)   # theta is zero in the axis
        R[i] = theta[i] - nu[i] 
        Q[i] = nu[i] + theta[i]
        C_minus[i] = np.tan(theta[i] - mi[i])
        C_plus[i] = np.tan(theta[i] + mi[i])
        print(f"i: {i}, x_p[i]: {x_p[i]}, y_p[i]: {y_p[i]}, theta[i]: {theta[i]}, nu[i]: {nu[i]}, mi[i]: {mi[i]}, R[i]:{R[i]}, Q[i]:{Q[i]} C_minus[i]: {C_minus[i]}, C_plus[i]: {C_plus[i]}")

    for i in range(z, points):
        
        print(f"i: {i-1}, x_p[i]: {x_p[i-1]}, y_p[i]: {y_p[i-1]}, theta[i]: {theta[i-1]}, nu[i]: {nu[i-1]}, mi[i]: {mi[i-1]}, C_minus[i]: {C_minus[i-1]}, C_plus[i]: {C_plus[i-1]}")
        if i in wall_indices(z, points):
            # Wall point
            c_plus = np.tan(theta[i-z+1] + mi[i-z+1]) #testar com o C_plus[i-z+1]
            print(f"i: {i}, z: {z}, x_p[i-z+1]: {x_p[i-z+1]}, y_p[i-z+1]: {y_p[i-z+1]}, c_plus: {c_plus}, [i-z+1]: {[i-z+1]}")
            x_p[i], y_p[i], theta[i] = find_intersection_with_contour(x, y, x_p[i-z+1], y_p[i-z+1], c_plus)
            R[i] = R[i-z+1]
            nu[i] = theta[i] + R[i]
            Q[i] = theta[i] + nu[i]
            M[i] = eq.inverse_prandtl_meyer(const.GAMMA, nu[i], 'newton')
            mi[i] = eq.mach_angle(M[i])
            C_minus[i] = 0
            C_plus[i] = 0
        elif i in axis_indices(z,points):
            #Axis point
            theta[i] = 0
            Q[i] = Q[i-z]
            nu[i] = Q[i]
            R[i] = nu[i] - theta[i]
            M[i] = eq.inverse_prandtl_meyer(const.GAMMA, nu[i], 'newton')
            mi[i] = eq.mach_angle(M[i])
            C_minus[i] = np.tan((theta[i-z]/2)-((mi[i]+mi[i-z])/2))
            C_plus[i] = np.tan((theta[i-z]/2)+((mi[i]+mi[i-z])/2))
            x_p[i], y_p[i] = calculate_x_axis_coordinates(C_minus[i], x_p[i-z], y_p[i-z])
        else:
            #internal point
            print("internal point: ",i)
            Q[i] = Q[i-z]
            R[i] = R[i-z-1]
            nu[i] = 0.5*(Q[i] + R[i])
            theta[i] = 0.5*(Q[i] - R[i])
            print(nu[i])
            M[i] = eq.inverse_prandtl_meyer(const.GAMMA, nu[i], 'newton')
            mi[i] = eq.mach_angle(M[i])
            C_minus[i] = np.tan(((theta[i]+theta[i-z])/2)-((mi[i]+mi[i-z])/2))
            C_plus[i] = np.tan(((theta[i]+theta[i-z+1])/2)+((mi[i]+mi[i-z+1])/2))
            x_p[i], y_p[i] = calculate_x_y_coordinates(C_minus[i], C_plus[i], x_p[i-z], y_p[i-z], x_p[i-z+1], y_p[i-z+1])
        if x_p[i-1] >= x[-1]:
            break
    return nu, R, theta, Q, M, mi, x_p, y_p, C_minus, C_plus