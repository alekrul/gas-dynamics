import numpy as np
import nozzle_geometry as geometry
import equations as eq
import constants as const
import MOC as moc
import plot as plot

z = 1000
x,y, throat_area, thetas, throat_location_x, throat_location_y = geometry.initialize_nozzle_geometry("nozzle_design/nozzle-geometry.csv")

#theta_geometry = [13.19,0,0,0,0,9.8925,0,0,0,6.595,0,0,3.2975,0,0]
#for i in range(len(theta_geometry)):
#    theta_geometry[i] = eq.deg2rad(theta_geometry[i])  # Convert angles to radians
#x = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
#y = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

def find_intersection_with_contour(x0, y0, slope, x_wall, y_wall):
    for i in range(len(x_wall)-1):
        x1, y1 = x_wall[i], y_wall[i]
        x2, y2 = x_wall[i+1], y_wall[i+1]
        wall_slope = (y2 - y1) / (x2 - x1) if (x2 - x1) != 0 else np.inf

        # Parametrize both lines and compute intersection
        denom = (x2 - x1) * slope - (y2 - y1)
        if denom == 0:
            continue  # Parallel

        # intersection from parametric form of lines
        t = ((x1 - x0) * slope - (y1 - y0)) / denom
        x_int = x1 + t * (x2 - x1)
        y_int = y1 + t * (y2 - y1)

        # Check if intersection lies on the wall segment only
        if min(x1, x2) <= x_int <= max(x1, x2) and min(y1, y2) <= y_int <= max(y1, y2):
            return x_int, y_int, wall_slope
    return None, None, None

nu, R, theta, Q, M, mi, x_p, y_p, C_minus, C_plus = moc.MOC(z, x, y, thetas, 1.0, 1, 288)

for i in range(len(nu)):
    print(f"Point {i}: nu = {nu[i]:.4f}, R = {R[i]:.4f}, theta = {theta[i]:.4f}, Q = {Q[i]:.4f}, M = {M[i]:.4f}, mi = {mi[i]:.4f}, x_p = {x_p[i]:.4f}, y_p = {y_p[i]:.4f}, C_minus = {C_minus[i]:.4f}, C_plus = {C_plus[i]:.4f}")

plt = plot.plot_nozzle_and_points(x, y, x_p, y_p)
plt.show()
