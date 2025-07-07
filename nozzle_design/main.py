import numpy as np
import nozzle_geometry as geometry
import equations as eq
import constants as const
import MOC_gradual as moc
import plot as plot

z = 6
x,y, A,throat_area, thetas, throat_location_x, throat_location_y = geometry.initialize_nozzle_geometry("nozzle_design/nozzle-geometry.csv")

#theta_geometry = [13.19,0,0,0,0,9.8925,0,0,0,6.595,0,0,3.2975,0,0]
#for i in range(len(theta_geometry)):
#    theta_geometry[i] = eq.deg2rad(theta_geometry[i])  # Convert angles to radians
#x = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0]
#y = [1,0,0,0,0,0,0,0,0,0,0,0,0,0,0]

nu, R, theta, Q, M, mi, x_p, y_p, C_minus, C_plus = moc.MOC(z, x, y, A, thetas)

for i in range(len(nu)):
    print(f"Point {i}: nu = {nu[i]:.4f}, R = {R[i]:.4f}, theta = {theta[i]:.4f}, Q = {Q[i]:.4f}, M = {M[i]:.4f}, mi = {mi[i]:.4f}, x_p = {x_p[i]:.4f}, y_p = {y_p[i]:.4f}, C_minus = {C_minus[i]:.4f}, C_plus = {C_plus[i]:.4f}")

plt = plot.plot_nozzle_and_points(x, y, x_p, y_p)
plt.show()

plt = plot.plot_mach_distribution(x, y, x_p, y_p, M)
plt.show()