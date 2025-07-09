import numpy as np
import nozzle_geometry as geometry
import equations as eq
import constants as const
import MOC_gradual as moc
import plot as plot

z = 5
x,y, A,throat_area, thetas, throat_location_x, throat_location_y = geometry.initialize_nozzle_geometry("nozzle_design/nozzle-geometry.csv", 10, 3)

nu, R, theta, Q, M, mi, x_p, y_p, C_minus, C_plus, P, T = moc.MOC(z, x, y, A, thetas, 0.9433, 297.65)

for i in range(len(nu)):
    print(f"Point {i}: x_p = {x_p[i]:.4f}, y_p = {y_p[i]:.4f}, \n\u03BD = {nu[i]:.4f}, R = {R[i]:.4f}, Q = {Q[i]:.4f}, \u03B8 = {theta[i]:.4f},  M = {M[i]:.4f}, \u03BC\u1D62 = {mi[i]:.4f}, p = {P[i]} C_minus = {C_minus[i]:.4f}, C_plus = {C_plus[i]:.4f}")


plt = plot.plot_mach_distribution(x, y, x_p[:-1], y_p[:-1], M[:-1])
plt.show()

plt = plot.plot_data_heatmap_interpolated(x, y, x_p[:-1], y_p[:-1], M[:-1],600, "Mach Number", "Interpolated Mach Number Distribution (Heatmap)")
plt.show()

plt = plot.plot_data_heatmap_interpolated(x, y, x_p[:-1], y_p[:-1], T[:-1],300, "Temperature", "Interpolated Temperature Distribution (Heatmap)")
plt.show()

plt = plot.plot_data_heatmap_interpolated(x, y, x_p[:-1], y_p[:-1], P[:-1],600, "Pressure", "Interpolated Pressure Distribution (Heatmap)")
plt.show()
