import numpy as np
import nozzle_geometry as geometry
import equations as eq
import constants as const
import plot as plot

x,y, A,throat_area, thetas, throat_location_x, throat_location_y = geometry.initialize_nozzle_geometry("nozzle_design/nozzle-geometry.csv", 10, 3, only_divergent=False)

M = np.zeros(len(x))
P = np.zeros(len(x))
P0_P = np.zeros(len(x))
P#0_P[0] = eq.p0_p(const.GAMMA, 1.0)
p0 = 0.9433
#P[0] = p0 / P0_P[0]  # Pressure at the throat

for i in range(len(x)):
    supersonic = False
    if x[i] >= throat_location_x:
        supersonic = True
    M[i] = eq.solve_mach(A[i]/throat_area, gamma=const.GAMMA, supersonic=supersonic)    
    P0_P[i] = eq.p0_p(const.GAMMA, M[i])
    P[i] = p0 / P0_P[i]  # Pressure at the current point

y_lower = np.zeros(len(x))
plt = plot.plot_quasi1D_heatmap(x, y, y_lower, M, resolution=400, title="Quasi-1D Mach Number Heatmap", label="Mach Number")
plt.show()

plt = plot.plot_quasi1D_heatmap(x, y, y_lower, P, resolution=400, title="Quasi-1D Pressure Heatmap", label="Pressure")
plt.show()

plt = plot.plot_quasi1D_heatmap(x, y, y_lower, 1/P0_P, resolution=400, title="Quasi-1D P/P0 Heatmap", label="P/P0")
plt.show()