import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from scipy.interpolate import griddata


def plot_nozzle_contour(x, y):
    """
    Plots the nozzle contour given x and y coordinates.

    Parameters:
    x (list): x-coordinates of the nozzle contour.
    y (list): y-coordinates of the nozzle contour.
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    plt.figure(figsize=(8, 4))
    plt.plot(x, y, label='Nozzle Contour')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.title('Nozzle Contour')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    return plt

def plot_points(x_p, y_p):
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(x_p[0], y_p[0], c="blue", marker="o", label="Origem", zorder=4)
    ax.scatter(x_p, y_p, c="black", marker="o", label="Pontos", zorder=4)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    plt.tight_layout()
    
    return plt

def plot_nozzle_and_points(x, y, x_p, y_p):
    """
    Plots the nozzle contour and MOC points on the same graph.

    Parameters:
    x (list): x-coordinates of the nozzle contour.
    y (list): y-coordinates of the nozzle contour.
    x_p (list): x-coordinates of the MOC points.
    y_p (list): y-coordinates of the MOC points.
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    if len(x_p) != len(y_p):
        raise ValueError("x_p and y_p must have the same length")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x, y, label='Nozzle Contour', color='tab:blue', linewidth=2, zorder=3)
    ax.scatter(x_p[0], y_p[0], c="red", marker="o", label="Origem", s=30, zorder=4)
    ax.scatter(x_p, y_p, c="black", marker="o", label="MOC Points", s=20, zorder=4)
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_title('Nozzle Contour with MOC Points')
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    return plt

def plot_mach_distribution(x, y, x_p, y_p, mach):
    """
    Plots the Mach number distribution along the nozzle using a colormap.

    Parameters:
    x (list): x-coordinates of the nozzle contour.
    y (list): y-coordinates of the nozzle contour.
    x_p (list): x-coordinates where Mach numbers are defined.
    y_p (list): y-coordinates where Mach numbers are defined.
    mach (list): Mach numbers at each (x_p, y_p) location.
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    if len(x_p) != len(y_p) or len(x_p) != len(mach):
        raise ValueError("x_p, y_p, and mach must have the same length")

    fig, ax = plt.subplots(figsize=(10, 6))
    # Plot nozzle contour
    ax.plot(x, y, label='Nozzle Contour', color='tab:blue', linewidth=2, zorder=2)
    # Scatter Mach numbers with colormap
    sc = ax.scatter(x_p, y_p, c=mach, cmap='viridis', vmin=1, vmax=max(mach), s=40, zorder=3)
    # Colorbar
    cbar = plt.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label('Mach Number')
    # Labels and title
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_title('Mach Number Distribution Along Nozzle')
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    return plt


def plot_data_heatmap_interpolated(x, y, x_p, y_p, data, resolution=300, label = "Mach number",title="Interpolated Mach Number Distribution (Heatmap)"):
    """
    Plota um heatmap interpolado da distribuição de Mach ao longo do bocal.

    Parameters:
    - x, y: coordenadas da parede do bocal
    - x_p, y_p: coordenadas dos pontos MOC
    - mach: valores de Mach nos pontos
    - resolution: resolução da malha de interpolação (quanto maior, mais suave)
    """
    # Gera malha regular no domínio dos pontos
    xi = np.linspace(min(x_p), max(x_p), resolution)
    yi = np.linspace(min(y_p), max(y_p), resolution)
    Xi, Yi = np.meshgrid(xi, yi)

    # Interpolação dos dados de Mach para a malha regular
    Mi = griddata((x_p, y_p), data, (Xi, Yi), method='linear', fill_value=np.nan)

    # Cria o gráfico
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plota o heatmap interpolado
    mask = ~np.isnan(Mi)
    cmap = plt.get_cmap('viridis')
    c = ax.contourf(Xi, Yi, Mi, levels=100, cmap=cmap)
    
    # Adiciona o contorno do bocal
    ax.plot(x, y, color='white', linewidth=2, label="Nozzle Contour")

    # Adiciona a barra de cores
    cbar = plt.colorbar(c, ax=ax)
    cbar.set_label(label)

    # Labels e formatação
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_title(title)
    ax.legend()
    ax.set_aspect(3)
    plt.tight_layout()
    return plt

import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
from scipy.interpolate import griddata


def plot_nozzle_contour(x, y):
    """
    Plots the nozzle contour given x and y coordinates.

    Parameters:
    x (list): x-coordinates of the nozzle contour.
    y (list): y-coordinates of the nozzle contour.
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    plt.figure(figsize=(8, 4))
    plt.plot(x, y, label='Nozzle Contour')
    plt.xlabel('x (mm)')
    plt.ylabel('y (mm)')
    plt.title('Nozzle Contour')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()

    return plt

def plot_points(x_p, y_p):
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.scatter(x_p[0], y_p[0], c="blue", marker="o", label="Origem", zorder=4)
    ax.scatter(x_p, y_p, c="black", marker="o", label="Pontos", zorder=4)
    ax.set_xlabel("x (mm)")
    ax.set_ylabel("y (mm)")
    plt.tight_layout()
    
    return plt

def plot_nozzle_and_points(x, y, x_p, y_p):
    """
    Plots the nozzle contour and MOC points on the same graph.

    Parameters:
    x (list): x-coordinates of the nozzle contour.
    y (list): y-coordinates of the nozzle contour.
    x_p (list): x-coordinates of the MOC points.
    y_p (list): y-coordinates of the MOC points.
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    if len(x_p) != len(y_p):
        raise ValueError("x_p and y_p must have the same length")
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(x, y, label='Nozzle Contour', color='tab:blue', linewidth=2, zorder=3)
    ax.scatter(x_p[0], y_p[0], c="red", marker="o", label="Origem", s=30, zorder=4)
    ax.scatter(x_p, y_p, c="black", marker="o", label="MOC Points", s=20, zorder=4)
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_title('Nozzle Contour with MOC Points')
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    return plt

def plot_mach_distribution(x, y, x_p, y_p, mach):
    """
    Plots the Mach number distribution along the nozzle using a colormap.

    Parameters:
    x (list): x-coordinates of the nozzle contour.
    y (list): y-coordinates of the nozzle contour.
    x_p (list): x-coordinates where Mach numbers are defined.
    y_p (list): y-coordinates where Mach numbers are defined.
    mach (list): Mach numbers at each (x_p, y_p) location.
    """
    if len(x) != len(y):
        raise ValueError("x and y must have the same length")
    if len(x_p) != len(y_p) or len(x_p) != len(mach):
        raise ValueError("x_p, y_p, and mach must have the same length")

    fig, ax = plt.subplots(figsize=(10, 6))
    # Plot nozzle contour
    ax.plot(x, y, label='Nozzle Contour', color='tab:blue', linewidth=2, zorder=2)
    # Scatter Mach numbers with colormap
    sc = ax.scatter(x_p, y_p, c=mach, cmap='viridis', vmin=1, vmax=max(mach), s=40, zorder=3)
    # Colorbar
    cbar = plt.colorbar(sc, ax=ax, pad=0.02)
    cbar.set_label('Mach Number')
    # Labels and title
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_title('Mach Number Distribution Along Nozzle')
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    return plt


def plot_data_heatmap_interpolated(x, y, x_p, y_p, data, resolution=300, label = "Mach number",title="Interpolated Mach Number Distribution (Heatmap)"):
    """
    Plota um heatmap interpolado da distribuição de Mach ao longo do bocal.

    Parameters:
    - x, y: coordenadas da parede do bocal
    - x_p, y_p: coordenadas dos pontos MOC
    - mach: valores de Mach nos pontos
    - resolution: resolução da malha de interpolação (quanto maior, mais suave)
    """
    # Gera malha regular no domínio dos pontos
    xi = np.linspace(min(x_p), max(x_p), resolution)
    yi = np.linspace(min(y_p), max(y_p), resolution)
    Xi, Yi = np.meshgrid(xi, yi)

    # Interpolação dos dados de Mach para a malha regular
    Mi = griddata((x_p, y_p), data, (Xi, Yi), method='linear', fill_value=np.nan)

    # Cria o gráfico
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Plota o heatmap interpolado
    mask = ~np.isnan(Mi)
    cmap = plt.get_cmap('viridis')
    c = ax.contourf(Xi, Yi, Mi, levels=100, cmap=cmap)
    
    # Adiciona o contorno do bocal
    ax.plot(x, y, color='white', linewidth=2, label="Nozzle Contour")

    # Adiciona a barra de cores
    cbar = plt.colorbar(c, ax=ax)
    cbar.set_label(label)

    # Labels e formatação
    ax.set_xlabel('x (mm)')
    ax.set_ylabel('y (mm)')
    ax.set_title(title)
    ax.legend()
    ax.set_aspect(3)
    plt.tight_layout()
    return plt


# New function: plot_quasi1D_mach_heatmap
def plot_quasi1D_heatmap(x, y_upper, y_lower, M, resolution=100, title="Quasi-1D Mach Number Heatmap", label="Mach Number"):
    """
    Plots a quasi-1D heatmap of Mach number distribution across a nozzle.

    Parameters:
    - x: array-like, x-coordinates along the nozzle centerline
    - y_upper: array-like, upper wall y-coordinates
    - y_lower: array-like, lower wall y-coordinates (typically zero)
    - M: array-like, Mach numbers at each x position
    - resolution: int, vertical resolution for the heatmap
    """
    x = np.array(x)
    y_upper = np.array(y_upper)
    y_lower = np.array(y_lower)
    M = np.array(M)

    # Create meshgrid
    xi = np.tile(x, (resolution, 1)).T  # Shape: (len(x), resolution)
    yi = np.linspace(0, 1, resolution)  # Normalized vertical coordinate

    # Scale each vertical slice according to local nozzle height
    nozzle_height = y_upper - y_lower
    yi_scaled = y_lower[:, None] + yi * nozzle_height[:, None]

    # Broadcast Mach values across vertical slices
    Mi = np.tile(M, (resolution, 1)).T

    # Plotting
    fig, ax = plt.subplots(figsize=(10, 4))
    c = ax.pcolormesh(xi, yi_scaled, Mi, cmap='viridis', shading='auto')

    # Plot wall contours
    ax.plot(x, y_upper, color='white', linewidth=2, label='Upper Wall')
    ax.plot(x, y_lower, color='white', linewidth=2, label='Lower Wall')

    # Labels and formatting
    cbar = plt.colorbar(c, ax=ax)
    cbar.set_label(label)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title(title)
    ax.legend()
    plt.tight_layout()

    return plt

