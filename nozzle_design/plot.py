import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict

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