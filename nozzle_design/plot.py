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