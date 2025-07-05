import matplotlib.pyplot as plt

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