import numpy as np
import matplotlib.pyplot as plt

from mpl_toolkits.mplot3d import Axes3D

def plot_spheric_3D(ax:Axes3D, theta:np.ndarray, phi:np.ndarray, color_value:np.ndarray, cmap='plasma'):
    """
    Plot a 3D spherical surface.
    
    Parameters
    ----------
    ax : plt.Axes
        The axes to plot on.
    theta, phi : np.ndarray
        2D meshgrid of the angles theta and phi.
    color_value : np.ndarray
        values of the function to plot on the surface.
    """
    x = np.sin(theta) * np.cos(phi)
    y = np.sin(theta) * np.sin(phi)
    z = np.cos(theta)

    cmap = plt.cm.get_cmap(cmap)
    ax.plot_surface(x, y, z, facecolors=cmap(color_value), shade=True, rstride=1, cstride=1, linewidth=0)

    # Set the axis labels and title
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_aspect('equal')
