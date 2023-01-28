# Libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl



def f(X,Y): # f(x,y,0)
    return 3-X**2-Y**2+np.sin(X*Y)


def simulation_heat_eq(rect,hs,BC, c, frames, eps = 1e-10):

    """ Simulation of the heat equation on a
    2D rectangular grid using the Finite Difference Method """

    X, Y = np.meshgrid(np.linspace(rect[0],rect[1],hs[0]),
                       np.linspace(rect[2],rect[3],hs[1])) # 2D meshgrid

    # Initial Conditions
    Z_init = f(X,Y) # f(x,y,0)
    zmax = max(Z_init.max(), BC[0], BC[1], BC[2], BC[3])
    zmin = min(Z_init.min(), BC[0], BC[1], BC[2], BC[3])

    # Boundary Conditions
    Z_init[0] = np.ones(len(Z_init[0])) * BC[0]
    Z_init[-1] = np.ones(len(Z_init[-1]))* BC[1]
    Z_init[:,0] =  np.ones(len(Z_init[:,0]))* BC[2]
    Z_init[:, -1] = np.ones(len(Z_init[:, -1])) * BC[3]

    # Figure settings
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.axes.set_xlim3d(rect[0] + eps, rect[1] - eps)
    ax.axes.set_ylim3d(rect[2] + eps, rect[3] - eps)
    ax.axes.set_zlim3d(zmin, zmax)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'

    # Plots for lines constituting the rectangle
    ax.plot([rect[0], rect[1]], [rect[2], rect[2]], [BC[0], BC[0]], color='black', linewidth=2)
    ax.plot([rect[1], rect[1]], [rect[2], rect[3]], [BC[3], BC[3]], color='black', linewidth=2)
    ax.plot([rect[0], rect[0]], [rect[2], rect[3]], [BC[2], BC[2]], color='black', linewidth=2)
    ax.plot([rect[0], rect[1]], [rect[3], rect[3]], [BC[1], BC[1]], color='black', linewidth=2)

    ax.plot([rect[0], rect[0]], [rect[2], rect[2]], [BC[0], BC[2]], color='black', linewidth=2)
    ax.plot([rect[1], rect[1]], [rect[2], rect[2]], [BC[0], BC[3]], color='black', linewidth=2)
    ax.plot([rect[0], rect[0]], [rect[3], rect[3]], [BC[1], BC[2]], color='black', linewidth=2)
    ax.plot([rect[1], rect[1]], [rect[3], rect[3]], [BC[1], BC[3]], color='black', linewidth=2)
    ax.set_title('Temperature development in a \n'
                 'rectangular room', fontsize = 18, fontname = 'STIXGeneral')

    # Infitesimals
    dx = (rect[1] - rect[0]) / (hs[0] - 1)
    dt = dx**2/(5*c**2)

    surf = ax.plot_surface(X, Y, Z_init, alpha=0.7, cmap=plt.cm.jet, vmin=zmin, vmax=zmax)
    cbar = fig.colorbar(surf)
    cbar.ax.set_ylabel('Temperature in ' + r'$^\circ\mathrm{C}$', rotation=270,fontsize = 14, labelpad=20)
    ax.set_axis_off()

    # Finite Difference Method
    for iteration in range(frames):

        surf.remove()
        for row in range(1, hs[0]-1):
            for col in range(1,hs[1]-1):
                Z_init[row,col] = Z_init[row,col] + c**2 * dt / (dx)**2  *(Z_init[row+1,col] - 4*Z_init[row,col]
                                               + Z_init[row-1,col] + Z_init[row,col+1]
                                               + Z_init[row,col-1])

        # Plotting surface with slight pause
        surf = ax.plot_surface(X, Y, Z_init, alpha = 0.7, cmap =plt.cm.jet, vmin = zmin, vmax = zmax)
        plt.pause(dt * 100)

simulation_heat_eq(rect = [-5,5,-5,5], hs = [25,25], BC = [-20,0,20,0], c = 5, frames = 1000)




