# Libraries
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib as mpl

def f(X,Y): # u(x,y,0) = f(x,y,0)
    return X**2+Y**2

def g(X,Y): # u_t(x,y,0) = g(x,y,0)
    return -np.exp(Y/1000)

def laplacian(arr, row,col, dx ): # used to find initial condition for Wave eq.
    return (arr[row + 1, col]
                            + arr[row - 1, col] + arr[row, col + 1]
                            + arr[row, col - 1]- 4*arr[row,col])/(dx**2)

def simulation_pdes(rect,hs,BC, c, frames, eq, eps = 1e-10):

    """ Simulation of the heat and wave equation on a
    2D rectangular grid using the Finite Difference Method """

    X, Y = np.meshgrid(np.linspace(rect[0],rect[1],hs[0]),
                       np.linspace(rect[2],rect[3],hs[1])) # 2D meshgrid

    # Initial Conditions
    Z_init = f(X,Y) # u(x,y,0) = f(x,y,0)
    Z_init_2 = f(X,Y)
    Z_dot_init = g(X,Y) # u_t(x,y,0) = g(x,y,0)
    zs = []

    zmax = max(Z_init.max(), BC[0], BC[1], BC[2], BC[3])
    zmin = min(Z_init.min(), BC[0], BC[1], BC[2], BC[3])

    # Boundary Conditions
    if eq == "Heat":
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

    # Lines constituting the rectangle
    ax.plot([rect[0], rect[1]], [rect[2], rect[2]], [BC[0], BC[0]], color='black', linewidth=2)
    ax.plot([rect[1], rect[1]], [rect[2], rect[3]], [BC[3], BC[3]], color='black', linewidth=2)
    ax.plot([rect[0], rect[0]], [rect[2], rect[3]], [BC[2], BC[2]], color='black', linewidth=2)
    ax.plot([rect[0], rect[1]], [rect[3], rect[3]], [BC[1], BC[1]], color='black', linewidth=2)

    ax.plot([rect[0], rect[0]], [rect[2], rect[2]], [BC[0], BC[2]], color='black', linewidth=2)
    ax.plot([rect[1], rect[1]], [rect[2], rect[2]], [BC[0], BC[3]], color='black', linewidth=2)
    ax.plot([rect[0], rect[0]], [rect[3], rect[3]], [BC[1], BC[2]], color='black', linewidth=2)
    ax.plot([rect[1], rect[1]], [rect[3], rect[3]], [BC[1], BC[3]], color='black', linewidth=2)

    # Infitesimals
    dx = (rect[1] - rect[0]) / ((hs[0] - 1))
    if eq == 'Heat':
        dt = dx**2/(10*c**2)
    else:
        dt = dx/(10*c)

    # Case separation for titles and labels
    if eq == "Heat":
        ax.set_title('Temperature development in a \n'
                     'rectangular room', fontsize=18, fontname='STIXGeneral')
        surf = ax.plot_surface(X, Y, Z_init, alpha=0.7, cmap=plt.cm.jet, vmin=zmin, vmax=zmax)
        cbar = fig.colorbar(surf)
        cbar.ax.set_ylabel('Temperature in ' + r'$^\circ\mathrm{C}$', rotation=270,fontsize = 14, labelpad=20)
        ax.set_axis_off()

    if eq == "Wave":
        ax.set_title('Simulation of water waves in a \n'
                     'rectangular room', fontsize=18, fontname='STIXGeneral')
        surf = ax.plot_surface(X, Y, Z_init, alpha=0.7, cmap="magma", vmin=zmin, vmax=zmax)
        cbar = fig.colorbar(surf)
        cbar.ax.set_ylabel('Water height in meters', rotation=270, fontsize=14, labelpad=20)
        ax.set_axis_off()

    # Finite Difference Method
    if eq == 'Wave':
        for row in range(1, hs[0] - 1):
            for col in range(1, hs[1] - 1):
                Z_init_2[row,col] = Z_init_2[row,col] - 2 * laplacian(Z_dot_init,row,col, dx) + 1/2 * c ** 2 * dt ** 2 / (
                    dx) ** 2 * (Z_init_2[row+1, col] - 4 * Z_init_2[row, col]
                                + Z_init_2[row-1, col] + Z_init_2[row, col+1]
                                + Z_init_2[row, col-1])

        zs.append(Z_init_2)
        zs.append(Z_init)

        for iteration in range(2,frames):

            surf.remove()
            # Creating temporary matrix
            Z_temp = a = np.zeros((hs[0], hs[1]))
            Z_temp[0] = np.ones(len(Z_temp[0])) * BC[0]
            Z_temp[-1] = np.ones(len(Z_temp[-1])) * BC[1]
            Z_temp[:, 0] = np.ones(len(Z_temp[:, 0])) * BC[2]
            Z_temp[:, -1] = np.ones(len(Z_temp[:, -1])) * BC[3]

            for row in range(1, hs[0]-1):
                for col in range(1,hs[1]-1):
                    Z_temp[row,col] += 2 * zs[iteration - 1][row, col] - zs[iteration - 2][row, col] + c ** 2 * dt ** 2 / (
                            dx**2) * (zs[iteration - 1][row + 1, col] - 4 * zs[iteration - 1][row, col]
                                        + zs[iteration - 1][row - 1, col] + zs[iteration - 1][row, col + 1]
                                        + zs[iteration - 1][row, col - 1])

            zs.append(Z_temp)
            surf = ax.plot_surface(X, Y, zs[iteration-2], alpha = 1, cmap ="magma", vmin = zmin, vmax = zmax)
            plt.pause(dt)

    if eq == 'Heat':

        for iteration in range(frames):

            surf.remove()
            for row in range(1, hs[0]-1):
                for col in range(1,hs[1]-1):
                    Z_init[row,col] = Z_init[row,col] + c**2 * dt / (dx)**2  *(Z_init[row+1,col] - 4*Z_init[row,col]
                                                   + Z_init[row-1,col] + Z_init[row,col+1]
                                                   + Z_init[row,col-1])

            # Plotting surface with slight pause
            surf = ax.plot_surface(X, Y, Z_init, alpha = 1, cmap =plt.cm.jet, vmin = zmin, vmax = zmax)
            plt.pause(dt)

simulation_pdes(rect = [-1,2,-5,1], hs = [40,40], BC = [0,0,0,0], c = 5, frames = 1000, eq = 'Wave')






