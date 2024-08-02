import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

def f(r, theta):  # u(r, theta, 0) = f(r, theta, 0)
    return np.sqrt(r)*np.sin(theta)

def g(r, theta):  # u_t(r, theta, 0) = g(r, theta, 0)
    return np.sin(theta)

def laplacian_polar(arr, r, row, col, dr, dtheta):
    term1 = (arr[row + 1, (col) % (arr.shape[1])] - 2 * arr[row, (col) % (arr.shape[1])] + arr[row - 1, (col) % (arr.shape[1])]) / dr**2
    term2 = (arr[row + 1, (col) % (arr.shape[1])] - arr[row - 1, (col) % (arr.shape[1])]) / (2 * r[row] * dr)
    term3 = (arr[row, (col + 1) % (arr.shape[1])] - 2 * arr[row, (col) % (arr.shape[1])] + arr[row, (col - 1) % (arr.shape[1])]) / (r[row]**2 * dtheta**2)
    return term1 + term2 + term3

def simulation_pdes_polar(r_bounds, theta_bounds, hs, BC, BC_THETA, c, frames, eq):
    r_inner, r_outer = r_bounds
    r = np.linspace(r_inner, r_outer, hs[0], endpoint=True)
    theta = np.linspace(theta_bounds[0], theta_bounds[1], hs[1], endpoint=True)
    R, Theta = np.meshgrid(r, theta, indexing='ij')
    X = R * np.cos(Theta)
    Y = R * np.sin(Theta)

    # Initial Conditions
    Z_init = f(R, Theta)
    Z_0 = f(R, Theta)
    Z_dot_init = g(R, Theta)
    zs = []

    # Max/min for colorbar
    zmax = max(Z_init.max(), BC[0], BC[1], BC_THETA[0], BC_THETA[1])
    zmin = min(Z_init.min(), BC[0], BC[1], BC_THETA[0], BC_THETA[1])

    print(Z_init)

    # Boundary Conditions
    if eq == "Heat":
        Z_init[0,:] = BC[0]  # Inner radius
        Z_init[-1,:] = BC[1]  # Outer radius
        if theta_bounds[0] !=0 or theta_bounds[1] != 2*np.pi:
            Z_init[:,0] = BC_THETA[0]  # theta = theta_bounds[0]
            Z_init[:,-1]= BC_THETA[1]  # theta = theta_bounds[1]

        # boundary on theta_bounds[0]
        # boudnary on theta_bounds[1]
        # where is Z_init has theta = 0

    

    # Figure settings
    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.set_zlim(zmin, zmax)
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['font.family'] = 'STIXGeneral'


    # Boundary outline
    theta_full = np.linspace(theta_bounds[0], theta_bounds[1], 100)
    r_inner_x = r_inner * np.cos(theta_full)
    r_inner_y = r_inner * np.sin(theta_full)
    r_outer_x = r_outer * np.cos(theta_full)
    r_outer_y = r_outer * np.sin(theta_full)
    ax.plot(r_inner_x, r_inner_y, BC[0], color='black', linewidth = 5)
    ax.plot(r_outer_x, r_outer_y, BC[1], color='black', linewidth = 5)

    # Radial lines
    radial_x = [r_inner * np.cos(theta_bounds[0]), r_outer * np.cos(theta_bounds[0])]
    radial_y = [r_inner * np.sin(theta_bounds[0]), r_outer * np.sin(theta_bounds[0])]
    ax.plot(radial_x, radial_y, [BC[0], BC[1]], color='black', linewidth = 5)

    radial_x = [r_inner * np.cos(theta_bounds[1]), r_outer * np.cos(theta_bounds[1])]
    radial_y = [r_inner * np.sin(theta_bounds[1]), r_outer * np.sin(theta_bounds[1])]
    ax.plot(radial_x, radial_y, [BC[0], BC[1]], color='black', linewidth = 5)

    # Infinitesimals
    dr = (r_outer - r_inner) / (hs[0] - 1)
    dtheta = 2 * np.pi / (hs[1] - 1)
    if eq == 'Heat':
        dt = dr**2 / (100 * c**2)
    else:
        dt = dr / (100 * c)

    # Case separation for titles and labels
    if eq == "Heat":
        ax.set_title('Temperature development in a circular domain', fontsize=18, fontname='STIXGeneral')
        surf = ax.plot_surface(X, Y, Z_init, alpha=0.7, cmap=plt.cm.jet, vmin=zmin, vmax=zmax)
        cbar = fig.colorbar(surf)
        cbar.ax.set_ylabel('Temperature in ' + r'$^\circ\mathrm{C}$', rotation=270, fontsize=14, labelpad=20)
        ax.set_axis_off()

    if eq == "Wave":
        ax.set_title('Simulation of waves in a circular domain', fontsize=18, fontname='STIXGeneral')
        surf = ax.plot_surface(X, Y, Z_init, alpha=0.7, cmap="magma", vmin=zmin, vmax=zmax)
        cbar = fig.colorbar(surf)
        cbar.ax.set_ylabel('Wave height in meters', rotation=270, fontsize=14, labelpad=20)
        ax.set_axis_off()

    # Finite Difference Method
    if eq == 'Wave':
        start_col = 0
        if theta_bounds[0] !=0 or theta_bounds[1] != 2*np.pi:
                start_col = 1
        for row in range(1, hs[0] - 1):
            for col in range(start_col, hs[1]-start_col):
                Z_0[row, col] = Z_0[row, col] - 2 * laplacian_polar(Z_dot_init, r, row, col, dr, dtheta) + \
                                1/2 * c**2 * dt**2 * laplacian_polar(Z_0, r, row, col, dr, dtheta)

        zs.append(Z_0.copy())
        zs.append(Z_init.copy())

        for iteration in range(2, frames):
            surf.remove()
            Z_temp = np.zeros((hs[0], hs[1]))
            Z_temp[:,0] = BC[0]  # Inner radius
            Z_temp[:,-1] = BC[1]  # Outer radius

            for row in range(1, hs[0]-1):
                for col in range(1, hs[1]):
                    Z_temp[row, col] += 2 * zs[iteration - 1][row, col] - zs[iteration - 2][row, col] + \
                                        c**2 * dt**2 * laplacian_polar(zs[iteration - 1], r, row, col, dr, dtheta)

            zs.append(Z_temp.copy())
            zmin, zmax = Z_temp.min(), Z_temp.max()
            surf = ax.plot_surface(X, Y, zs[iteration - 2], alpha=1, cmap="magma", vmin=zmin, vmax=zmax)
            cbar.update_normal(surf)
            plt.pause(dt)

    if eq == 'Heat':
        start_col = 0
        if theta_bounds[0] !=0 or theta_bounds[1] != 2*np.pi:
                start_col = 1
        for iteration in range(frames):
            surf.remove()
            
            for row in range(1, hs[0] - 1):
                for col in range(start_col, hs[1]-start_col):
                    Z_init[row, col] = Z_init[row, col] + c**2 * dt * laplacian_polar(Z_init, r, row, col, dr, dtheta)
            zmin, zmax = Z_init.min(), Z_init.max()
            surf = ax.plot_surface(X, Y, Z_init, alpha=1, cmap=plt.cm.jet, vmin=zmin, vmax=zmax)
            cbar.update_normal(surf)
            plt.pause(dt)

simulation_pdes_polar(r_bounds=[5, 10], theta_bounds = [0,np.pi/3], hs=[60, 60], BC=[0,0], BC_THETA = [0,0], c=10, frames=20000, eq='Heat')
