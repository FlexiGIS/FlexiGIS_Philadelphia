import os

import numpy as np
from scipy.interpolate import griddata
import matplotlib
import matplotlib.pyplot as plt 

def compute_optimum(load, wind, solar):

    # Calculate the covariance matrix for the three timeseries
    combdata = np.row_stack((load, wind, solar))
    covmat = np.cov(combdata)

    # Calculate all weights of wind and solar we wish to model, and the corresponding grid share.
    step = 0.001
    a = np.linspace(0, 1, round(1/step))
    b = np.linspace(0, 1, round(1/step))
    aa, bb = np.meshgrid(a, b)
    cc = 1-aa-bb

    # Grid share can only be positive or zero, so downselect
    a = aa[cc >= 0]
    b = bb[cc >= 0]
    c = cc[cc >= 0]

    # Build a 3d matrix of the coefficients to multiply by the covariance matrix based on load-a*wind-b*solar:
    coeffmat = np.zeros([3, 3, len(a)])
    coeffmat[0][0][:] = 1
    coeffmat[1][1][:] = a**2
    coeffmat[2][2][:] = b**2
    coeffmat[0][1][:] = coeffmat[1][0][:] = -a
    coeffmat[0][2][:] = coeffmat[2][0][:] = -b
    coeffmat[1][2][:] = coeffmat[2][1][:] = -a*-b

    # Each tiled 3x3 matrix represents the covariance matrix for load,a*wind,b*solar
    mult = np.tile(covmat[:, :, None], len(a)) * coeffmat
    # Sum the scaled covariance matrix to compute the standard deviation of combination
    stdev = np.sqrt(np.sum(mult, axis=(0, 1)))

    # # Alternate Method that uses variances explicitly
    # var_l = covmat[0][0]
    # var_w = covmat[1][1]
    # var_s = covmat[2][2]
    # covar_lw = covmat[0][1]
    # covar_ls = covmat[0][2]
    # covar_ws = covmat[1][2]
    # stdev = np.sqrt(var_l + a**2*var_w + b**2*var_s - 2*a*covar_lw - 2*b*covar_ls + 2*a*b*covar_ws)

    # Combine data for output
    outdata = np.column_stack((stdev, a, b, c))

    # Some additional statistics
    sd_min = np.min(stdev)
    sd_min_a = a[np.argmin(stdev)]
    sd_min_b = b[np.argmin(stdev)]
    sd_min_c = c[np.argmin(stdev)]
    rl_min = sd_min_a * wind + sd_min_b * solar + sd_min_c - load

    print("Minimum: a = ", sd_min_a, "b = ", sd_min_b, "c = ", sd_min_c, "stdev = ", sd_min)
    print("balancing energy", np.sum(rl_min[rl_min < 0]))

    return outdata


def genplot(datfile):
    # Extract data from file
    data=np.genfromtxt(datfile, delimiter=',')
    S = data[:, 0]
    a = data[:, 1]
    b = data[:, 2]
    c = data[:, 3]

    # Compute the wind fraction among the renewables
    AB = a / (a + b)
    ab = np.nan_to_num(AB)

    # Subtract the baseline load variability
    S = S - S[np.bitwise_and(a==0,b==0)]

    # Generate gridded data for the plot axes
    grid_c, grid_ab = np.mgrid[0:1:100j, 0:1:100j]
    grid_ab1 = griddata((c, ab), S, (grid_c, grid_ab), method='nearest')
    clevels = np.linspace(0,1.2,13) 

    # Create a figure, plot the data
    plt.figure(figsize=(6,4.5))
    plt.imshow(grid_ab1.T, extent=[0,1,0,1], origin='lower',aspect="auto", cmap="jet", vmin=0, vmax=1.2)
    cbar = plt.colorbar()
    CS = plt.contour(grid_c, grid_ab,grid_ab1, clevels, colors='white',linewidths=1.5)

    # Apply tick formatting
    matplotlib.rcParams['xtick.direction'] = 'out'
    matplotlib.rcParams['ytick.direction'] = 'out'
    ax = plt.gca()
    ax.tick_params(axis='both', which='major', labelsize=10)
    cbar.set_label('Required Balancing Energy [a.u.]', rotation=90,fontsize=11)
    cbar.ax.tick_params(axis='both', which='major', labelsize=10)
    plt.ylim(0,1)

    # Labels
    plt.title("(Configurations of renewable power generation and the required balancing energy from Strorage)")
    plt.xlabel('Imported Energy',fontsize=12)
    plt.ylabel('Solar vs. Wind',fontsize=12)
    plt.clabel(CS, clevels[::2], inline=1, fontsize=8, fmt='%0.2f', colors='white')

    plt.show()


if __name__ == "__main__":
    # Define the input and output filenames
    loadfile = 'data/norm_load.csv'
    solarfile = 'data/norm_solar.csv'
    windfile = 'data/norm_wind.csv'
    outfile = 'output/standard_deviation.csv'

    # Read the input files
    l = np.loadtxt(loadfile,skiprows=1,usecols=1,delimiter=',')
    s = np.loadtxt(solarfile,skiprows=1,usecols=1,delimiter=',')
    w = np.loadtxt(windfile,skiprows=1,usecols=1,delimiter=',')

    alldata = compute_optimum(l, w, s)

    # Create the output directory if necessary
    dirname = os.path.dirname(outfile)
    if not os.path.exists(dirname):
        os.makedirs(dirname)

    # Save the output file
    np.savetxt(outfile, np.asarray(alldata), delimiter=",")

    # Plot the data
    genplot(outfile)
