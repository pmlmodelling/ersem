import argparse
from matplotlib import pyplot as plt
import numpy as np
import os
import sys


parser = argparse.ArgumentParser()
parser.add_argument('-p', '--model-path', type=str, required=True,
                    help='Path to GOTM-FABM-ERSEM output')
args, _ = parser.parse_known_args()
model_path = args.model_path

try:
    import netCDF4 as nc
except ImportError:
    print("Please install the Python interface to netCDF")
    sys.exit()

def plot_time_series(axes, data, depth, model_var_name):
    times = data.variables['time']
    dates = nc.num2date(times[:],
                        units=times.units,
                        calendar=times.calendar)

    z = data.variables['z'][:].squeeze()
    zi = data.variables['zi'][:].squeeze()
    var = data.variables[model_var_name]

    # Interpolate variable data to the given depth below the moving free surface
    var_time_series = []
    for i in range(var.shape[0]):
        depth_offset = depth + zi[i, -1]  # Remove offset introduced by the moving free surface

        var_time_series.append(np.interp(depth_offset, z[i, :], var[i, :].squeeze()))
    dates = [str(d).split(" ")[0] for d in dates]
    axes.plot(dates, var_time_series)
    axes.set_ylabel('{} ({})'.format(model_var_name, var.units))
    for i, t in enumerate(axes.get_xticklabels()):
        if (i % np.floor(len(dates)/10)) != 0:
            t.set_visible(False)
    axes.xaxis.set_ticks_position('none')


model_file_name = "L4_time_daily_mean_16.06.nc"

fname = os.path.join(model_path, model_file_name)
data = nc.Dataset(fname, 'r')
data.set_auto_maskandscale(True)

# Depth at which to compare the model and data, given as relative to the moving free surface.
depth = 0.0

fig, ax_arr = plt.subplots(3,1)
plot_time_series(ax_arr[0], data, depth, 'N1_p')

# Nitrate and nitrite
plot_time_series(ax_arr[1], data, depth, 'N3_n')

# Silicate
plot_time_series(ax_arr[2], data, depth, 'N5_s')

# Plot formatting
for ax in ax_arr:
    ax.tick_params(axis='both', which='major')
    ax.tick_params(axis='both', which='minor')
    ax.set_ylim(0)

plt.show()
