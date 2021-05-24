import argparse
import matplotlib.pyplot as plt
from matplotlib import ticker
import numpy as np
import os
import sys
import re

try:
    import netCDF4 as nc
except ImportError:
    print("Please install the Python interface to netCDF")
    sys.exit()


def main(model_path):
    """
    Run GOTM-ERSEM tutorial

    :param model_path: Full path to netCDF model output
    :type model_path: str

    :return: Dictionary containing output model variables for systests
    :rtype: dict
    """

    # Dictionary used to save output variables for testing
    output_vars = {}
    def plot_time_series(axes, data, depth, model_var_name):
        # sets the maximum number of labels on the x axis to be 10
        xticks = ticker.MaxNLocator(10)

        times = data.variables['time']
        dates = nc.num2date(times[:],
                            units=times.units,
                            calendar=times.calendar)

        z = data.variables['z'][:].squeeze()
        zi = data.variables['zi'][:].squeeze()
        var = data.variables[model_var_name]

        long_label = '{} ({})'.format(var.long_name, var.units)
        y_label = re.sub(r'(\s\S*?)\s', r'\1\n',long_label)

        # Interpolate variable data to the given depth below the moving free surface
        var_time_series = []
        for i in range(var.shape[0]):
            # Remove offset introduced by the moving free surface
            depth_offset = depth + zi[i, -1]

            var_time_series.append(np.interp(depth_offset, z[i, :], var[i, :].squeeze()))

        dates = [str(d).split(" ")[0] for d in dates]
        axes.plot(dates, var_time_series)
        axes.set_ylabel(y_label)
        axes.xaxis.set_major_locator(xticks)
        output_vars[model_var_name] = var_time_series
        output_vars['dates'] = dates

    data = nc.Dataset(model_path, 'r')
    data.set_auto_maskandscale(True)

    # Depth at which to compare the model and data, given as relative to the moving free surface.
    depth = 0.0

    fig, ax_arr = plt.subplots(3,1)
    DPI = fig.get_dpi()
    fig.set_size_inches(1200.0/float(DPI),800.0/float(DPI))

    names = ['N1_p', 'N3_n', 'N5_s']
    for ax, name in zip(ax_arr, names):
        plot_time_series(ax, data, depth, name)

        ax.tick_params(axis='both', which='major')
        ax.tick_params(axis='both', which='minor')
        ax.set_ylim(0)

    plt.show()

    return output_vars

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--model-path', type=str, required=True,
                        help='Path to GOTM-FABM-ERSEM output')
    args, _ = parser.parse_known_args()
    model_path = args.model_path

    main(model_path)
