import argparse
import matplotlib.pylab as plt
from matplotlib import ticker
import numpy as np
import datetime
import re
import sys


try:
    import netCDF4 as nc
except ImportError:
    print("Please install the Python interface to netCDF")
    sys.exit()

def plot_var(axes, x, y, label):
    # sets the maximum number of labels on the x axis to be 10
    xticks = ticker.MaxNLocator(10)

    var = data.variables[label]
    y = var[:].squeeze()

    long_label = '{} ({})'.format(var.long_name, var.units)
    y_label = re.sub(r'(\s\S*?)\s', r'\1\n',long_label)

    axes.plot(x, y)
    axes.set_ylabel(y_label)
    axes.xaxis.set_major_locator(xticks)

parser = argparse.ArgumentParser()
parser.add_argument('-p', '--model-path', type=str, required=True,
                    help='Path to FABM0D-ERSEM output')
args, _ = parser.parse_known_args()
model_path = args.model_path

data = nc.Dataset(model_path, 'r')

times = data.variables['time']
dates = nc.num2date(times[:],
                    units=times.units,
                    calendar=times.calendar)
dates = [str(d).split(" ")[0] for d in dates]

input_data = ["light_parEIR", "temp", "salt"]
fig, ax_arr = plt.subplots(len(input_data), 1)
DPI = fig.get_dpi()
fig.set_size_inches(1200.0/float(DPI),800.0/float(DPI))
for ax, label in zip(ax_arr, input_data):
    plot_var(ax, dates, data, label)
fig.suptitle("Input data")

input_data = ["N1_p", "N3_n", "B1_c", "P2_c"]
fig, ax_arr = plt.subplots(len(input_data), 1)
DPI = fig.get_dpi()
fig.set_size_inches(1200.0/float(DPI),800.0/float(DPI))
for ax, label in zip(ax_arr, input_data):
    plot_var(ax, dates, data, label)
fig.suptitle("Output data")
plt.show()