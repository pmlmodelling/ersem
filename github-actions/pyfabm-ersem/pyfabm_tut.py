import matplotlib.pylab as plt
import numpy as np
import os
import pyfabm
import yaml
import tempfile

def main():
    """
    Run pyfabm-ERSEM tutorial

    :param model_path: Full path to netCDF model output
    :type model_path: str

    :return: list of oxygen saturation concentrations
    :rtype: list
    """
    dir_path = os.path.dirname(os.path.realpath(__file__))
    ersem_dir = os.path.dirname(os.path.dirname(dir_path))

    # Path to ERSEM yaml file
    ersem_yaml_file = os.path.join(ersem_dir,
                                   'testcases',
                                   'fabm-ersem-26.02-dvm.yaml')

    if not os.path.isfile(ersem_yaml_file):
        raise RuntimeError("Could not find Ersem yaml file with the "
                           "{}".format(ersem_yaml_file))


    # Removing DVM from pyfabm tests
    with open(ersem_yaml_file) as f:
        raw = f.read()

    data = yaml.safe_load(raw.replace('\t', '  '))
    del data["instances"]["dm"]
    del data["instances"]["dw"]
    del data["instances"]["db"]

    with tempfile.NamedTemporaryFile(mode="w", suffix=".yaml", delete=False) as tmp:
        yaml.safe_dump(data, tmp, sort_keys=False)
        temp_path = tmp.name

    # Create model
    model = pyfabm.Model(temp_path)

    # Configure the environment
    model.findDependency('longitude').value = -4.15
    model.findDependency('latitude').value = 50.25
    model.findDependency('number_of_days_since_start_of_the_year').value = 0.
    model.findDependency('temperature').value = 10.
    model.findDependency('wind_speed').value = 1.
    model.findDependency('surface_downwelling_shortwave_flux').value = 50.
    model.findDependency('practical_salinity').value = 35.
    model.findDependency('pressure').value = 10.
    model.findDependency('density').value = 1035.
    model.findDependency('mole_fraction_of_carbon_dioxide_in_air').value = 280.
    model.findDependency('absorption_of_silt').value = 0.07
    model.findDependency('bottom_stress').value = 0.
    model.findDependency('cell_thickness').value = 1.
    model.setCellThickness(1)

    # Verify the model is ready to be used
    assert model.checkReady(), 'One or more model dependencies have not been fulfilled.'

    # Define ranges over which temperature and salinity will be varied
    n_points = 50
    temperature_array = np.linspace(5, 35, n_points)
    salinity_array = np.linspace(25, 45, n_points)

    # Create an array in which to store oxygen_saturation_concentrations
    oxygen_saturation_concentration = np.empty((n_points, n_points), dtype=float)

    # Calculate oxygen saturation concentrations using ERSEM
    for t_idx, t in enumerate(temperature_array):
        model.findDependency('temperature').value = t
        for s_idx, s in enumerate(salinity_array):
            model.findDependency('practical_salinity').value = s
            _ = model.getRates()
            oxygen_saturation_concentration[t_idx, s_idx] = \
                    model.findDiagnosticVariable('O2/osat').value

    # Create the figure
    figure = plt.figure()
    axes = plt.gca()

    # Set color map
    cmap = plt.colormaps.get_cmap('YlOrRd')

    # Plot
    plot = axes.pcolormesh(salinity_array,
                           temperature_array,
                           oxygen_saturation_concentration,
                           shading='auto',
                           cmap=cmap)

    axes.set_xlabel('Salinity (psu)')
    axes.set_ylabel('Temperature (deg. C)')



    # Add colour bar
    cbar = figure.colorbar(plot)
    cbar.set_label('O$_{2}$ (mmol m$^{-3}$)')

    plt.show()

    return oxygen_saturation_concentration

if __name__ == "__main__":
    main()
