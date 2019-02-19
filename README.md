----------------------------------------------------------
Created by Stefanie Falk - February 2019
----------------------------------------------------------
A collection of my past and present analsysis scripts written in python2.7.

Directory structure:
o General purpose:
  - mytools
    > color_maps.py
    > line_styles.py
    > met_tools.pyx     |
    > root_tools.pyx    |*
    > netcdf_tools.py
  | Need to be compiled using cython!
  | 1. Remove __init__*
  | 2. Insert file which shall be compiled in setup.py
  | 3. > python setup.py build_ext --inplace
  | 4. touch __init__.py
  |* Not relevant. Does not work without pre-installed ROOT (CERN software).
  
o Basic examples and test cases:
  - examples
    > bar_test.py
    > cartopy_example_read_openifs.py   |
    > cartopy_example_subplots.py
    > cartopy_example_tick_labels.py
    > fitting.py
    > hist2d_with_labels.py
    > test_fft.py
  | For starters: This is probably the most relevant script,
  | since it shows how to plot ECWMF OpenIFS data using
  | cartopy map projections.
o Special but limited purpose:
  - shipping
    > plot_shipping_emissions.py
    > plot_shipping_fc.pc
    > EMEP_shipping_emission_data.txt   | Shipping data in bad ASCII format
  - regridding
    > generate_cdo_grid_describtion.py
    > prepare_grid.sh
    > hardacre_grid.txt                 |
    > osloctm3_grid.txt                 |
  | Generated files to be used in cdo regrid precedure.
  | 1. > python generate_cdo_grid_describtion.py
  | 2. > ./prepare_grid.sh
o Analysis scrips used in publications
  - BrXplo                              |
  - OsloCTM3dd                          |
  | I will not list all the content here.
o Current working directoy
  - ozone_metrics
  
