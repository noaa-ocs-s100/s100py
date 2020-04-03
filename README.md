s100py
======
[![Build Status](https://travis-ci.com/noaa-ocs-s100/s100py.svg?branch=master)](https://travis-ci.com/noaa-ocs-s100/s100py)

Python Utilities for Working with IHO S-100 Data Formats

Overview
--------

This python package provides utilities for encoding hydrographic
datasets in the International Hydrographic Organization (IHO) S-100
format.

Background
----------

The IHO S-100 standard is a data framework for digital products and
services for hydrographic, maritime, and GIS communities, comprised of
multiple data encoding formats designed for interoperability with
Electronic Navigational Charts (ENCs).

The initial focus of this package is on three of the S-100 encoding
formats:

-   S-102 Bathymetry
-   S-104 Water Level Information for Surface Navigation
-   S-111 Surface Currents

However, support for additional formats will likely be added in the
future.

For further information about S-100 formats, see the [IHO
website](http://s100.iho.int/S100/).

Features
--------

-   Create and modify S-111 compliant HDF5 files in all four data coding
    formats:


        1.  Time-series at fixed station
        2.  Regularly-gridded arrays
        3.  Ungeorectified gridded arrays (i.e. irregular grid)
        4.  Time series for moving platform

-   Chop output into multiple subgrids (i.e. tiles), each written to a
    distinct S-111 file, to reduce file sizes

-   Create and modify HDF5 S-100/S-111 metadata

Requirements
------------

This codebase is written for Python 3 and relies on the following python
packages:

-   h5py
-   numpy
-   [thyme](https://github.com/noaa-ocs-modeling/thyme)

Installation
------------

This package relies on [thyme](https://github.com/noaa-ocs-modeling/thyme),
which requires the GDAL Python bindings be present, so it usually can\'t 
just be installed using `pip install gdal`. We recommend installing GDAL 
either through a package manager (e.g. `conda`, `apt`, `yum`, `pacman`) 
or by compiling from scratch. [Miniconda](https://docs.conda.io/en/latest/miniconda.html) 
is probably the easiest method.

Once `gdal` has been installed, s100py can be installed using `pip`:

```bash
pip install s100py
```

Example Usage
-------------

**Create an S-111 File (Type 1):**

```python
from s100py import s111

data_coding_format = 1

# meters below sea surface (0 = at/near surface)
current_depth = 0

file_metadata = s111.S111Metadata(
        "Gulf of Mexico",  # region
        "harmonic_current_predictions",  # product description
        3,  # current type code for astronomical prediction
        "US",  # producer code
        "station1234",  # station id
        None)  # model identifier

input_data = [s111.S111TimeSeries(
            longitude,  # 1D `numpy.ndarray` containing longitude values
            latitude,  # 1D `numpy.ndarray` containing latitude values
            speed,  # 1D `numpy.ndarray` containing speed values in knots
            direction,  # 1D `numpy.ndarray` containing Direction values in arc-degrees
            datetime_values)]  # List containing a `datetime.datetime` for each observation in the series

s111.time_series_to_s111(
        input_data,
        '/path/to/s111_directory',
        file_metadata,
        data_coding_format,
        current_depth)
```

**Create an S-111 File (Type 2)**

NOS Chesapeake Bay Operational Forecast System file valid at 7/9/2019
0000 UTC:

```python
import datetime
from s100py import s111
from thyme.model import roms

data_coding_format = 2
target_cellsize = 500 # meters

# meters below sea surface (0 = at/near surface), default = 4.5 m
target_depth = 0

file_metadata = s111.S111Metadata(
        "Chesapeake Bay",  # region
        "ROMS_Hydrodynamic_Model_Forecasts",  # product type description
        6,  # current data type for hydrodynamic forecast
        "US",  # producer code
        None,  # station id
        "CBOFS")  # model identifier

native_model_file = roms.ROMSFile('/path/to/nos.cbofs.fields.f001.20190709.t00z.nc')
model_index_file = roms.ROMSIndexFile('/path/to/create/index_file.nc')

try:
    native_model_file.open()
    model_index_file.open()
    model_index_file.init_nc(native_model_file, target_cellsize, 'my_roms_model', '/path/to/shoreline_shapefile.shp')

finally:
    model_index_file.close()
    native_model_file.close()

s111.model_to_s111(
        model_index_file,
        [native_model_file],
        '/path/to/s111_directory',
        datetime.datetime(2019, 7, 9, 0, 0),
        file_metadata,
        data_coding_format,
        target_depth)
```
**Create an S-111 File (Type 3)**

NOS Chesapeake Bay Operational Forecast System file valid at 7/9/2019
0000 UTC:

```python
import datetime
from s100py import s111
from thyme.model import roms

data_coding_format = 3

# meters below sea surface (0 = at/near surface), default = 4.5 m
target_depth = 0

file_metadata = s111.S111Metadata(
        "Chesapeake Bay",  # region
        "ROMS_Hydrodynamic_Model_Forecasts",  # product type description
        6,  # current data type for hydrodynamic forecast
        "US",  # producer code
        None,  # station id
        "CBOFS")  # model identifier

native_model_file = roms.ROMSFile('/path/to/nos.cbofs.fields.f001.20190709.t00z.nc')

s111.model_to_s111(
        None,
        [native_model_file],
        '/path/to/s111_directory',
        datetime.datetime(2019, 7, 9, 0, 0),
        file_metadata,
        data_coding_format,
        target_depth)

```
**Create an S-111 file (Type 4):**

```python
from s100py import s111

data_coding_format = 4

# meters below sea surface (0 = at/near surface)
current_depth = 15

file_metadata = s111.S111Metadata(
        "Western_N_Pacific_Ocean_Philippine_Sea",  # region
        "argos_lagrangian_drifter_12hr_interpolated",  # product type description
        4,  # current type code for analysis or hybrid method
        "US",  # producer code
        None,  # station id
        None)  # model identifier

input_data = [s111.S111TimeSeries(
            longitude,  # 1D `numpy.ndarray` containing longitude values
            latitude,  # 1D `numpy.ndarray` containing latitude values
            speed,  # 1D `numpy.ndarray` containing speed values in knots
            direction,  # 1D `numpy.ndarray` containing Direction values in arc-degrees
            datetime_values)]  # List containing a `datetime.datetime` for each observation in the series

s111.time_series_to_s111(
        input_data,
        '/path/to/s111_directory',
        file_metadata,
        data_coding_format,
        current_depth)
```

Authors
-------

-   Erin Nagel (UCAR), <erin.nagel@noaa.gov>
-   Jason Greenlaw (ERT), <jason.greenlaw@noaa.gov>
-   Barry Gallagher, <barry.gallagher@noaa.gov>
-   Glen Rice, <glen.rice@noaa.gov>

License
-------

This work, as a whole, falls under Creative Commons Zero (see
[LICENSE](LICENSE)).

Disclaimer
----------

This repository is a scientific product and is not official
communication of the National Oceanic and Atmospheric Administration, or
the United States Department of Commerce. All NOAA GitHub project code
is provided on an 'as is' basis and the user assumes responsibility for
its use. Any claims against the Department of Commerce or Department of
Commerce bureaus stemming from the use of this GitHub project will be
governed by all applicable Federal law. Any reference to specific
commercial products, processes, or services by service mark, trademark,
manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce.
The Department of Commerce seal and logo, or the seal and logo of a DOC
bureau, shall not be used in any manner to imply endorsement of any
commercial product or activity by DOC or the United States Government.

Acknowledgments
---------------

This software has been developed by the National Oceanic and Atmospheric
Administration (NOAA)/National Ocean Service (NOS)/Office of Coast
Survey (OCS)/Coast Survey Development Lab (CSDL) for use by the
scientific and oceanographic communities.

CSDL wishes to thank the following entities for their assistance:

-   NOAA/NOS/Center for Operational Oceanographic Products and Services
    (CO-OPS)
