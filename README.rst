######
s100py
######

Python Utilities for Working with IHO S-100 Data Formats


Overview
========
This python package provides utilities for encoding hydrographic datasets in
the International Hydrographic Organization (IHO) S-100 format.


Background
==========
The IHO S-100 standard is a data framework for digital products and services
for hydrographic, maritime, and GIS communities, comprised of multiple data
encoding formats designed for interoperability with Electronic Navigational
Charts (ENCs).

The initial focus of this package is on two of the S-100 encoding formats:

- S-104 Water Level Information for Surface Navigation
- S-111 Surface Currents

However, support for additional formats will likely be added in the future.

**S-111 Data Types**

    1 - Time series data at one or more fixed stations

    2 - Regularly-gridded data at one or more times

    3 - Ungeorectified gridded data or point set data at one or more times

    4 - Time series data for one moving platform

For further information about S-100 formats, see the
`IHO website <http://s100.iho.int/S100/>`_.

Features
========
- Create and modify S-111 compliant HDF5 files in all four data coding formats:
  1. Time-series at fixed station
  2. Regularly-gridded arrays
  3. Ungeorectified gridded arrays (i.e. irregular grid)
  4. Time series for moving platform
- Chop output into multiple subgrids (i.e. tiles), each written to a distinct
  S-111 file, in order to minimize file sizes
- Create and modify HDF5 S-100/S-111 metadata


Requirements
============
This codebase is written for Python 3 and relies on the following python
packages:

- h5py
- numpy
- `thyme <https://github.com/noaa-ocs-modeling/thyme>`_


Example Usage
=============
To create an S-111 regular grid (Type 2) file from the NOS Chesapeake Bay
Operational Forecast System file valid at 7/9/2019 0000 UTC:

.. code-block:: python

    import datetime
    from s100py import s111
    from thyme.model import roms
    data_coding_format = 2
    file_metadata = s111.S111Metadata(
            "Chesapeake Bay",
            "ROMS",
            data_coding_format,
            "US",
            4.5,
            None,
            "CBOFS")
    native_model_file = roms.ROMSFile('/path/to/cbofs_file.nc')
    model_index_file = roms.ROMSIndexFile('/path/to/cbofs_index_file.nc')
    s111.model_to_s111(
            model_index_file,
            [native_model_file],
            '/path/to/s111_directory',
            datetime.datetime(2019,7,8,0,0),
            file_metadata,
            data_coding_format)

To create an S-111 file containing moving platform (Type 4) real-time
observation data:

.. code-block:: python

    from s100py import s111
    file_metadata = s111.S111Metadata(
            "Gulf of Mexico",
            "Moving Platform",
            2, # 2 = Real-time observation
            "US",
            0, # meters below sea surface (0 = at/near surface)
            "station1234")
    input_data = []
    input_data.append(
        s111.S111TimeSeries(
                longitude, # 1D `numpy.ndarray` containing longitude values
                latitude, # 1D `numpy.ndarray` containing latitude values
                speed, # 1D `numpy.ndarray` containing speed values in knots
                direction, # 1D `numpy.ndarray` containing Direction values in arc-degrees
                datetime_values)) # List containing a `datetime.datetime` for each observation in the series
    s111.time_series_to_s111(
            input_data,
            '/path/to/s111_directory',
            file_metadata,
            4) # 4 = Time series for moving Platform

Authors
=======
- Erin Nagel (UCAR), erin.nagel@noaa.gov
- Jason Greenlaw (ERT), jason.greenlaw@noaa.gov

License
=======
This work, as a whole, is licensed under the BSD 2-Clause License (see
`LICENSE <LICENSE>`_), however it contains major contributions from the U.S.
National Oceanic and Atmospheric Administration (NOAA), 2017 - 2019, which are
individually dedicated to the public domain.

Disclaimer
==========
This repository is a scientific product and is not official communication of
the National Oceanic and Atmospheric Administration, or the United States
Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’
basis and the user assumes responsibility for its use. Any claims against the
Department of Commerce or Department of Commerce bureaus stemming from the use
of this GitHub project will be governed by all applicable Federal law. Any
reference to specific commercial products, processes, or services by service
mark, trademark, manufacturer, or otherwise, does not constitute or imply their
endorsement, recommendation or favoring by the Department of Commerce. The
Department of Commerce seal and logo, or the seal and logo of a DOC bureau,
shall not be used in any manner to imply endorsement of any commercial product
or activity by DOC or the United States Government.

Acknowledgments
===============
This software has been developed by the National Oceanic and Atmospheric
Administration (NOAA)/National Ocean Service (NOS)/Office of Coast Survey
(OCS)/Coast Survey Development Lab (CSDL) for use by the scientific and
oceanographic communities.

CSDL wishes to thank the following entities for their assistance:

- NOAA/NOS/Center for Operational Oceanographic Products and Services (CO-OPS)
- Canadian Hydrographic Service (CHS)
- Teledyne CARIS

