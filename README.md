s100py
======
[![Build Status](https://travis-ci.com/noaa-ocs-s100/s100py.svg?branch=master)](https://travis-ci.com/noaa-ocs-s100/s100py)

Python API and Utilities for Working with IHO S-100 HDF5 Data Formats

Overview
--------

This python package provides utilities for encoding hydrographic
datasets in the International Hydrographic Organization (IHO) S-100
HDF5 format.

Background
----------

The IHO S-100 standard is a data framework for digital products and
services for hydrographic, maritime, and GIS communities, comprised of
multiple data encoding formats designed for interoperability with
Electronic Navigational Charts (ENCs).

The initial focus of this package is on three of the S-100 encoding
formats:

-   S-102 Bathymetric Surface 
-   S-104 Water Level Information for Surface Navigation
-   S-111 Water Currents for Surface Navigation

However, support for additional formats will likely be added in the
future.

For further information about S-100 formats, see the [IHO
website](http://s100.iho.int/S100/).

Create S100 Products
--------------------

- Follow [S-102](https://s100py.readthedocs.io/en/latest/s102.html#) examples to create an S-102 File
- Follow [S-104](https://s100py.readthedocs.io/en/latest/s104.html#) examples to create an S-104 File
- Follow [S-111](https://s100py.readthedocs.io/en/latest/s111.html#) examples to create an S-111 File

For S100 API Developers
-----------------------

- [Extending the API](https://s100py.readthedocs.io/en/latest/extending_the_api.html)
- [Example API](https://s100py.readthedocs.io/en/latest/sample_api.html)
- [Using Example API](https://s100py.readthedocs.io/en/latest/using_sample_api.html)
- [Advanced Example API Usage](https://s100py.readthedocs.io/en/latest/more_sample_api.html)

Requirements
------------

This codebase is written for Python 3 and relies on the following python
packages:

-   h5py
-   numpy
-   gdal


Installation
------------

This package requires the GDAL Python bindings be present, so it usually can\'t 
just be installed using `pip install gdal`. We recommend installing GDAL 
either through a package manager (e.g. `conda`, `apt`, `yum`, `pacman`) 
or by compiling from scratch. [Miniconda](https://docs.conda.io/en/latest/miniconda.html) 
is probably the easiest method.

Once `gdal` has been installed, s100py can be installed using `pip`:

```bash
pip install s100py
```

Release Notes
-------------
**Version 1.0.0-rc.1 (2021-02-11)**

*Note: All minor releases have been deprecated*

Authors
-------

-   Barry Gallagher (NOAA), <barry.gallagher@noaa.gov>
-   Erin Nagel (UCAR), <erin.nagel@noaa.gov>
-   Glen Rice (NOAA), <glen.rice@noaa.gov>
-   Jason Greenlaw


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
