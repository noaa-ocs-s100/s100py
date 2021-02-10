s100py
======

Python API and Utilities for Working with IHO S-100 HDF5 Data Formats

Overview
--------

The IHO S-100 standard is a data framework for digital products and
services for hydrographic, maritime, and GIS communities, comprised of
multiple data encoding formats designed for interoperability with
Electronic Navigational Charts (ENCs).

The initial focus of this package is on three of the S-100 encoding
formats:

-    S-102 Bathymetry
-    S-104 Water Level Information for Surface Navigation
-    S-111 Surface Water Currents

Installation
------------

This package requires the GDAL Python bindings be present, so it usually can\'t
just be installed using :any:`pip install gdal`. We recommend installing GDAL
either through a package manager (e.g. :any:`conda`, :any:`apt`, :any:`yum`, :any:`pacman`)
or by compiling from scratch.  `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_
is probably the easiest method.

Once :any:`gdal` has been installed, s100py can be installed using :any:`pip`: ::

    pip install s100py
::

S-100 Product Modules
---------------------

.. toctree::
    :maxdepth: 4

    s111
    s104
    s102
    s100
    s1xx
    for_api_developers

