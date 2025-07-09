s100py
======
[![Tests](https://github.com/noaa-ocs-s100/s100py/actions/workflows/s100py-tests.yaml/badge.svg?event=push)](https://github.com/noaa-ocs-s100/s100py/actions/workflows/s100py-tests.yaml)

API and Utilities for Working with IHO S-100 HDF5 Data Formats

Overview
--------

This python package includes an API and utilities for encoding hydrographic
datasets in the International Hydrographic Organization (IHO) S-100
HDF5 format.

See [s100py Read the Docs](https://s100py.readthedocs.io/en/latest/) for more information.

Background
----------

The IHO S-100 standard is a data framework for digital products and
services for hydrographic, maritime, and GIS communities, comprised of
multiple data encoding formats designed for interoperability with
Electronic Navigational Charts (ENCs).

This package includes the following S-100 encoding formats:

-   S-102 Bathymetric Surface 
-   S-104 Water Level Information for Surface Navigation
-   S-111 Water Currents for Surface Navigation

However, support for additional formats as they become available will
be added in the future.  For further information about S-100 formats,
see the [IHO website](http://s100.iho.int/S100/).

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

This codebase is written for Python >=3.8 and relies on the following python
packages:

-   h5py
-   numpy
-   gdal


Installation
------------

Conda
-----
Installing s100py from the conda-forge channel can be achieved by adding conda-forge to your channels with:
```bash
conda config --add channels conda-forge
conda config --set channel_priority strict
```
Once the conda-forge channel has been enabled, s100py can be installed with conda:
```bash
conda install s100py
```
or with mamba:
```bash
mamba install s100py
```
It is possible to list all the versions of s100py available on your platform with conda:
```bash
conda search s100py --channel conda-forge
```
or with mamba:
```bash
mamba search s100py --channel conda-forge
```
Alternatively, mamba repoquery may provide more information:

**Search all versions available on your platform:**

```bash
mamba repoquery search s100py --channel conda-forge
```
**List packages depending on `s100py`:**
```bash
mamba repoquery whoneeds s100py --channel conda-forge
```
**List dependencies of `s100py`:**

```bash
mamba repoquery depends s100py --channel conda-forge
```

PIP
---

If installing with PIP, this package requires GDAL Python bindings be present,
so it usually can't be installed just using `pip install gdal`, but for more information
visit https://pypi.org/project/GDAL/.

Once `gdal` has been installed, s100py can be installed using `pip`:

```bash
pip install s100py
```

Release Notes
-------------
**Version 2.0.0 (2025-07-09)**
- This is a major release, which includes several new APIs that encapsulate the latest S-100 standard (S-100 Edition 5.2.0) and S-100 Product Specifications (S-102 Edition 3.0.0, S-111 Edition 2.0.0, S-104 Edition 2.0.0)
- Additional convenience utilities have been added to generate:
  * S-102 Edition 3.0.0 datasets
    * **Not Supported in s100py v2.0.0**:
      * Option for multiple VerticalDatums with multiple BathymetryCoverages
  * S-111 Edition 2.0.0 datasets
    * **Breaking Changes in s100py v2.0.0**:
      * Must call `add_surface_current_instance()` before calling `add_data_from_arrays()` to generate datasets
      * `update_metadata()` parameter requirements have changed and only dataDynamicity and dateTimeOfFirstRecord are required in the metadata dictionary parameter
    * **New Features in s100py v2.0.0**:
      * Optional speed_uncertainty and direction_uncertainty grids now supported, if adding an uncertainty grids, speed_uncertainty and direction_uncertainty must be set to True in `create_s111()` and uncertainty data grids must be added to `add_data_from_arrays()`
      * `add_surface_current_instance()` can be used to create multiple surface current feature instances in a single dataset
    * **Not Supported in s100py v2.0.0**:
      * Utilities for generating S-111 Data Coding Formats 4 and 8
  * S-104 Edition 2.0.0 datasets
    * **Breaking Changes in s100py v2.0.0**:
      * Must call `add_water_level_instance()` before calling `add_data_from_arrays()` to generate datasets
      * `update_metadata()` parameter requirements have changed and only dataDynamicity and dateTimeOfFirstRecord are required in the metadata dictionary parameter
    * **New Features in s100py v2.0.0**:
      * Optional water level height uncertainty grid now supported, if adding an uncertainty grid , uncertainty must be set to True in `create_s104()` and uncertainty data grids must be added to `add_data_from_arrays()`
      * `add_water_level_instance()` can be used to create multiple water level feature instances in a single dataset
    *    * **Not Supported in s100py v2.0.0**:
      * Option for multiple VerticalDatums with multiple coverages
- S-102 Edition 3.0.0 (S-100 Edition 5.2.0) significant changes include:
    * Make Uncertainty optional in BathymetryCoverage
    * Renamed QualityOfSurvey to QualityOfBathymetryCoverage
    * Added option for multiple VerticalDatums with multiple BathymetryCoverages
    * dataOffsetCode now required in s102 and must have value of 5 (barycentric)
    * BathymetryCoverage now is dataCodingFormat is 2, regular grid
    * QualityOfBathymetryCoverage dataCodingFormat remains 9, feature oriented
    * BathymetryCoverage Container
      * DataCodingFormat revised from 9 to 2 (regularGrid)
      * common_point_rule revised from 1 to 2 (low)
      * data_offset_code = 5  (barycenter)
- S-111 Edition 2.0.0 (S-100 Edition 5.2.0) significant changes include:
    * Added directionUncertainty and speedUncertainty added to the values record as optional attributes
    * New fill value for date-time attribute
    * Added UTM zones and newer WGS84 epochs
    * Added verticalCoordinateBase embedded metadata for S-100 consistency
    * Removed “DateTime” as UoM name for surface current time attribute in Group_F
    * Removed ISO metadata files
    * Added restriction on length of string attributes in metadata
    * The following DataDynamicity classification have been added:
        * 6: observedMinusPredicted - Observation minus astronomical prediction
        * 7: observedMinusAnalysis - Observation minus analysis or hybrid
        * 8: observedMinusHindcast - Observation minus hydrodynamic hindcast
        * 9: observedMinusForecast - Observation minus hydrodynamic forecast
        * 10: forecastMinusPredicted - Hydrodynamic forecast minus astronomical prediction
    * Added optional dataOffsetCode enumeration for Offset of data point in cell data, mandatory if data points are at grid cell centres
- S-104 Edition 2.0.0 (S-100 Edition 5.2.0) significant changes include:
    * Scope reduced to use only the regular grid spatial type
    * The only allowed interpolation type for DCF2 is 1 (nearestneighbor)
    * The following s104_DataDynamicity classifications are not allowed in S-104 Edition 2.0:
        * 4: hydrodynamicHindcast - Observation minus astronomical prediction
        * 6: observedMinusPredicted - Observation minus analysis or hybrid
        * 7: observedMinusAnalysis - Observation minus hydrodynamic hindcast
        * 8: observedMinusHindcast - Observation minus hydrodynamic hindcast
        * 9: observedMinusForecast - Observation minus hydrodynamic forecast
        * 10: forecastMinusPredicted - Hydrodynamic forecast minus astronomical prediction
    * verticalCoordinateBase attribute added the only allowed value is verticalDatum (2), the attribute is now mandatory
    * Instance metadata constraints adjusted for Water Level adjustment compatibility
    * Restricted maximum length of HDF5 string attributes
    * Optional uncertainty attribute added to values record, represents the uncertainty at a particular grid point and may be omitted if the uncertainty is the same at all grid points
    * Added UTM zones and newer WGS84 realizations
    * Removed ISO metadata files
    * Optional Group F waterLevelTime attribute has been removed
    * Optional Group F uncertainty attribute has been added
    * Extended format to include grids with datum jumps (multiple vertical datums)
    * Domain extent polygon added to Feature Instance if and only if the feature covers an area with a different vertical datum from the root group
    * Added optional dataOffsetCode enumeration for offset of data point in cell data, mandatory if data points are at grid cell centres
    * The values seaFloor (47), seaSurface (48), and hydrographicZero (49) from S100 allowable vertical and sounding datums are now not allowed in S-104 Edition 2.0

**Version 1.0.0 (2023-10-27)**
- This is a major release, which includes several new APIs that encapsulate the latest S-100 standard (S-100 Edition 5.0.0)
and S-100 Product Specifications (S-102 Edition 2.2.0, S-111 Edition 1.2.0, S-104 Edition 1.1.0) and adds support for
multiple S-100 Editions and Product Specification versions
- Additional convenience utilities have been added to generate:
  * S-102 Edition 2.1.0 & 2.2.0 datasets
  * S-111 Edition 1.0.0 & 1.2.0 datasets
  * S-104 Edition 1.0.0 & 1.1.0 datasets
  * S-102 Data Coding Format 9 (Feature Oriented Regular Grid) datasets
  * S-104 Data Coding Format 3 (Irregular Grid) datasets
  * S-104 Data Coding Format 7 (TIN) datasets
- S-102 Edition 2.1.0 (S-100 Edition 4.0) significant changes include:
    * The depth orientation has changed from positive up to positive down
- S-102 Edition 2.2.0 (S-100 Edition 5.0) significant changes include:
    * General metadata _horizontalDatumReference_ and _horizontalDatumValue_ have been removed
    * General metadata attributes _horizontalCRS_ and _typeOfHorizontalCRS_ were added
    * Stricter datatypes have been set for all attributes
    * Added the feature attribute table _QualityOfSurvey_ which provides the following additional survey metadata:
      * The categorization of the assessment level of bathymetric data for an area
      * If the least depth of detected features in an area was measured
      * If significant features have or have not been detected in the course of a survey
      * The size of detected bathymetric features in an area
      * Percentage of depth that a feature of such size could be detected
      * If full seafloor coverage has been achieved in the area by hydrographic surveys
      * Flag for bathy coverage nodes populated by interpolation
      * The best estimate of the fixed horizontal or vertical accuracy component for positions, depths, heights, vertical
      distances, and vertical clearances
      * The factor to be applied to the variable component of uncertainty equation to provide the best estimate of the
      variable horizontal or vertically accuracy component for positions, depths, heights, vertical distances, and
      vertical clearances
      * The start date of the period of the hydrographic survey
      * The end date of the period of the hydrographic survey
      * The survey filename or ID
      * The authority which as responsible for the hydrographic survey
      * An estimate of the magnitude of the difference between true and estimated bathymetric depth, after all appropriate
      corrections are made
- S-111 Edition 1.2.0 (S-100 Edition 5.0.0) significant changes include:
  * General metadata _horizontalDatumReference_ and _horizontalDatumValue_ have been removed
  * General metadata attributes _horizontalCRS_ and _typeOfHorizontalCRS_ were added
  * Stricter datatypes have been set for all attributes
  * General metadata attribute _depthTypeIndex_ has been modified to only allow _heightOrDepth_ or _layerAverage_
  * General metadata _datasetDeliveryInterval_ was added, describing the expected time interval between availability of
  successive datasets for time-varying data, formatted as an ISO 8601 duration (e.g. 'PnYnMnDTnHnMnS')
  * The enumeration attribute for the classification of data according the relationship between the time of its collection,
  generation or calculation of generation parameters, in relation the time of publication has been changed from
  _typeOfCurrentData_ to _dataDynamicity_
- S-104 Edition 1.1.0 (S-100 Edition 5.0.0) significant changes include:
  * General metadata _horizontalDatumReference_ and _horizontalDatumValue_ have been removed
  * General metadata attributes _horizontalCRS_ and _typeOfHorizontalCRS_ were added
  * Stricter datatypes have been set for all attributes
  * General metadata _datasetDeliveryInterval_ was added, describing the expected time interval between availability of
  successive datasets for time-varying data, formatted as an ISO 8601 duration (e.g. 'PnYnMnDTnHnMnS')
  * The enumeration attribute for the classification of data according the relationship between the time of its collection,
  generation or calculation of generation parameters, in relation the time of publication has been changed from
  _typeOfWaterLevelData_ to _dataDynamicity_
  * General metadata _trendInterval_ was added, describing the interval in minutes over which trend at a particular time
  is calculated

**Version 1.0.0-rc.1 (2021-02-11)**
- This is a major release, which includes a new API that encapsulates the data specifications to allow
  introspection with Python to determine what data is available or should be and what data types would
  be acceptable. Convenience utilities are available to convert data into S102/S104/S111 so detailed
  knowledge of the S100 specs and APIs is not required in most cases
- Support for s102 (bathymetry) and s104 (water levels) has been added
- A consistent API for S100 data structures was added and is used to encode S102, S104, and S111
- The previous S111 library has been migrated to this general S100 API and therefore any
  code written against the previous s100py library will no longer work
- Examples of using the new S111 API are available [here](https://s100py.readthedocs.io/en/latest/s111.html#example-usage)

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
