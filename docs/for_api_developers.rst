For S100 API Developers
=======================

.. toctree::
    :maxdepth: 4

    extending_the_api
    example_api
    using_example_api
    advanced_example_api_usage

S100 is a family of specifications where each data product has it's own extensions to the basic format dictated
by S100.  S100 files are HDF5 where the data's types, names and structure is in the aforementioned specs.
s100py uses pep8 naming for the api itself.

The :any:`S1XXFile` class handles the top level data management and is derived from the h5py.File class.
There is a "root" data member for the file object which is where the S100/HDF5 data is held.
Each data specification should create it's own S100 type file object but there is a basic s100.S100File class
which can read the general top level metadata that pertains to all data specifications.

The three primary data types in HDF5 are attributes, datasets and groups.  s100py encapsulates groups
(which can have attributes) with the :any:`S1xxObject` class.
Datasets are derived from either :any:`S1xxDatasetBase` or :any:`S1xxGridsBase`.
S100 also adds groups that have a trailing number and can have an arbitrary number of occurrences.
These objects are managed with the :any:`S1xxCollection` class.
Classes derived from :any:`S1xxDatasetBase` and  :any:`S1xxCollection` will act as
lists of groups (:any:`S1xxObject`)

Each S100 class has data and some methods to determine what data is expected to be contained therein.
The function get_standard_properties() will show what child data is referenced in the specs.
The :meth:`~s100py.s1xx.S1xxObject.initialize_properties` method will create default values for all expected child data.
Using initialize_properties(True) will recurse the data spec and create a skeleton of all children and their children.
The get_standard_properties_mapping() show the HDF5 names and they pythonic pep8 names that are used with them.

Inside of each S100 data class are data values and functions related to each data value.
For a "dataname" value there are two properties, __dataname_type__ and __dataname_hdf_name__,
and a function dataname_create.

Note, with S100 v5.0 the product specs S102, S104, S111 starting stricter types.
Specifying a numpy type in the __dataname_type__ property will write the data as that specific type in the HDF5.
The python type hints can still be "float" but the "__dataname_type__" would be "numpy.float32".

When the File is being written any data that has not been set or initialized will be omitted.
This way optional data will not appear in the HDF5 file at all.