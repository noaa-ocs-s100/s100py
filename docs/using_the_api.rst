Using the API
=============

S100 is a family of specifications where each data product has it's own extensions to the basic format dictated
by S100.  S100 files are HDF5 where the data's types, names and structure is in the aforementioned specs.
s100py uses pep8 naming for the api itself.

The :any:`S1XXFile` class handles the top level data management and is derived from the h5py.File class.
There is a "root" data member for the file object which is where the S100/HDF5 data is held.
Each data specification should create it's own S100 type file object but there is a basic s100.S100File class
which can read the general top level metadata that pertains to all data specifications.

The three primary data types in HDF5 are attributes, datasets and groups.  s100py encapsulates groups
(which can have attributes) with the :any:`S1xxAttributesBase` class.
Datasets are derived from either :any:`S1xxDatasetBase` or :any:`S1xxGridsBase`.
S100 also adds groups that have a trailing number and can have an arbitrary number of occurrences.
These objects are managed with the :any:`S1xxMetadataListBase` class.
Classes derived from :any:`S1xxDatasetBase` and  :any:`S1xxMetadataListBase` will act as
lists of groups (:any:`S1xxAttributesBase`)

Each S100 class has data and some methods to determine what data is expected to be contained therein.
The function get_standard_properties() will show what child data is referenced in the specs.
The :meth:`~s100py.s1xx.S1xxAttributesBase.initialize_properties` method will create default values for all expected child data.
Using initialize_properties(True) will recurse the data spec and create a skeleton of all children and their children.
The get_standard_properties_mapping() show the HDF5 names and they pythonic pep8 names that are used with them.

Inside of each S100 data class are data values and functions related to each data value.
For a "dataname" value there are two properties, __dataname_type__ and __dataname_attribute_name__,
and a function dataname_create.

When the File is being written any data that has not been set or initialized will be omitted.
This way optional data will not appear in the HDF5 file at all.

Code Sample for an S102 file
----------------------------

Create a blank file and see the root object is empty ::

    >>> from s100py import s102
    >>> f = s102.S102File("c:\\temp\\test.s102.h5")
    >>> print(f.root)
    <class 's100py.s102.S102Root'>
    OrderedDict()

What properties are allowed in the root?  Use the :any:`get_standard_properties` method which will include expected
HDF5 attributes, subgroups and dataset names (if any)::

    >>> f.root.get_standard_properties()
    ['bathymetry_coverage', 'east_bound_longitude', 'epoch', 'extent_type_code', ..., 'west_bound_longitude']

Knowing the spec, east_bound_longitude is a floating point attribute while bathymetry_coverage is a subgroup which
will eventually hold the depth and uncertainty values.  How could you tell that without reading the specs?
Each piece of data, like east_bound_longitude, will have an associated python property with the same
name followed by "_type".::

    >>> f.root.__east_bound_longitude_type__
    float

There is also a method for creating a default value, again with the properties name followed by "_create".
This is useful if there is a valid default or the type is another S100 class.
It is also used by the api to automatically create multiple data properties at once.::

    >>> f.root.east_bound_longitude_create()  # makes a default value - which in this case isn't very useful
    >>> f.root.east_bound_longitude  # let's see the default
    0.0

You can always just set the value for any data.  Since a default for east_bound_longitude would never make much sense
we could have jsut set the value right away like this::

    >>> f.root.east_bound_longitude = -124.5

As mentioned before, the data is stored in HDF5 (although other formats could easily be added such as XML).
The names in HDF5 are determined from the S100 specs which don't adhere to Python's pep8 standards so there is
a translation that occurs.  Similar to the _type and _create there is a class attribute ending in "_attribute_name"
which let's us know what the name would be in HDF5, if that is important to you.

What are the east_bound_longiture and bathymetry_ceverage names in S102 nomenclature ::

    >>> f.root.__east_bound_longitude_attribute_name__
    'eastBoundLongitude'
    >>> f.root.__bathymetry_coverage_attribute_name__
    'BathymetryCoverage'

So now we will dive into the bathymetry_coverage, what is its value? ::

    >>> f.root.bathymetry_coverage
    Traceback (most recent call last):
      File "C:\PydroTrunk\Miniconda36\NOAA\site-packages\Python3\s100py\s102.py", line 1042, in bathymetry_coverage
        return self._attributes[self.__bathymetry_coverage_attribute_name__]
    KeyError: 'BathymetryCoverage'

Well, you have to create it first!  Since we didn't use bathymetry_coverage_create() or
:meth:`~s100py.s1xx.S1xxAttributesBase.initialize_properties`  yet then we need to make the data.
bathymetry_coverage is a :any:`BathymetryContainer` but you don't really want to look that up in the api,
so let the create method do it for you.  Saying f.root.bathymetry_coverage is also a bit long, so let's make a reference
shortcut called "bathy"
(of course this wouldn't work if bathymetry_coverage was just a float and not an S100 data class).::

    >>> f.root.bathymetry_coverage_create()
    >>> bathy = f.root.bathymetry_coverage  # grab a reference to the data

Next is a little more complex, S100 says you can have multiple items and they would be named Name_NNN where _NNN
is a zero padded number.  The BathmetryCoverage inside the BathymetryCoverage (yes, they duplicated the names)
is one of these, so in HDF5 it's going to be BathymetryCoverage/BathymetryCoverage.001.  Note the dot in the name --
because S102 also is different than S100 which uses an underscore.
You can see this gets encoded as a 'BathymetryContainer' which is really a "List" type.::

    >>> bathy.__bathymetry_coverage_type__  # see what type this was (though we don't really need to)
    s100py.s102.api.BathymetryContainer
    >>> bathy.__bathymetry_coverage_attribute_name__  # We don't need to know this either
    'BathymetryCoverage[\\._]\\d+'

So, what we need to do is _create() the list then populate it with a new item.
If you let the api do the work, :meth:`~s100py.s1xx.S1xxMetadataListBase.append_new_item` will make an instance of
the right class for you.  We'll call ours "bathy_01" (FYI, it'll be named BathymetryCoverage.001 in the HDF5 file)::

    >>> bathy.bathymetry_coverage_create()
    >>> bathy_list = bathy.bathymetry_coverage
    >>> bathy_01 = bathy_list.append_new_item()

So, what can go in bathy_01?  Let's check the standard_properties so we can put in values into its properties::

    >>> bathy_01.get_standard_properties()
    ['bathymetry_group', 'date_time_of_first_record', 'date_time_of_last_record',
    'east_bound_longitude', 'west_bound_longitude', 'north_bound_latitude', 'south_bound_latitude',
    'extent_type_code',
    'grid_origin_latitude', 'grid_origin_longitude', 'grid_origin_vertical',
    'grid_spacing_latitudinal', 'grid_spacing_longitudinal', 'grid_spacing_vertical',
    'instance_chunking',
    'num_grp', 'num_points_latitudinal', 'num_points_longitudinal', 'num_points_vertical',
    'number_of_times', 'start_sequence', 'time_record_interval',
    'vertical_extent_maximum_z', 'vertical_extent_minimum_z',
     ]

Ugh, thats a lot.  Let's just make north_bound_latitude. ::

    >>> bathy_01.north_bound_latitude = 10.5
    >>> bathy_01
    <class 's100py.s102.api.BathymetryFeatureInstance'>
    OrderedDict([('northBoundLatitude', 10.5)])

You can see the latitude above and we made a reference named bathy_01,
but you can also access the data using list notation (it is derived from a list really) ::

    >>> bathy_list[0]  # look, here's our data!
    <class 's100py.s102.BathymetryFeatureInstance'>
    OrderedDict([('northBoundLatitude', 10.5)])

Finally let's initialize everything in this bathy coverage,  NOTE it overwrites our north latitude,
so we should have initialized first. ::

    >>> bathy_01.initialize_properties()
    >>> print(bathy_01)
    <class 's100py.s102.BathymetryFeatureInstance'>
    OrderedDict([('northBoundLatitude', 0.0), ('Group[\\._]\\d+', []), ('eastBoundLongitude', 0.0),
    ('extentTypeCode', False), ('gridOriginLatitude', 0.0), ('gridOriginLongitude', 0.0), ...])

Do you see that weird 'Group[\\._]\\d+' -- I happen to know that is another list (named Group.01, Group.02 etc)
inside the BathymetryCoverage.001.  That Group object is actually where the data grids for depth and uncertainty
would go.  But the point here is that a method to find those lists exists too,
:meth:`~s100py.s1xx.S1xxAttributesBase.get_standard_list_properties` will tell you the HDF5 naming and the python name::

    bathy_01.get_standard_list_properties()
    {'Group[\\._]\\d+': 'bathymetry_group'}

