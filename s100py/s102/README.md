s102
======

Python Utilities and API for Working with IHO S-102 Data Formats

Overview
--------

This python package provides an api and utilities for encoding  S-102
hydrographic datasets in the International Hydrographic Organization (IHO) 
S-100 format.


Build an Executable
-------------------

create a nomkl environment (mkl takes hundreds of megs compressed -- far more than the program warrants)
copy the s100py folder into the lib/site-packages
then use pyinstaller to make an executable

`conda create -n nomkl python=3.6 nomkl pyinstaller gdal h5py numpy`

Copy 
`cp s100py/ envs/nomkl/lib/site-packages/s100py`
change directory
`cd envs/nomkl/lib/site-packages/s100py`

`pyinstaller --onefile -n make_s102 utils.py `
(if proj.db isn't being found then add the data and modify the utils.py script to find it)
(to debug if the files are included use --onedir instead of --onefile)
`pyinstaller --onefile --add-data C:\PydroXL_19_Dev\envs\nomkl2\Library\share\proj;Library\share\proj -n make_s102 utils.py`

The resulting executable will be located at `dist/make_s102.exe`

Usage
-----
Odds are that the only function you will need is from_gdal(). 
You call it with the input path (some geocoded raster, like a BAG)
and the output path (for an HDF5 file in S102 format). You can override 
some of the metadata elements with a dictionary as well.

If you are converting a set of arrays or a format not supported by GDAL 
then use the from_arrays_with_metadata() method. You supply it numpy arrays
and information about the horizontal and vertical coordinate systems and it does the rest.

Create a blank file and see the root object is empty

```python
from s100py import s102
f = s102.S102File("c:\\temp\\test.s102.h5")
print(f.root)
<class 's100py.s102.S102Root'>
OrderedDict()
```

What properties are allowed in the root? Use the get_standard_properties method which will include expected HDF5 attributes, subgroups and dataset names (if any):

```python
f.root.get_standard_properties()
['bathymetry_coverage', 'east_bound_longitude', 'epoch', 'extent_type_code', ..., 'west_bound_longitude']
```

Knowing the spec, east_bound_longitude is a floating point attribute while bathymetry_coverage is a subgroup which will eventually hold the depth and uncertainty values. How could you tell that without reading the specs? Each piece of data, like east_bound_longitude, will have an associated python property with the same name followed by “_type”.:

```python
f.root.east_bound_longitude_type
float
```
There is also a method for creating a default value, again with the properties name followed by “_create”. This is useful if there is a valid default or the type is another S100 class. It is also used by the api to automatically create multiple data properties at once.:
```python
f.root.east_bound_longitude_create()  # makes a default value - which in this case isn't very useful
f.root.east_bound_longitude  # let's see the default
0.0
```
You can always just set the value for any data. Since a default for east_bound_longitude would never make much sense we could have jsut set the value right away like this:
```python
f.root.east_bound_longitude = -124.5
```
As mentioned before, the data is stored in HDF5 (although other formats could easily be added such as XML). The names in HDF5 are determined from the S100 specs which don’t adhere to Python’s pep8 standards so there is a translation that occurs. Similar to the _type and _create there is a class attribute ending in “_attribute_name” which let’s us know what the name would be in HDF5, if that is important to you.

What are the east_bound_longitude and bathymetry_coverage names in S102 nomenclature
```python
f.root.east_bound_longitude_attribute_name
'eastBoundLongitude'
f.root.bathymetry_coverage_attribute_name
'BathymetryCoverage'
```
So now we will dive into the bathymetry_coverage, what is its value?
```python
f.root.bathymetry_coverage
Traceback (most recent call last):
  File "C:\PydroTrunk\Miniconda36\NOAA\site-packages\Python3\s100py\s102.py", line 1042, in bathymetry_coverage
    return self._attributes[self.bathymetry_coverage_attribute_name]
KeyError: 'BathymetryCoverage'
```
Well, you have to create it first! Since we didn’t use bathymetry_coverage_create() or initialize_properties() yet then we need to make the data. bathymetry_coverage is a BathymetryContainer but you don’t really want to look that up in the api, so let the create method do it for you. Saying f.root.bathymetry_coverage is also a bit long, so let’s make a reference shortcut called “bathy” (of course this wouldn’t work if bathymetry_coverage was just a float and not an S100 data class).:
```python
f.root.bathymetry_coverage_create()
bathy = f.root.bathymetry_coverage  # grab a reference to the data
```
Next is a little more complex, S100 says you can have multiple items and they would be named Name_NNN where _NNN is a zero padded number. The BathmetryCoverage inside the BathymetryCoverage (yes, they duplicated the names) is one of these, so in HDF5 it’s going to be BathymetryCoverage/BathymetryCoverage.001. Note the dot in the name – because S102 also is different than S100 which uses an underscore. You can see this gets encoded as a ‘BathymetryContainer’ which is really a “List” type.:
```python
bathy.bathymetry_coverage_type  # see what type this was (though we don't really need to)
s100py.s102.api.BathymetryContainer
>>> bathy.bathymetry_coverage_attribute_name  # We don't need to know this either
'BathymetryCoverage[\\._]\\d+'
```
So, what we need to do is _create() the list then populate it with a new item. If you let the api do the work, append_new_item() will make an instance of the right class for you. We’ll call ours “bathy_01” (FYI, it’ll be named BathymetryCoverage.001 in the HDF5 file):
```python
bathy.bathymetry_coverage_create()
bathy_list = bathy.bathymetry_coverage
bathy_01 = bathy_list.append_new_item()
```
So, what can go in bathy_01? Let’s check the standard_properties so we can put in values into its properties:
```python
bathy_01.get_standard_properties()
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
```

Ugh, thats a lot. Let’s just make north_bound_latitude.

```python
bathy_01.north_bound_latitude = 10.5
bathy_01
<class 's100py.s102.api.BathymetryFeatureInstance'>
OrderedDict([('northBoundLatitude', 10.5)])
```
You can see the latitude above and we made a reference named bathy_01, but you can also access the data using list notation (it is derived from a list really)
```python
bathy_list[0]  # look, here's our data!
<class 's100py.s102.BathymetryFeatureInstance'>
OrderedDict([('northBoundLatitude', 10.5)])
```
Finally let’s initialize everything in this bathy coverage, NOTE it overwrites our north latitude, so we should have initialized first.

```python
bathy_01.initialize_properties()
print(bathy_01)
<class 's100py.s102.BathymetryFeatureInstance'>
OrderedDict([('northBoundLatitude', 0.0), ('Group[\\._]\\d+', []), ('eastBoundLongitude', 0.0),
('extentTypeCode', False), ('gridOriginLatitude', 0.0), ('gridOriginLongitude', 0.0), ...])
```

Do you see that weird ‘Group[._]\d+’ – I happen to know that is another list (named Group.01, Group.02 etc) inside the BathymetryCoverage.001. That Group object is actually where the data grids for depth and uncertainty would go. But the point here is that a method to find those lists exists too, get_standard_list_properties() will tell you the HDF5 naming and the python name:
```python
bathy_01.get_standard_list_properties()
{'Group[\\._]\\d+': 'bathymetry_group'}
```
