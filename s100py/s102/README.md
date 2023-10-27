s102
====

Python API and Utilities for Working with IHO S-102 Data Formats

Overview
--------

S-102 is an IHO standard outlining formats for storing and sending bathymetric
data and metadata.

Odds are that the only function you will need is from_gdal(). 
You call it with the input path (some geocoded raster, like a BAG) and the output 
path (for an HDF5 file in S102 format). You can override some of the metadata elements
with a dictionary as well.

If you are converting a set of arrays or a format not supported by GDAL 
then use the from_arrays_with_metadata() method. You supply it numpy arrays
and information about the horizontal and vertical coordinate systems and it 
does the rest.

```python
    from s100py import s102
    s102.from_gdal(input_tif, output_h5)
    s102.to_geotiff(input_h5, output_tif)
    # or
    from s100py.s102 import utils
    utils.from_gdal(input_tif, output_h5)
    utils.to_geotiff(input_h5, output_tif)
```

Example Usage
-------------
If you need more involved interaction with the data then import s102 from s100py.
Methods and classes from both the api.py and utils.py will be available.
Older versions of the S100 specs are also available, for example "from s100py.s102 import v2_1".

Most times if someone is accessing the data they will want the depths, uncertainties, and metadata.
Rather than diving into the entire spec to find the data, a few convenience properties are available.

```python

>>> from s100py import s102
>>> bluetopo = s102.S102File(r'C:\data\Bluetopo\RATs\BlueTopo_BC25M26L_20221102b.2_2.h5')
>>> bluetopo.depth[:5, :5]
array([[1000000.  , 1000000.  , 1000000.  , 1000000.  , 1000000.  ],
       [1000000.  , 1000000.  , 1000000.  , 1000000.  , 1000000.  ],
       [  -1198.53,   -1219.43,   -1243.89,   -1293.43,   -1332.8 ],
       [  -1091.72,   -1114.74,   -1154.59,   -1208.1 ,   -1257.05],
       [  -1077.18,   -1090.34,   -1107.16,   -1165.56,   -1215.53]],
      dtype=float32)
>>> bluetopo.uncertainty[:5, :5]
array([[1.000e+06, 1.000e+06, 1.000e+06, 1.000e+06, 1.000e+06],
       [1.000e+06, 1.000e+06, 1.000e+06, 1.000e+06, 1.000e+06],
       [6.293e+01, 6.398e+01, 6.520e+01, 6.768e+01, 6.965e+01],
       [5.759e+01, 5.874e+01, 6.073e+01, 6.341e+01, 6.586e+01],
       [5.686e+01, 5.752e+01, 5.836e+01, 6.128e+01, 6.378e+01]],
      dtype=float32)
>>> bluetopo.grid_origin_latitude, bluetopo.grid_origin_longitude
(2788956.1599999997, 198254.24)
>>> bluetopo.grid_spacing_latitudinal, bluetopo.grid_spacing_longitudinal
(1352.32, 1228.48)
>>> bluetopo.epsg  # note that if this is -1 then the coordinate reference is user defined using other S100 values 
32610
# NOTE: the values above are convenience properties that return the values from the bathymetry_coverage
# To directly access or create the values you could do the following
>>> coverage = bluetopo.root.bathymetry_coverage.bathymetry_coverage[0]
>>> coverage.bathymetry_group[0].values.depth[:5, :5]
>>> coverage.bathymetry_group[0].values.uncertainty[:5, :5]
>>> coverage.grid_origin_latitude, coverage.grid_origin_longitude
>>> coverage.grid_spacing_latitudinal, coverage.grid_spacing_longitudinal
```

In S102 version 2.2 the Feature Attribute Table was added.  
It holds a n array matching the depth and uncertainty arrays in shape and geolocation.
For each value in the array there is a table to lookup metadata about the depth and uncertainty values.
Using the same bluetopo data opened above we can access the Feature Attribute Table.
```python
# Use the convenience property to get the quality of survey data
>>> bluetopo.quality[:5, :5]
array([[      0,       0,       0,       0,       0],
       [      0,       0,       0,       0,       0],
       [1188902, 1188902, 1188902, 1188902, 1188902],
       [1188902, 1188902, 1188902, 1188902, 1188902],
       [1188902, 1188902, 1188902, 1188902, 1188902]], dtype=uint32)

# Use the convenience property to get the FeatureAttributeDataset, but notice it returns a FeatureAttributeDataset which acts like a list
>>> bluetopo.feature_attribute_table
<class 's100py.s102.v2_2.api.FeatureAttributeDataset'>
OrderedDict()<class 's100py.s102.v2_2.api.FeatureAttributeRecord'>
OrderedDict([('id', 11134), ('dataAssessment', 1), ('featuresDetected.leastDepthOfDetectedFeaturesMeasured', 64), ('featuresDetected.significantFeaturesDetected', 0), ('featuresDetected.sizeOfFeaturesDetected', 1000000.0), ('featureSizeVar', 0.0), ('fullSeafloorCoverageAchieved', 1), ('bathyCoverage', 1), ('zoneOfConfidence.horizontalPositionUncertainty.uncertaintyFixed', 500.0), ('zoneOfConfidence.horizontalPositionUncertainty.uncertaintyVariableFactor', 0.0), ('surveyDateRange.dateStart', '1939-01-01'), ('surveyDateRange.dateEnd', '1939-01-01'), ('sourceSurveyID', 'H06499A'), ('surveyAuthority', 'DOC/NOAA/NOS/OCS -- Office of Coast Survey'), ('bathymetricUncertaintyType', <BATHYMETRIC_UNCERTAINTY_TYPE.unknown: 0>)])<class 's100py.s102.v2_2.api.FeatureAttributeRecord'>
OrderedDict([('id', 20908), ('dataAssessment', 1), ('featuresDetected.leastDepthOfDetectedFeaturesMeasured', 64), ('featuresDetected.significantFeaturesDetected', 0), ('featuresDetected.sizeOfFeaturesDetected', 1000000.0), ('featureSizeVar', 0.0), ('fullSeafloorCoverageAchieved', 1), ('bathyCoverage', 1), ('zoneOfConfidence.horizontalPositionUncertainty.uncertaintyFixed', 50.0), ('zoneOfConfidence.horizontalPositionUncertainty.uncertaintyVariableFactor', 0.0), ('surveyDateRange.dateStart', '1991-01-01'), ('surveyDateRange.dateEnd', '1991-01-01'), ('sourceSurveyID', 'B00265'), ('surveyAuthority', 'DOC/NOAA/NOS/OCS -- Office of Coast Survey'), ('bathymetricUncertaintyType', <BATHYMETRIC_UNCERTAINTY_TYPE.unknown: 0>)])<class 's100py.s102.v2_2.api.FeatureAttributeRecord'>
...
OrderedDict([('id', 1188907), ('dataAssessment', 1), ('featuresDetected.leastDepthOfDetectedFeaturesMeasured', 64), ('featuresDetected.significantFeaturesDetected', 0), ('featuresDetected.sizeOfFeaturesDetected', 1000000.0), ('featureSizeVar', 0.0), ('fullSeafloorCoverageAchieved', 0), ('bathyCoverage', 0), ('zoneOfConfidence.horizontalPositionUncertainty.uncertaintyFixed', 500.0), ('zoneOfConfidence.horizontalPositionUncertainty.uncertaintyVariableFactor', 0.0), ('surveyDateRange.dateStart', '2022-10-28'), ('surveyDateRange.dateEnd', '2022-10-28'), ('sourceSurveyID', 'GMRT_61m_US3TX1ZC'), ('surveyAuthority', 'Global Multi-Resolution Topography Data Synthesis (GMRT)'), ('bathymetricUncertaintyType', <BATHYMETRIC_UNCERTAINTY_TYPE.unknown: 0>)])<class 's100py.s102.v2_2.api.FeatureAttributeRecord'>
OrderedDict([('id', 1188908), ('dataAssessment', 1), ('featuresDetected.leastDepthOfDetectedFeaturesMeasured', 64), ('featuresDetected.significantFeaturesDetected', 0), ('featuresDetected.sizeOfFeaturesDetected', 1000000.0), ('featureSizeVar', 0.0), ('fullSeafloorCoverageAchieved', 0), ('bathyCoverage', 0), ('zoneOfConfidence.horizontalPositionUncertainty.uncertaintyFixed', 500.0), ('zoneOfConfidence.horizontalPositionUncertainty.uncertaintyVariableFactor', 0.0), ('surveyDateRange.dateStart', '2022-10-28'), ('surveyDateRange.dateEnd', '2022-10-28'), ('sourceSurveyID', 'GMRT_61m_US3TX1ZD'), ('surveyAuthority', 'Global Multi-Resolution Topography Data Synthesis (GMRT)'), ('bathymetricUncertaintyType', <BATHYMETRIC_UNCERTAINTY_TYPE.unknown: 0>)])

# Get the feature attribute record for a specific feature using the feature id
>>> d = bluetopo.get_feature_attribute_dict()
>>> d[1188902].date_start, d[1188902].survey_authority
('2022-03-28', 'Global Multi-Resolution Topography Data Synthesis (GMRT)')

# NOTE: similar to the depth, uncertainty and grid conveniences, the feature_attribute_table is a convenience property
#    which could use the code below to get the same results
>>> bluetopo.root.quality_of_survey.quality_of_survey[0].quality_group[0].values[:5, :5]
>>> bluetopo.root.quality_of_survey.feature_attribute_table
```

Deeper Example Usage
--------------------

Create a blank file and see the root object is empty

```python

>>> from s100py import s102
>>> f = s102.S102File("c:\\temp\\test.s102.h5", "w")
>>> print(f.root)
<class 's100py.s102.S102Root'>
OrderedDict()
```

What properties are allowed in the root? Use the get_standard_properties method which will include expected HDF5 attributes, subgroups and dataset names (if any):

```python
>>> f.root.get_standard_properties()
['bathymetry_coverage', 'east_bound_longitude', 'epoch', 'extent_type_code', ..., 'west_bound_longitude']
```

Knowing the spec, east_bound_longitude is a floating point attribute while bathymetry_coverage is a subgroup which will eventually hold the depth and uncertainty values. 
How could you tell that without reading the specs? Each piece of data, like east_bound_longitude, will have an associated python property with the same name followed by “_type__”.:

```python
>>> f.root.__east_bound_longitude_type__
float
```
There is also a method for creating a default value, again with the properties name followed by “_create”. 
This is useful if there is a valid default or the type is another S100 class. It is also used by the api to automatically create multiple data properties at once.:
```python
>>> f.root.east_bound_longitude_create()  # makes a default value - which in this case isn't very useful
>>> f.root.east_bound_longitude  # let's see the default
0.0
```
You can always just set the value for any data. Since a default for east_bound_longitude would never make much sense we could have just set the value right away like this:
```python
>>> f.root.east_bound_longitude = -124.5
```
As mentioned before, the data is stored in HDF5 (although other formats could easily be added such as XML). 
The names in HDF5 are determined from the S100 specs which don’t adhere to Python’s pep8 standards so there is a translation that occurs. Similar to the _type__ and _create there is a class attribute ending in “_hdf_name__” which let’s us know what the name would be in HDF5, if that is important to you.

What are the east_bound_longitude and bathymetry_coverage names in S102 nomenclature
```python
>>> f.root.__east_bound_longitude_hdf_name__
'eastBoundLongitude'
>>> f.root.__bathymetry_coverage_hdf_name__
'BathymetryCoverage'
```
So now we will dive into the bathymetry_coverage, what is its value?
```python
>>> f.root.bathymetry_coverage
Traceback (most recent call last):
  File "C:\PydroTrunk\Miniconda36\NOAA\site-packages\Python3\s100py\s102.py", line 1042, in bathymetry_coverage
    return self._attributes[self.__bathymetry_coverage_hdf_name__]
KeyError: 'BathymetryCoverage'
```
Well, you have to create it first! Since we didn’t use bathymetry_coverage_create() or initialize_properties() yet then we need to make the data. bathymetry_coverage is a BathymetryContainer but you don’t really want to look that up in the api, so let the create method do it for you. Saying f.root.bathymetry_coverage is also a bit long, so let’s make a reference shortcut called “bathy” (of course this wouldn’t work if bathymetry_coverage was just a float and not an S100 data class).:
```python
>>> f.root.bathymetry_coverage_create()
>>> bathy = f.root.bathymetry_coverage  # grab a reference to the data
```
Next is a little more complex, S100 says you can have multiple items and they would be named Name_NNN where _NNN is a zero padded number. The BathmetryCoverage inside the BathymetryCoverage (yes, they duplicated the names) is one of these, so in HDF5 it’s going to be BathymetryCoverage/BathymetryCoverage.001. Note the dot in the name – because S102 also is different than S100 which uses an underscore. You can see this gets encoded as a ‘BathymetryContainer’ which is really a “List” type.:
```python
>>> bathy.__bathymetry_coverage_type__  # see what type this was (though we don't really need to)
s100py.s102.api.BathymetryContainer
>>> bathy.__bathymetry_coverage_hdf_name__  # We don't need to know this either
'BathymetryCoverage[\\._]\\d+'
```
So, what we need to do is _create() the list then populate it with a new item. If you let the api do the work, append_new_item() will make an instance of the right class for you. We’ll call ours “bathy_01” (FYI, it’ll be named BathymetryCoverage.001 in the HDF5 file):
```python
>>> bathy.bathymetry_coverage_create()
>>> bathy_list = bathy.bathymetry_coverage
>>> bathy_01 = bathy_list.append_new_item()
```
So, what can go in bathy_01? Let’s check the standard_properties so we can put in values into its properties:
```python
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
```

Ugh, thats a lot. Let’s just make north_bound_latitude.

```python
>>> bathy_01.north_bound_latitude = 10.5
>>> bathy_01
<class 's100py.s102.api.BathymetryFeatureInstance'>
OrderedDict([('northBoundLatitude', 10.5)])
```
You can see the latitude above and we made a reference named bathy_01, 
but you can also access the data using list notation (it is derived from a list really)
```python
>>> bathy_list[0]  # look, here's our data!
<class 's100py.s102.BathymetryFeatureInstance'>
OrderedDict([('northBoundLatitude', 10.5)])
```
Finally let’s initialize everything in this bathy coverage, NOTE it overwrites our north latitude, 
so we should have initialized first.

```python
>>> bathy_01.initialize_properties()
>>> print(bathy_01)
<class 's100py.s102.BathymetryFeatureInstance'>
OrderedDict([('northBoundLatitude', 0.0), ('Group[\\._]\\d+', []), ('eastBoundLongitude', 0.0),
('extentTypeCode', False), ('gridOriginLatitude', 0.0), ('gridOriginLongitude', 0.0), ...])
```

Do you see that weird ‘Group[._]\d+’ – I happen to know that is another list (named Group_001 for S102, in other S100 specs it could be Group.01, Group.02 etc) 
inside the BathymetryCoverage.001. That Group object is actually where the data grids for depth and uncertainty would go. 
But the point here is that a method to find those lists exists too, get_standard_list_properties() will tell you the HDF5 naming and the python name:
```python
>>> bathy_01.get_standard_list_properties()
{'Group[\\._]\\d+': 'bathymetry_group'}
```

For S-102 Developers
--------------------

- [S-102 Module Documentation](https://s100py.readthedocs.io/en/latest/s102.html#s102-module-docs)
- [S-100 Module Documentation](https://s100py.readthedocs.io/en/latest/s100.html)

Authors
-------

-   Barry Gallagher (NOAA), <barry.gallagher@noaa.gov>
-   Glen Rice (NOAA), <glen.rice@noaa.gov>