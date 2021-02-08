Using The Example API
=====================

To use the artificial sample api in :any:`example_api` we can follow the steps below.
The code below is also in the "s100py/tests" folder in the "sample_api_test.py" file which should be runnable
assuming you have a compliant, working python environment.

Open a new file in the local directory ::

    write_to_file = S999File("test.sample.h5")

Create our first, basic object a "MyObject" and give it a string value and add it to the file ::

    the_first_object = MyObject()
    the_first_object.data_value = "A sample string"
    write_to_file.root.my_first_object = the_first_object

We could have set the attribute via the constructor as well ::

    my_other_first_object = MyObject(data_value="A sample string")

Create a MyLocation object and add it to the file using the root object ::

    a_location = MyLocation()
    # create values for all child attributes and by passing in True it would recurse all grandchildren and beyond too
    a_location.initialize_properties(True)
    a_location.utm_zone = 18
    a_location.east_bound_longitude = 33.5

Let's tack on a non-S999 attribute and make up one called "extraData".  This shouldn't exist per the spec
but for extensibility it is supported. ::

    a_location.add_metadata("extraData", 12345)

The api supplies some introspection functionality too.  You can use a python IDE to see the function names
but some alternatives are available to see the python data names are and
what S100 specifies the names in the file and what type of data it should be.

get_standard_properties() gives the python names of attributes based on the S100+ specs ::

    print("Show what attributes are held inside our location class, can use either the instance or the class name itself")
    print(a_location.get_standard_properties())
    print(MyLocation.get_standard_properties())
    # both would return --
    # ['east_bound_longitude', 'extent_type_code', 'north_bound_latitude', 'south_bound_latitude', 'utm_zone', 'west_bound_longitude']

get_standard_keys() returns a list of the S100+ spec names as will be written to the HDF5 file ::

    # get_standard_keys will return the names expected from an HDF5 file based on the S100 specs
    print(a_location.get_standard_keys())
    # returns --
    # ['eastBoundLongitude',
    #  'extentTypeCode',
    #  'northBoundLatitude',
    #  'southBoundLatitude',
    #  'utmZone',
    #  'westBoundLongitude']

get_standard_properties_mapping() returns a dictionary of S100+ spec name (seen in HDF5) to python name,
basically tying together the two previous functions.
The mapping to see the HDF5 naming only works on an instance - not the class unfortunately. ::

    print(a_location.get_standard_properties_mapping())
    # returns --
    # {'eastBoundLongitude': 'east_bound_longitude',
    #  'extentTypeCode': 'extent_type_code',
    #  'northBoundLatitude': 'north_bound_latitude',
    #  'southBoundLatitude': 'south_bound_latitude',
    #  'utmZone': 'utm_zone',
    #  'westBoundLongitude': 'west_bound_longitude'}

get_all_keys will show all the HDF5 data names held in the object -- this can include non-standard data
If there is no additional data then it will be the same as get_standard_keys()::

    print(a_location.get_all_keys())
    # returns -- (notice the extraData in the list.
    # {'eastBoundLongitude',
    #  'extentTypeCode',
    #  'extraData',
    #  'northBoundLatitude',
    #  'southBoundLatitude',
    #  'utmZone',
    #  'westBoundLongitude'}

We can use the standard python 'del' command to remove a piece of data we don't want ::

    # get rid of west -- it shouldn't show up in the HDF5 file after this
    del a_location.west_bound_longitude

Now we will add the location to the file::

    write_to_file.root.my_location_group = a_location

Moving on the the "data_group" which will be a 'list' of data with _001 after it's name in HDF5.
In python we can just treat it as a standard list and not worry about the implementation detail of naming.
We will use the _create() helper to make the list. ::

    write_to_file.root.data_group_create()  # this makes the DataGroups which is a list container for the DataGroupObject

We could also have said root.data_group = DataGroups().

Let's pretend we are in a prompt and are trying to make this without looking at the docs.
Use the get_standard_properties() we saw above and we see the 'data_group' from the root.
However when we get the standard properties of the data_group it returns an empty list.
This is because it is a special list class and it holds a list of one type of data which
can be seen by checking the metadata_type(). ::

    # Introspect the data group to figure out what it wants without reading the docs :)
    print(write_to_file.root.get_standard_properties())
    # ['data_group', 'dataset_with_names', 'my_first_object', 'my_location_group']
    print(write_to_file.root.data_group.get_standard_properties())
    # []
    print(type(write_to_file.root.data_group.metadata_type()))
    # <class 'sample_api_test.DataGroupObject'>

Once we realize the data_group is actually a list and what it holds, then we can see what would be inside that data. ::

    print(type(write_to_file.root.data_group.metadata_type()).get_standard_properties())
    # ['data_grid', 'name_of_data']

So let's make some DataGroupObjects.  The potentially tricky part is that the data_grid will be a HDF5.dataset
while the name_of_data was set to be an enumeration.

Actually both are treated like any other string or numeric attribute.
The data_grid (dataset) just needs to be a numpy array or hdf5 dataset and it will work.
The enumeration can be set using either the strings or numbers that the S100+ spec describes.

Being a list we can make an arbitrary amount of them.  Let's plan on three so we can make the enumerations
in different ways and have different shaped datasets. ::

    data_1 = DataGroupObject()
    data_1.name_of_data = "spam"
    data_1.data_grid = numpy.zeros([2, 5])
    write_to_file.root.data_group.append(data_1)
    del data_1.name_of_data

    data_2 = DataGroupObject()
    data_2.name_of_data = 2
    data_2.data_grid = numpy.ones([3, 4])
    write_to_file.root.data_group.append(data_2)

    data_3 = write_to_file.root.data_group.append_new_item()
    data_3.data_grid = numpy.arange(0, 10, .75)
    data_3.name_of_data = data_3.__name_of_data_type__["cheese"]

Ok, let's make a fourth element for the list.  There is a append_new_item() which creates and returns the proper
dataype.  We'll give it a name by passing in an enumeration value but no grid data.  All on one line no less!! ::

    write_to_file.root.data_group.append_new_item().name_of_data = MONTY(2)

Compound datasets are lists of objects that will eventually be held as HDF5 datasets.
They work as lists of attribute classes in s100py so you never need to know they are datasets in reality.

Let's make our first entry and let it take on default values using the initialize_properties(). ::

    attr_1 = datasetWithNames()
    attr_1.initialize_properties(True)

Next we'll make a second set of attributes ::

    attr_2 = datasetWithNames()
    attr_2.initialize_properties(True)
    attr_2.attr_str = "A custom string this time"
    attr_2.attr_int = 27
    attr_2.attr_float = 35.0

And we'll add it to the file by creating a DatasetWithNames_List and passing the attr_1, attr_2 to it's constructor ::

    write_to_file.root.dataset_with_names = DatasetWithNames_List((attr_1, attr_2))  # also could have used _create and append/append_new_item

Now we'll save it all to disk ::

    write_to_file.write()
    write_to_file.close()

Now let's open that data file and spot check against the values we thought we wrote::

    read_from_file = S999File("test.sample.h5", "r")
    assert read_from_file.root.dataset_with_names[1].attr_int == 27
    assert read_from_file.root.dataset_with_names[0].attr_str in b"used a default string"
    assert write_to_file.root.data_group[1].name_of_data == MONTY(2)
    assert write_to_file.root.data_group[2].data_grid[1] == 0.75  # the second element of the range
    assert write_to_file.root.my_location_group.east_bound_longitude == 33.5
    try:
        write_to_file.root.my_location_group.west_bound_longitude
    except:
        pass  # all good, the value shouldn't exist
    else:
        assert write_to_file.root.my_location_group.west_bound_longitude is None  # this should not exist, even as None

Finally, let's make a copy of the data on disk and modify one of the values.
You can use the HDFView app to confirm everything worked. ::

    copy_of_file = S999File("test.rewrite.h5")
    copy_of_file.root = read_from_file.root
    # this shows how to initialize on creation
    copy_of_file.root.my_location_group = MyLocation(utm_zone=22, east_bound_longitude=11, extra_attr="This shouldn't even be here, but it works")
    copy_of_file.write()

