More Usage of the  API
==========================

Using the data and api shown in :any:`sample_api` and :any:`using_sample_api` we can look at some advanced usage.

The S100+ files are built on HDF5 and use h5py for access.  You can access data directly and use the api to read it
by passing any api class the appropriate h5py object.

We'll use a h5py.File object in the samples below, create one like this: ::

    h5file = h5py.File(filename)

So, for our first example let's read the "MyObject" that would be found in the HDF under the group named "myFirstObject" ::

    obj_location = "/myFirstObject"
    data = MyObject()
    data.read(h5file[obj_location])
    print(data.data_value)  # prints:  'A sample string' based on the previous examples

We can also edit the data any write it back to that same place, rather than writing an entire file at once. ::

    data.data_value = "A revised simple string"
    data.write(h5file[obj_location])  # overwrite the data in the HDF5 file

We could even put it in the another place in the file.
This wouldn'y be valid S100+ data but you could do it for compatibility with other software for example. ::


    h5file.require_group(obj_location + "_new")  # make a different subgroup called myFirstObject_new
    data.write(h5file[obj_location + "_new"])  # and write it to the file

Similarly we can access the enumerations of the datasets from the DataGroupObject we defined previously.  ::

    data = DataGroupObject()
    data.read(h5file["/dataGroup_003"])
    print(data.name_of_data, data.data_grid)  # should show original data

Now we can change it and write it back out::

    data.name_of_data = MONTY['spam']
    data.data_grid[3] = 99
    data.write(h5file["/dataGroup_003"])  # overwrite the existing data

And make a duplicate record if wherever we want ::

    # duplicate to a new name -- we will need to create a group for it to write into or put it into an existing group
    h5file.create_group("/dataGroup_003_new")
    data.write(h5file["/dataGroup_003_new"])

It works for entire lists too.
Our list was located at the root and unlike the classes above, a list needs to start at the parent node.
The lists have to name themselves (Group_001, Group_002) so need to start at it's parent level::

    data = DataGroups()
    data.read(h5file)  # side note: you can use h5file or h5file["/"] interchangeably

    data[2].data_grid[4] = 35
    data.write(h5file)  # overwrite

We can create a copy of the list, but we can't make a copy at the root where it already exists,
so we have to create a subgroup and tell it to write itself there.
We will end up with "/extra_list/Group_001", "/extra_list/Group_002"... ::

    # duplicate
    h5file.create_group("extra_list")
    data.write(h5file["extra_list"])

Compound arrays work the same way ::

    data = DatasetWithNames_List()
    data.read(h5file["/datasetWithNames"])
    data[0].attr_int = 5
    # the dataset type needs to create a dataset with it's own name, so we pass in the parent.
    data.write(h5file["/"])

    # duplicate
    data[0].attr_str = "duplicated"
    h5file.require_group("/new_compound_array")
    data.write(h5file["/new_compound_array"])

    h5file.close()

def test_change_names_on_new_data(revised_filename):
    """ This plays some games with the attribute_name.  Because the data is held in a dictionary based on the hdf5 names,
    changing the mapping between python name and HDF5 name can have consequences.  """
    h5file = h5py.File(revised_filename)

    # set up a standard object but store it in a non-standard group
    obj_with_standard_name = MyObject()
    obj_with_standard_name.data_value = "standard"
    assert obj_with_standard_name.data_value_attribute_name == "dataValue"
    h5file.require_group("/test_standard_name")
    obj_with_standard_name.write(h5file["/test_standard_name"])

    # change just the instance's name for HDF5, doing this BEFORE adding data works fine
    obj_with_non_standard_name = MyObject()
    obj_with_non_standard_name.data_value_attribute_name = "Change_instance_name"
    obj_with_non_standard_name.data_value = "Testing just the curreent instance"
    h5file.require_group("/test_instance_names")
    obj_with_non_standard_name.write(h5file["/test_instance_names"])

    assert obj_with_non_standard_name.data_value == "Testing just the curreent instance"

    # Change the class definition, which can be easier if ALL the data you ever want to read uses that different naming
    MyObject.data_value_attribute_name = "Change_all_classes"
    changed_class_obj = MyObject()
    changed_class_obj.data_value = "Change_the_class_itself"
    h5file.require_group("/test_class_names")
    changed_class_obj.write(h5file["/test_class_names"])

    assert changed_class_obj.data_value_attribute_name == "Change_all_classes"
    assert obj_with_non_standard_name.data_value_attribute_name == "Change_instance_name"

    # but watch out, existing data will also get the new name (but the one we changed just the instance of will be unaffected).
    # our standard name data will now have data that is orphaned and adding/changing the data via the api will only use the new names
    obj_with_standard_name.data_value = "still standard?"
    h5file.require_group("/test_standard_whoa")
    obj_with_standard_name.write(h5file["/test_standard_whoa"])
    assert obj_with_standard_name.data_value_attribute_name == "Change_all_classes"
    assert obj_with_non_standard_name.data_value == "Testing just the curreent instance"

    h5file.close()

def test_changing_names_on_existing_data(revised_filename):
    """ Change the attribute names in existing data, this requires re-mapping the old data to the new name or deleting the old."""
    h5file = h5py.File(revised_filename)
    obj_location = "/datasetWithNames"  # this is the root of the file
    data = DatasetWithNames_List()
    data.read(h5file[obj_location])

    # change the names and values of exising data
    for index, compund_arr in enumerate(data):
        del compund_arr.attr_int  #delete the old data before we rename
        compund_arr.attr_int_attribute_name = "changed_individual_int"
        compund_arr.attr_int = (index + 5) * 2

    # change the name in all the classes in existence at once.  This could corrupt other data in memory, in theory!
    old_name = datasetWithNames.attr_float_attribute_name
    datasetWithNames.attr_float_attribute_name = "changed_class_float"
    for index, compund_arr in enumerate(data):
        compund_arr.__delattr__(old_name)
        compund_arr.attr_float = (index + 6) * 3

    h5file.require_group("/compound_array_changed_names")
    data.write(h5file["/compound_array_changed_names"])

    h5file.close()

