Advanced Example API Usage
==========================

Using the data and api shown in :any:`example_api` and :any:`using_example_api` we can look at some advanced usage.
The sample code below is also in the s100py/tests/sample_api_test.py where it should run with pytest.

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
This wouldn't be valid S100+ data but you could do it for compatibility with other software for example. ::


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

Compound arrays similarly, when using them you supply it the parent as it need to call create_dataset. ::

    data = DatasetWithNames_List()
    data.read(h5file["/"])
    data[0].attr_int = 5
    data.write(h5file["/"])

We can put a copy in a new location but it will create a dataset (named "datasetWithNames" in this case)
under the parent location we supply. ::

    data[0].attr_str = "duplicated"
    h5file.require_group("/new_compound_array")
    data.write(h5file["/new_compound_array"])

Now to really abuse the system we can change names of the data but this is dangerous and not recommended.
It would both not adhere to the S100+ specs and also potentially be error prone.

First we'll change one instance of an object.
Remember the MyObject has one string attribute that should be named "dataValue"::

    # set up a standard object but store it in a non-standard group
    obj_with_standard_name = MyObject()
    obj_with_standard_name.data_value = "standard"
    assert obj_with_standard_name.__data_value_hdf_name__ == "dataValue"
    h5file.require_group("/test_standard_name")
    obj_with_standard_name.write(h5file["/test_standard_name"])

Now let's make another copy of MyObject and change __data_value_hdf_name__
which defines the mapping from python name to S100+.
Doing this BEFORE adding data works fine. ::

    obj_with_non_standard_name = MyObject()
    obj_with_non_standard_name.__data_value_hdf_name__ = "Change_instance_name"
    obj_with_non_standard_name.data_value = "Testing just the current instance"
    h5file.require_group("/test_instance_names")
    obj_with_non_standard_name.write(h5file["/test_instance_names"])

If you want to get in trouble then you can change the class variable __data_value_hdf_name__ which will then affect
ALL the future and existing instances of MyObject. ::

    MyObject.__data_value_hdf_name__ = "Change_all_classes"
    changed_class_obj = MyObject()
    changed_class_obj.data_value = "Change_the_class_itself"
    h5file.require_group("/test_class_names")
    changed_class_obj.write(h5file["/test_class_names"])

And here is where the weird stuff happens, the obj_with_standard_name we made just above will also now write into
that new location too.  Our standard name data will now have data that is orphaned and adding/changing
the data via the api will only use the new names.

This will end up having the old data under the old name and the new data under the new name -- definitely not
what someone probably wants.::

    obj_with_standard_name.data_value = "still standard?"
    h5file.require_group("/test_standard_whoa")
    obj_with_standard_name.write(h5file["/test_standard_whoa"])

But, if you need to change some existing data, you can do it.
Changing the attribute names in existing data will require re-mapping the old data to the new name
or deleting the old data.

Here we will change some of the items in the compound array.  It had attr_int, attr_float and attr_str.
First we'll change each instance's atrr_int naming and delete the old data and set new data.::

    data = DatasetWithNames_List()
    data.read(h5file["/"])

    # change the names and values of exising data
    for index, compound_arr in enumerate(data):
        del compound_arr.attr_int  #delete the old data before we rename
        compound_arr.__attr_int_hdf_name__ = "changed_individual_int"
        compound_arr.attr_int = (index + 5) * 2

Then we'll change the attr_float naming for the whole class (and any other existing data in the processes memory)::

    old_name = datasetWithNames.__attr_float_hdf_name__
    datasetWithNames.__attr_float_hdf_name__ = "changed_class_float"
    for index, compound_arr in enumerate(data):
        compound_arr.__delattr__(old_name)
        compound_arr.attr_float = (index + 6) * 3

    h5file.require_group("/compound_array_changed_names")
    data.write(h5file["/compound_array_changed_names"])


