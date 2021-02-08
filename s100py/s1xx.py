""" This module contains base classes.  they should not be used directly but only used by code
implementing one of the S100 family of specifications.

See :any:`extending_the_api` for further details about using the classes to create or modify an api.
"""

import collections
from abc import ABC, abstractmethod
from typing import Callable, Iterator, Union, Optional, List, Type
import re
import logging
import inspect
import traceback
import datetime
from enum import Enum

import h5py
# @todo - consider removing the numpy dependence
import numpy

Record = s1xx_sequence = Union[numpy.ndarray, h5py.Dataset]
s1xx_sequence_types = s1xx_sequence.__args__

try:
    h5py_string_dtype = h5py.special_dtype(vlen=str)
except:
    h5py_string_dtype = h5py.string_dtype(encoding='utf-8', length=None)


def is_sub_class(cls, clsinfo):
    """ Python 3.7+ changed the behavior of issubclass to raise an exception if the cls object is not a class.
    So when a function is passed in (numpy.array vs numpy.ndarray) it raises an exception that we don't want. 
    """
    try:
        return issubclass(cls, clsinfo)
    except TypeError as e:
        return False

    
class FixedTimeZones(datetime.tzinfo):
    """Fixed offset in minutes east from UTC."""

    def __init__(self, hours, minutes=0, name=""):
        self.__offset = datetime.timedelta(hours=hours, minutes=minutes)
        self.__name = name

    def utcoffset(self, dt):
        return self.__offset

    def tzname(self, dt):
        return self.__name

    def dst(self, dt):
        return datetime.timedelta(0)


def convert_numpy_strings_to_h5py(vals, names=None):
    """ change numpy arrays with "U" into array using the h5py special string_dtype that translates to utf-8 in the file.

    Parameters
    ----------
    vals
        a list or array
    names
        optional list of names to use with the record array

    Returns
    -------
    A new numpy array which will have h5py special types embedded for the strings

    """
    rec_array = numpy.core.records.fromarrays(vals, names=names)
    new_dtype = []
    for i in range(len(rec_array.dtype)):
        nm = rec_array.dtype.names[i]
        if rec_array.dtype[i].type is numpy.unicode_:
            try:  # h5py <=2.9
                dt = h5py.special_dtype(vlen=str)
            except:  # h5py >=2.10
                dt = h5py.string_dtype(encoding='utf-8', length=None)
        else:
            dt = rec_array.dtype[i].descr[0][1]
        new_dtype.append((nm, dt))
    rec_array_revised = numpy.core.records.fromarrays(vals, dtype=new_dtype)
    return rec_array_revised


class S1xxAttributesBase(ABC):
    """ This class implements a general hdf5 group object that has attributes, dataset or sub-groups.
    Works with S1xxMetadataListBase if the subgroups have multiple occurences (like Group.01, Group.02)
    Works with S1xxDatasetBase for things that are stored like a numpy array (dataset) in hdf5

    __version__ must be overridden.
    To call a base class property use super().property, e.g. super().__version__
    This base class is built from the version 2.0.0 that was eventually published Nov. 2019
    """
    _attr_name_suffix = "_hdf_name__"

    def __init__(self, recursively_create_children=False, **kywrds):
        self._hdf5_path = ""
        self._attributes = collections.OrderedDict()
        if recursively_create_children:
            self.initialize_properties(recursively_create_children)
        for ky, val in kywrds.items():
            if ky in self.get_standard_properties():
                exec("self.{} = val".format(ky))
            else:
                self.add_data(ky, val)
        # self._child_groups = {}

    @property
    @abstractmethod
    def __version__(self) -> int:
        return -1

    def __repr__(self):
        """ When converting to a string for display, shows the class and what data is currently set internally."""
        s = str(self.__class__) + "\n"
        s += str(self._attributes)
        return s

    def __delattr__(self, item):
        """ Delete an attribute from the current data.  Does the conversion from python names to S100+ names,
        If those names are not found, will delete a non-standard attibute.  Will raise an AttributeError if not found.
        """
        # mapping = self.get_standard_properties_mapping()
        if item in self._attributes or item in self.get_standard_properties_mapping():
            del self._attributes[item]
        elif item in self.get_standard_properties():
            del self._attributes[eval("self.__{}_hdf_name__".format(item))]
        else:
            del self.__dict__[item]

    def read_simple_attributes(self, group_object):
        """ Reads the standard simple types (strings, ints, floats, datetimes, enums) from the given group as specified by the class specs.

        Parameters
        ----------
        group_object
            The group (an h5py.File is a group too) to read from.

        Returns
        -------
        None

        """
        logging.debug("Reading attributes" + str(self))
        self._hdf5_path = group_object.name
        expected_items = self.get_standard_properties_mapping()
        # basic attributes -- should be simple types so just set them
        for attr_name in group_object.attrs:
            if attr_name not in expected_items:
                logging.info(" The attr/val: " + attr_name + "/" + str(
                    group_object.attrs[attr_name]) + " was in the group_object but not found in the standard attributes")
                self._attributes[attr_name] = group_object.attrs[attr_name]
            else:
                use_type = self.__getattribute__("__" + expected_items[attr_name] + "_type__")
                if is_sub_class(use_type, Enum):
                    logging.debug(" Enumerated attr/val: " + attr_name + "/" + str(group_object.attrs[attr_name]) + " found and read")
                    self.set_enum_attribute(group_object.attrs[attr_name], attr_name, use_type)
                elif is_sub_class(use_type, (datetime.date, datetime.datetime, datetime.time)):
                    logging.debug(" datetime string: " + attr_name + "/" + str(group_object.attrs[attr_name]) + " found and read")
                    self.set_datetime_attribute(group_object.attrs[attr_name], attr_name, use_type)
                else:
                    logging.debug(" Standard attr/val: " + attr_name + "/" + str(group_object.attrs[attr_name]) + " found and read")
                    setattr(self, expected_items[attr_name], group_object.attrs[attr_name])

    def read(self, group_object):
        """ Given an h5py.File or a h5py group then read the data based on the encoded S100+ spec.

        Parameters
        ----------
        group_object
            The group (an h5py.File is a group too) to read from.

        Returns
        -------
        None

        """

        logging.debug("Reading " + str(self))
        self.read_simple_attributes(group_object)

        expected_items = self.get_standard_properties_mapping()

        # keys are HDF5 groups or datasets
        group_lists = self.get_standard_list_properties()
        all_keys = list(group_object.keys())
        all_keys.sort()
        basic_keys = []  # basic keys will be a list of the s102 group names to be directly imported
        list_type_keys = {}  # list_type_keys will be a dictionary where the key is the attribute name to fill and the value is a list of the S102 group names that would be found
        # separate out the keys that belong to a list of values (BathymetryCoverage.01, BathymetryCoverage.02 etc)
        for data_key in all_keys:
            use_key = None
            for list_key, attr_name in group_lists.items():
                if re.match(list_key, data_key):
                    use_key = attr_name
            if use_key:
                lk = list_type_keys.setdefault(use_key, [])
                lk.append(data_key)
            elif data_key in expected_items:
                basic_keys.append(data_key)
            else:
                logging.warning(
                    data_key + " is an HDF5 group and was in the group_object but not found in the standard attributes, SKIPPING!")

        # now read the basic keys and the keys having list data
        # -- only need to call the read once for things that are lists as the Metadata_List class will find all of its occurrences
        for key in basic_keys:
            # create/clear the data
            logging.debug(" Standard attr/val: " + key + " found and reading")

            # read in the HDF5 attributes etc from the group
            use_type = self.__getattribute__("__" + expected_items[key] + "_type__")
            if is_sub_class(use_type, S1xxAttributesBase):
                self.__getattribute__(expected_items[key] + "_create")()
                data = self.__getattribute__(expected_items[key])
                if is_sub_class(use_type, S1xxWritesOwnGroupBase):  # pass the parent
                    data.read(group_object)
                else:  # pass the exact location for the data
                    data.read(group_object[key])
            else:
                self.__setattr__(expected_items[key], group_object[key])
        for list_type_group in list_type_keys:
            # create/clear the data
            logging.debug(" Standard LIST based attr/val: " + list_type_group + " found and reading")
            self.__getattribute__(list_type_group + "_create")()
            # read in the HDF5 attributes etc from the group
            o = self.__getattribute__(list_type_group)
            o.read(group_object)

    def write_simple_attributes(self, group_object):
        # this is for all the types that can be attributes of a group or dataset in HDF5
        # these simple hdf5 attributes can't have subgroups or datasets

        # iterate through all the key/values in _attributes and write to hdf5
        # if a value is an enum then translate to the correct Enum class
        # if a value is a date, time - convert to character string per S100, section 10C-7 table 10C-1
        # otherwise write as a simple attribute and simple type
        self._hdf5_path = group_object.name

        for key, val in self._attributes.items():
            if isinstance(val, s1xx_sequence_types):
                continue  # skip these types for now
            elif isinstance(val, S1xxWritesOwnGroupBase):
                continue  # skip these types for now
            elif isinstance(val, S1xxAttributesBase):
                continue  # skip these types for now
            elif isinstance(val, (datetime.date, datetime.datetime, datetime.time)):
                logging.debug(key + " datetime: {}", val)
                group_object.attrs[key] = val.isoformat()
            elif isinstance(val, Enum):
                logging.debug(key + " enumeration: " + str(val))
                enum_as_dict = collections.OrderedDict([[item.name, item.value] for item in type(val)])
                int_type = numpy.uint8
                try:  # enum_dtype is added in h5py 2.10
                    enumtype = h5py.enum_dtype(enum_as_dict, int_type)
                except AttributeError:  # special_dtype is for h5py <= 2.9
                    enumtype = h5py.special_dtype(enum=(int_type, enum_as_dict))
                try:
                    group_object.attrs.create(key, val.value, dtype=enumtype)
                except TypeError:  # h5py isn't accepting OrderedDict, convert to dict
                    try:
                        enumtype = h5py.enum_dtype(dict(enum_as_dict), int_type)
                    except AttributeError:
                        enumtype = h5py.special_dtype(enum=(int_type, dict(enum_as_dict)))
                    group_object.attrs.create(key, val.value, dtype=enumtype)

            else:
                logging.debug(key + " simple type: " + str(val))
                group_object.attrs[key] = val

    def write(self, group_object):
        """ write the contained data and all it's children into an HDF5 file using h5py.

        Parameters
        ----------
        group_object
            An h5py.File or an h5py group object

        Returns
        -------
        None

        """
        logging.debug("Writing" + " " + str(self))

        self.write_simple_attributes(group_object)

        # iterate through all the key/values in _attributes and write to hdf5
        # if a value is numpy.array then convert to h5py.Dataset
        # if a value is a S1xxAttributesBase instance then create a group and call it's write function
        # if a value is a S1xxWritesOwnGroupBase instance then let it create the group and tell it to write into the current group_object
        for key, val in self._attributes.items():
            if isinstance(val, s1xx_sequence_types):  # this looks inside the typing.Union to see what arrays should be treated like this
                logging.debug(key + " array: " + str(val.shape))
                # convert any strings to bytes or h5py will fail to write the unicode
                # converted_vals = [v if not isinstance(v, str) else v.encode("utf-8") for v in val]
                converted_vals = [v if not isinstance(v, bytes) else v.decode() for v in val]
                try:
                    del group_object[key]
                except KeyError:
                    pass  # didn't exist, no error
                # We are only supporting single column data, check the dtype that the convert function found
                revised_vals = convert_numpy_strings_to_h5py([converted_vals])
                # Now re-package so it shows up as a single column in hdf5 - not sure why the revised_vals array would have
                # to be accessed with two dimensions otherwise -- like val[0][0] instead of val[0]
                revised_2 = numpy.array(converted_vals, dtype=revised_vals.dtype[0])
                try:
                    new_dataset = group_object.create_dataset(key, data=revised_2)
                except Exception as e:
                    raise e
            elif isinstance(val, S1xxWritesOwnGroupBase):
                # things that either create a dataset and have to combine data into it or make multiple sub groups that the parent can't predict
                logging.debug("{}  S100 object - writing itself now...".format(key))
                val.write(group_object)
            elif isinstance(val, S1xxAttributesBase):
                logging.debug(key + " S100 object - writing itself now...")
                new_group = group_object.require_group(key)
                val.write(new_group)

    def write_as_xml(self, etree_object):
        # basically add a flag to read/write functions, then everywhere a group, dataset or attribute is written either use xml or hdf5
        pass  # this pass makes the warnings in PyCharm go away, must be a pattern of if just NotImplemented treats it as abstract
        raise NotImplementedError("flesh this out if we want an xml representation of S100+ file")

    def get_data(self, key):
        return self._attributes[key]

    def add_data(self, key, value):
        self._attributes[key] = value

    def get_s1xx_attr(self, s1xx_name):
        expected_items = self.get_standard_properties_mapping()
        return self.__getattribute__(expected_items[s1xx_name])

    def set_s1xx_attr(self, s1xx_name, val):
        expected_items = self.get_standard_properties_mapping()
        setattr(self, expected_items[s1xx_name], val)

    def get_write_order(self):
        """ Override this method if the write order of attributes/groups/dataset items is important

        Returns
        -------
        A list of key names if order is important, None otherwise.
        """
        return None

    def get_compound_dtype(self):
        """ Override this method if the dtype of compound dataset items is important

        Returns
        -------
        A list of dtype, None otherwise.
        """
        return None

    def get_all_keys(self):  # this includes custom keys
        """ Gets all the non-standard keys that are contained in this object currently as well as all the
        standard keys (S100/HDF5 style names) that could be added.

        Returns
        -------
        list
        """
        current_keys = set(self._attributes.keys())
        current_keys.update(self.get_standard_keys())
        return current_keys

    # def get_standard_keys(self):  # this is only things that have properties associated
    #     props = [p[0] for p in inspect.getmembers(self.__class__, lambda x: isinstance(x, property))]
    #     implemented_properties = [self.__getattribute__(p) for p in props if p.endswith(self._attr_name_suffix) and p[2:-len(self._attr_name_suffix)] in props]
    #     return implemented_properties

    def get_standard_keys(self):  # this is only things that have properties associated
        """  Returns the S100 (HDF5) names for the things that are listed in the specs for this class.

        Returns
        -------
        list
            The S102 HDF5 group/attribute/dataset names from this object that will be written or read from an HDF5 file.
            e.g. BathymetryCoverage or westBoundLongitude etc.

            For the class "Root":
            ['BathymetryCoverage', 'Group_F', 'TrackingListCoverage']
        """
        return list(self.get_standard_properties_mapping().keys())

    def get_standard_list_properties(self):
        """ Returns a list of properties that are lists (children based on S1xxMetadataListBase).
        Basically a way of finding which items will be named <name>_001, <name>_002 etc

        Returns
        -------
        list
            The property names that will have auto-generated names based on their index in a list.

        """
        s100_to_property = self.get_standard_properties_mapping()
        s100_to_property_for_lists = {}
        for s100_attr, prop in s100_to_property.items():
            if is_sub_class(self.__getattribute__("__" + prop + "_type__"), S1xxMetadataListBase):
                s100_to_property_for_lists[s100_attr] = prop
        return s100_to_property_for_lists

    # @classmethod
    # @todo question - if we make all the __xxx_hdf_name__ as staticmethods or class variables then this could be a classmethod.
    # @todo Do we want to allow flexibility of changing the HDF5 name for an instance?
    # a little testing shows that using class variables and letting the user change it would fail in this function.
    # so if we change the _hdf_name to class variables then the user could change the name on demand - for better or worse
    #  but this would not be a classmethod.
    # if we want this to be a classmethod then all the hdf_name would bee to become staticmethods and become functions - meaning adding ()
    def get_standard_properties_mapping(self):
        """ This function autodetermines the HDF5 or xml names and their associated property names.
        Keys are the s100 (HDF5 spelling) strings and the values are the python style naming.

        Returns
        -------
        dict
            dictionary of xml element names as keys and property names as values.

            For the class "Root":
            {'BathymetryCoverage': 'bathymetry_coverage',
            'Group_F': 'feature_information',
            'TrackingListCoverage': 'tracking_list_coverage'}
         """
        s100_to_property = {}
        for prop in self.get_standard_properties():
            s100_to_property[self.__getattribute__("__" + prop + self._attr_name_suffix)] = prop
        return s100_to_property

    def initialize_properties(self, recursively_create_children=False, overwrite=True):
        """ Calls the create function for all the properties of the class.
        Default values will be created for each attribute that is expected to be contained in this object.

        For example, say a class has simple attributes of ESPG code (int) and locatilty (str) and then a class  made from S1xxAttributesBase
        called "extents" which has east and west inside it.

        Calling initialize_properties(recursively_create_children=False) would result in EPSG=0, locality="" and
        an instance of the "extents" class but NO value (nothing would be written to HDF5) for east, west.

        Calling initialize_properties(recursively_create_children=True)  would result in EPSG=0, locality="" and
        an instance of the "extents" class but with east=0.0 and wesst=0.0 as well.

        Calling initialize_properties(recursively_create_children=True, overwrite=False) with an esiting dataset, say locality="test"
        would result in EPSG=0 being made, locality="test" being retained and an instance of the "extents" class with east=0.0 and wesst=0.0 as well.

        Parameters
        ----------
        recursively_create_children
            True = Create children for any child data that would have other children
            False = Only create data for immediate children of this instance
        overwrite
            True = Overwrite existing data with new default data
            False = Keep existing data if it exists but create new data otherwise

        Returns
        -------
        None

        """
        for prop in self.get_standard_properties():
            if overwrite or not self.__getattribute__(prop):
                exec("self.{}_create()".format(prop))
                o = eval("self.{}".format(prop))
                if recursively_create_children and isinstance(o, S1xxAttributesBase):
                    o.initialize_properties(recursively_create_children, overwrite)

    @classmethod
    def get_standard_properties(cls):
        """  This function autodetermines the properties implemented (which have get/set @properties and _hdf_name associated)

        Returns
        -------
        list
            Names of the properties implemented.

            For eample class "Root" might have (for S102):
            ['bathymetry_coverage', 'feature_information', 'tracking_list_coverage']
        """
        # allow for properties or class attributes (the str check does this)
        props = [p[0] for p in inspect.getmembers(cls, lambda x: isinstance(x, (property, str)))]
        # remove the leading double underscore and the _hdf_name__ suffix
        implemented_properties = [p[2:-len(cls._attr_name_suffix)] for p in props if
                                  p.endswith(cls._attr_name_suffix) and p[2:-len(cls._attr_name_suffix)] in props]
        return implemented_properties

    def set_enum_attribute(self, val, hdf_name, enum_type):
        """ Function to set an attribute that is an enumeration type using either it's string or numeric value
        or enumeration instance.

        Parameters
        ----------
        val
            The value as a string, int or Enum().
        hdf_name
            The S100 name (hdf5 spelling).
        enum_type
            The class of enumeration to use if an instance needs to be created.

        Returns
        -------
        None
        """
        if isinstance(val, str):
            val = enum_type[val]
        if isinstance(val, (int, numpy.integer)):
            val = enum_type(val)
        self._attributes[hdf_name] = val

    def set_datetime_attribute(self, val, hdf_name, date_type):
        """
        A DateTime is a combination of a date and a time type. Character encoding of a
        DateTime must follow ISO 8601:2004  ( :2004 took away partial dates with two digit year or just month/day)
        EXAMPLES
        19850412T101530
        2001-07-17T04:50:00
        2012-11-01T00:44:00+10:30
        2001-07-17
        19850412T101530.44
        19850412T101530.44Z
        19850412T101530.44+10
        19850412T101530.44+1030
        19850412T10:15:30Z

        Parameters
        ----------
        val
        hdf_name
        date_type

        Returns
        -------

        """
        re_date = r"(?P<year>\d{4})[-]?(?P<month>\d{2})[-]?(?P<day>\d{2})"
        re_time = r"(?P<hour>\d{2})[: -]?(?P<minute>\d{2})[:-]?(?P<second>\d{2})(?P<decimal_sec>\.\d+)?"
        re_timezone = r"(?P<tz>(Z|(?P<tz_hr>[+-]\d{2})[:]?(?P<tz_min>\d{2})?))?"
        re_time_with_zone = re_time + re_timezone
        re_full_datetime = re_date + "T?" + re_time_with_zone
        re_date_optional_time = re_date + "T?(" + re_time_with_zone + ")?"
        FixedTimeZones

        def _tz(match_obj):
            if match_obj['tz']:
                if match_obj['tz'] == "Z":
                    z = FixedTimeZones(0)
                else:
                    tz_min = int(match_obj['tz_min']) if match_obj['tz_min'] else 0
                    z = FixedTimeZones(int(match_obj['tz_hr']), tz_min)
            else:
                z = None
            return z

        if isinstance(val, str):
            # turns out the python fromisofomrat only reads the same format it would write using datetime.isoformat()
            # try:  # python 3.7+ has fromisoformat() builtin
            #    val = datetime.datetime.fromisoformat(val)
            # except AttributeError:
            # read as a full datetime first.
            match = re.match(re_full_datetime, val)
            if match:
                decimal_sec = int(float(match['decimal_sec']) * 1000000) if match['decimal_sec'] else 0
                zone = _tz(match)
                val = datetime.datetime(int(match['year']), int(match['month']), int(match['day']),
                                        int(match['hour']), int(match['minute']), int(match['second']),
                                        decimal_sec, tzinfo=zone)
            else:
                if is_sub_class(date_type, datetime.date):
                    match = re.match(re_date, val)
                    if match:
                        val = datetime.date(int(match['year']), int(match['month']), int(match['day']))
                elif is_sub_class(date_type, datetime.time):
                    match = re.match(re_time_with_zone, val)
                    if match:
                        decimal_sec = int(float(match['decimal_sec']) * 1000000) if match['decimal_sec'] else 0
                        zone = _tz(match)
                        val = datetime.time(int(match['hour']), int(match['minute']), int(match['second']),
                                            decimal_sec, tzinfo=zone)

            if not match:
                print("failed to parse date and/or time from '" + val + "' storing as string in ", hdf_name)

        if isinstance(val, (datetime.datetime, datetime.date, datetime.time)):
            if isinstance(val, datetime.datetime):
                if is_sub_class(date_type, datetime.date):
                    val = val.date()
                elif is_sub_class(date_type, datetime.time):
                    tz = val.tzinfo
                    val = val.time()
                    # getting the time this way drops the timezone info, so we stored it and now replace it
                    if tz:
                        val = val.replace(tzinfo=tz)

        self._attributes[hdf_name] = val

    def get_hdf5_from_file(self, file_obj):
        if self._hdf5_path:
            obj = file_obj[self._hdf5_path]
        else:
            obj = None
        return obj


class S1xxWritesOwnGroupBase(S1xxAttributesBase):
    """ Derive things that either create a dataset and have to combine data into it or make multiple sub groups that the parent can't predict
    The S1xxAttributesBase will call the derived class' writer without pre-making group for it.
    i.e. the derived class can specify its own group name or dataset and apply specialized logic as needed.
    """
    pass


class S1xxMetadataListBase(list, S1xxWritesOwnGroupBase):
    """ This class represents arrays (noted in UML as *, 1..*, 0..* etc) which is not really part of HDF5.
    The S100 spec is using a atttribute.NNN to repreent this type of record.
    This class takes the supplied name and type and will make it act like a list in python and read/write the data in HDF5 like S102 wants.

    """
    read_re_pattern = r"[_\.](\d+)$"
    write_format_str = "_%03d"

    def __init__(self, *args, **opts):
        # initialize the list in case data was passed in.
        super().__init__(*args, **opts)  # standard init for lists
        S1xxAttributesBase.__init__(self)  # initialize the s100 class

    @property
    @abstractmethod
    def metadata_name(self) -> str:
        raise NotImplementedError()

    @property
    @abstractmethod
    def metadata_type(self) -> type:
        raise NotImplementedError()

    def append_new_item(self):
        self.append(self.metadata_type())
        return self[-1]

    def read(self, group_object):
        logging.debug("Reading " + str(self))

        # keys are HDF5 groups or datasets
        self._hdf5_path = group_object.name

        all_keys = list(group_object.keys())
        keys_to_process = []
        for data_key in all_keys:
            m = re.match(self.metadata_name + self.read_re_pattern, data_key)
            if m:
                keys_to_process.append([int(m.groups()[0]), data_key])
        # sort in order of the integers encoded in the string.  Don't trust that the zero padding is enough or allow string comparison to sort
        # did this in the regular expression instead
        # key=lambda name: int(name[name.find(".")+1:])
        keys_to_process.sort()
        for index, data_key in keys_to_process:
            obj = self.metadata_type()
            obj.read(group_object[data_key])
            self.append(obj)

    def write(self, group_object):
        # Iterate through the values in the list and write to hdf5
        # Write each item as self.metadata_name + ".%03d" % index
        # They should all be the same type but we aren't sure what type they are.
        # Most likely to be another S100_Attribute type.
        # if a value is numpy.array then convert to h5py.Dataset
        # if a value is a S1xxAttributesBase instance then create a group and call it's write function
        # if a value is a date, time - convert to character string per S100, section 10C-7 table 10C-1
        # otherwise write as a simple attribute and simple type
        self._hdf5_path = group_object.name

        logging.debug("Writing" + " " + str(self))
        # create N new group objects named as metadata_name.NNN
        for index, val in enumerate(self):
            name = self.metadata_name + self.write_format_str % (index + 1)
            if isinstance(val, s1xx_sequence_types):  # this looks inside the typing.Union to see what arrays should be treated like this
                raise NotImplementedError()
            # elif isinstance(val, S1xxMetadataListBase):
            #     raise NotImplementedError("Nested Lists")
            elif isinstance(val, S1xxAttributesBase):  # Attributes, List and Datasets all work the same (for now)
                new_group = group_object.require_group(name)
                val.write(new_group)
            elif isinstance(val, (datetime.date, datetime.datetime, datetime.time)):
                logging.debug(name + " datetime: {}", val)
                # @TODO: figure out how to write datetimes
                raise NotImplementedError("DateTimes not supported yet")
            else:
                logging.debug(name + " simple type: " + str(val))
                group_object[name] = val


class S1xxDatasetBase(list, S1xxWritesOwnGroupBase):
    """ The S102 spec stores some things as attributes that could (or should) be stored as attributes.
    This class reads/writes datasets but stores/accesses them as a list of class instances.
    Data access should then be used as object[index].attribute
    So for the FeatureInformation class that would be feat[0].name = "depth" and feat[1].name = "uncertainty"
    """

    def __init__(self, *args, **opts):
        # initialize the list in case data was passed in.
        super().__init__(*args, **opts)  # standard init for lists
        S1xxAttributesBase.__init__(self)  # initialize the s102 class

    @property
    @abstractmethod
    def metadata_name(self) -> str:
        raise NotImplementedError()

    @property
    @abstractmethod
    def metadata_type(self) -> type:
        raise NotImplementedError()

    def append_new_item(self):
        self.append(self.metadata_type())
        return self[-1]

    def __repr__(self):
        s = S1xxAttributesBase.__repr__(self)
        for data_object in self:
            s += str(data_object)
        return s

    def read(self, group_object_parent):
        group_object = group_object_parent[self.metadata_name]
        self.read_simple_attributes(group_object)  # put any attributes from the dataset object into the overall _attributes
        # @fixme -- I think below is more correct, but should be functionally the same as above
        # S1xxAttributesBase.read(self, group_object)
        # del self._attributes[self.metadata_name]  # remove the dataset since we are going to interpret it as a list instead

        list_length = group_object.shape[0]
        has_extra_dimension = False  # NAVO (and maybe the spec) wrote data as shape (2,1) instead of (2,) so we have to use an extra index - data[0][0]
        if len(group_object.shape) > 1:
            if group_object.shape[1] == 1:
                has_extra_dimension = True
        self.clear()
        # expected_items = self.metadata_type().get_standard_properties_mapping()
        for i in range(list_length):
            current_obj = self.metadata_type()
            for data_name in group_object.dtype.names:
                val = group_object[data_name, i]
                if has_extra_dimension:  # the data was in an array of length one while we want the value.
                    val = val[0]
                # setattr(self, expected_items[data_name], val)
                try:
                    current_obj.set_s1xx_attr(data_name, val)
                except KeyError:  # data not expected per S102 spec in the dataset
                    # store any additional data in the _attribute dictionary with each item in the list
                    current_obj._attributes[data_name] = val
            self.append(current_obj)

    def write(self, group_object):
        # @todo - is there a bug here if some instances are missing attributes leading to a mismatched array?
        """ Write out the dataset using order specified with any extra values as unordered but named at the end.

        Parameters
        ----------
        group_object
            HDF5 object to write into

        Returns
        -------
        HDF5 dataset created during the write method

        """
        # First determine the write order of the keys
        logging.debug("Writing" + " " + str(self))
        dataset = None
        if len(self) > 0:
            val = self[0]
            write_keys = []
            if val.get_write_order():  # @todo I think bathycoverage and trackingcoverage in the feature information may want to be ordered
                write_keys.extend(val.get_write_order())

            # to preserve order of other keys - iterate instead of using set logic
            for key in val._attributes:
                if key not in write_keys:
                    write_keys.append(key)
            # write_keys.extend(set(self._attributes.keys()).difference(write_keys))
            write_array = []
            for val in self:
                list_vals = []
                for key in write_keys:
                    try:
                        v = val._attributes[key]
                    except KeyError as key_err:
                        raise KeyError(
                            "{} in {} is missing data, this would give a mismatched array \n  please fill all data {} for all items in the list/dataset".format(
                                key_err.args[0], self.metadata_name, str(write_keys)))

                    if isinstance(v, bytes):
                        # no longer doing this--
                        # convert unicode strings into ascii since HDF5 doesn't like the unicode strings that numpy will produce
                        # v = v.encode("utf-8")

                        # convert bytes strings to strings which will be encoded as utf8 later
                        v = v.decode()
                    elif isinstance(v, Enum):  # convert Enums to integars
                        v = v.value
                    list_vals.append(v)

                # list_vals = [v if not isinstance(v, str) else v.encode("utf-8") for v in list_vals]
                write_array.append(list_vals)

            # hdf5 needs names to the columns which is done in a record array or structured array.
            # but to create that without specifying type we need to transpose first then call 'fromarrays'
            transposed_array = list(map(list, zip(*write_array)))
            if write_keys:
                rec_array_revised = convert_numpy_strings_to_h5py(transposed_array, write_keys)
            else:
                rec_array_revised = h5py.Empty("")
                raise ValueError(self.metadata_name + " had no data fields defined to write - this would create an h5py.Empty dataset")
            try:
                del group_object[self.metadata_name]
            except KeyError:
                pass  # didn't exist, no error
            dataset = group_object.create_dataset(self.metadata_name, data=rec_array_revised)
            self.write_simple_attributes(dataset)
        return dataset


class S1xxGridsBase(S1xxWritesOwnGroupBase):
    @property
    @abstractmethod
    def metadata_name(self) -> str:
        raise NotImplementedError()

    def read(self, group_object_parent):
        group_object = group_object_parent[self.metadata_name]
        logging.debug("reading grids")
        self.read_simple_attributes(group_object)
        # for attr in self.get_standard_properties():
        #    setattr(self, attr, group_object[getattr(self, "__" + attr + self._attr_name_suffix)])
        for name in group_object.dtype.names:
            self._attributes[name] = group_object[name]

    def write(self, group_object):
        # @todo - is there a bug here if some instances are missing attributes leading to a mismatched array?
        """ Write out the dataset using order specified with any extra values as unordered but named at the end.

        Parameters
        ----------
        group_object
            HDF5 object to write into

        Returns
        -------
        HDF5 dataset created during the write method
        """

        # First determine the write order of the keys
        logging.debug("Writing" + " " + str(self))

        dataset = None

        write_keys = []
        if self.get_write_order():  # @todo I think bathycoverage and trackingcoverage in the feature information may want to be ordered
            write_keys.extend(self.get_write_order())

        # to preserve order of other keys - iterate instead of using set logic
        for key, val in self._attributes.items():
            if key not in write_keys and isinstance(val, s1xx_sequence_types):
                write_keys.append(key)
        # write_keys.extend(set(self._attributes.keys()).difference(write_keys))
        write_array = [self._attributes[key] for key in write_keys]

        write_compound_dtype = []
        if self.get_compound_dtype():
            write_compound_dtype.extend(self.get_compound_dtype())
        if len(write_keys) != len(write_compound_dtype):
            raise Exception("write keys and write_compound_dtype must be same length {} vs {}".format(write_keys, write_compound_dtype))
        # hdf5 needs names to the columns which is done in a record array or structured array.
        # but to create that without specifying type we need to transpose first then call 'fromarrays'

        # numpy.array is coming out with wrong (at least different) shape and fromarrays is working -- not sure why right now.
        # rec_array = numpy.array(write_array, dtype=[(name, 'f4') for name in write_keys])
        rec_array = numpy.core.records.fromarrays(write_array, dtype=[(name, dtype) for name, dtype in zip(write_keys, write_compound_dtype)])
        dataset = group_object.create_dataset(self.metadata_name, data=rec_array, chunks=True, compression='gzip', compression_opts=9)
        #         # noinspection PyAttributeOutsideInit
        # pylint: disable=attribute-defined-outside-init
        self.write_simple_attributes(dataset)


class S1XXFile(h5py.File):
    """
    hdf5 files have primary creation methods of
    create_dataset  to insert array data
    attrs           a dictionary-like to add/read metadata about the current group
    create_group    to make a group containing datasets and/or metadata
    """

    def __init__(self, *args, **kywrds):
        # @TODO: This is the NAVO default setting, have to decide if that is best and handle other options too.
        kywrds.setdefault('root', None)
        self.root = None
        self.__root_type__ = kywrds.pop('root')
        if "driver" in kywrds:
            if kywrds['driver'] == 'family':  # @todo @fixme -- this is from the NAVO files, figure how to set memb_size automatically.
                kywrds.setdefault('memb_size', 681574400)
        super().__init__(*args, **kywrds)
        # initialize with the s102 data if the file already exists.
        # if this is an empty file or opening for write then this is essentially a no-op
        if self.__root_type__:
            self.read()

    def read(self):
        self.root = self.__root_type__()
        self.root._hdf5_path = "/"
        self.root.read(self)

    def write(self):
        self.root._hdf5_path = "/"
        self.root.write(self)

    def create_empty_metadata(self):
        self.root = self.__root_type__(True)

    def show_keys(self, obj, indent=0):
        try:  # print attributes of dataset or group
            print("    " * indent + "ATTRS: " + str(list(obj.attrs.items())))
        except:
            print("    " * indent + "No attributes")
        if hasattr(obj, "keys"):
            for k in obj.keys():
                print("    " * indent + k)
                self.show_keys(obj[k], indent + 1)
        else:
            print("    " * indent + str(obj))
            indent = indent + 1
            try:  # print out any dataset arrays
                for n in obj.dtype.names:
                    try:
                        s = str(obj[n][:10])
                        s = "    " * (indent + 1) + s.replace("\n", "\n" + "    " * (indent + 1))
                        print("    " * indent, n, obj[n].shape)
                        print(s)
                    except:
                        traceback.print_exc()
            except:
                try:
                    s = str(obj[:])
                    s = "    " * (indent + 1) + s.replace("\n", "\n" + "    " * (indent + 1))
                    print("    " * indent, obj.shape)
                    print(s)
                except:
                    print("    " * indent + "dtype not understood")
