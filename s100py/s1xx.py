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


class S1XX_Attributes_base(ABC):
    """ This class implements a general hdf5 group object that has attributes, dataset or sub-groups.
    Works with S1XX_MetadataList_base if the subgroups have multiple occurences (like Group.01, Group.02)
    Works with S1XX_Dataset_base for things that are stored like a numpy array (dataset) in hdf5

    __version__ and S1XX_version must be overridden.
    To call a base class property use super().property, e.g. super().__version__
    This base class is built from the version 2.0.0 that was eventually published Nov. 2019
    """
    _attr_name_suffix = "_attribute_name"

    def __init__(self, fill_empty=False, **kywrds):
        self._attributes = collections.OrderedDict()
        if fill_empty:
            self.initialize_properties(fill_empty)
        for ky, val in kywrds.items():
            if ky in self.get_standard_properties():
                exec("self.{} = val".format(ky))
            else:
                self.add_metadata(ky, val)
        # self._child_groups = {}

    @property
    @abstractmethod
    def __version__(self) -> int:
        return -1

    @property
    def S1XX_version(self) -> tuple:
        return (2, 0, 0)

    def read_hdf5_attributes(self, group_object, indent=0):
        indentstr = "    " * indent
        logging.debug(indentstr + "Reading attributes" + str(self))
        expected_items = self.get_standard_properties_mapping()
        # basic attributes -- should be simple types so just set them
        for attr_name in group_object.attrs:
            if attr_name not in expected_items:
                logging.info(indentstr + " The attr/val: " + attr_name + "/" + str(
                    group_object.attrs[attr_name]) + " was in the group_object but not found in the standard attributes")
                self._attributes[attr_name] = group_object.attrs[attr_name]
            else:
                use_type = self.__getattribute__(expected_items[attr_name] + "_type")
                if issubclass(use_type, Enum):
                    logging.debug(indentstr + " Enumerated attr/val: " + attr_name + "/" + str(group_object.attrs[attr_name]) + " found and read")
                    self.set_enum_attribute(group_object.attrs[attr_name], attr_name, use_type)
                else:
                    logging.debug(indentstr + " Standard attr/val: " + attr_name + "/" + str(group_object.attrs[attr_name]) + " found and read")
                    setattr(self, expected_items[attr_name], group_object.attrs[attr_name])

    def __repr__(self):
        s = str(self.__class__) + "\n"
        s += str(self._attributes)
        return s

    # def _remove_attr(self, attr_name):
    #     try:
    #         self._attributes.pop(attr_name)
    #     except KeyError:
    #         pass

    def __delattr__(self, item):
        # mapping = self.get_standard_properties_mapping()
        if item in self._attributes or item in self.get_standard_properties_mapping():
            del self._attributes[item]
        elif item in self.get_standard_properties():
            del self._attributes[eval("self.{}_attribute_name".format(item))]
        else:
            del self.__dict__[item]

    def read(self, group_object, indent=0):
        """
        Parameters
        ----------
        group_object
            The group (an h5py.File is a group too) to read from.

        Returns
        -------

        """

        indentstr = "    " * indent
        logging.debug(indentstr + "Reading " + str(self))
        self.read_hdf5_attributes(group_object, indent)

        expected_items = self.get_standard_properties_mapping()

        # keys are HDF5 groups or datasets
        group_lists = self.get_standard_list_properties()
        try:
            all_keys = list(group_object.keys())
        except:
            pass
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
                    indentstr + data_key + " is an HDF5 group and was in the group_object but not found in the standard attributes, SKIPPING!")

        # now read the basic keys and the keys having list data
        # -- only need to call the read once for things that are lists as the Metadata_List class will find all of its occurrences
        for key in basic_keys:
            # create/clear the data
            logging.debug(indentstr + " Standard attr/val: " + key + " found and reading")

            # read in the HDF5 attributes etc from the group
            use_type = self.__getattribute__(expected_items[key] + "_type")
            if issubclass(use_type, S1XX_Attributes_base):
                self.__getattribute__(expected_items[key] + "_create")()
                self.__getattribute__(expected_items[key]).read(group_object[key], indent=indent + 1)
            else:
                self.__setattr__(expected_items[key], group_object[key])
        for list_type_group in list_type_keys:
            # create/clear the data
            logging.debug(indentstr + " Standard LIST based attr/val: " + list_type_group + " found and reading")
            self.__getattribute__(list_type_group + "_create")()
            # read in the HDF5 attributes etc from the group
            o = self.__getattribute__(list_type_group)
            o.read(group_object, indent=indent + 1)

    def write(self, group_object, indent=0):
        try:
            # iterate through all the key/values in _attributes and write to hdf5
            # if a value is numpy.array then convert to h5py.Dataset
            # if a value is a S1XX_Attributes_base instance then create a group and call it's write function
            # if a value is a date, time - convert to character string per S100, section 10C-7 table 10C-1
            # otherwise write as a simple attribute and simple type
            logging.debug(indent * "  " + "Writing" + " " + str(self))
            for key, val in self._attributes.items():
                key_repr = indent * "  " + key
                if isinstance(val, s1xx_sequence.__args__):  # this looks inside the typing.Union to see what arrays should be treated like this
                    logging.debug(key_repr + " array: " + str(val.shape))
                    # convert any strings to bytes of h5py will fail to write the unicode
                    converted_vals = [v if not isinstance(v, str) else v.encode("utf-8") for v in val]
                    # transposed_array = list(map(list, zip(*write_array)))
                    # rec_array = numpy.core.records.fromarrays(transposed_array, names=write_keys)
                    new_dataset = group_object.create_dataset(key, data=converted_vals)
                elif isinstance(val, S1XX_WritesOwnGroup_base):
                    # things that either create a dataset and have to combine data into it or make multiple sub groups that the parent can't predict
                    logging.debug("{}  S100 object - writing itself now...".format(key_repr))
                    val.write(group_object, indent=indent + 1)
                elif isinstance(val, S1XX_Attributes_base):
                    logging.debug(key_repr + " S100 object - writing itself now...")
                    new_group = group_object.create_group(key)
                    val.write(new_group, indent=indent + 1)
                elif isinstance(val, (datetime.date, datetime.datetime, datetime.time)):
                    logging.debug(key_repr + " datetime: {}", val)
                    # @TODO: figure out how to write datetimes
                    raise NotImplementedError("DateTimes not supported yet")
                elif isinstance(val, Enum):
                    logging.debug(key_repr + " enumeration: " + str(val))
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
                    logging.debug(key_repr + " simple type: " + str(val))
                    group_object.attrs[key] = val
            # raise NotImplementedError()
        except Exception as e:
            raise e

    def write_as_xml(self, etree_object):
        raise NotImplementedError("flesh this out if we want an xml representation of S100+ file")
        # basically add a flag to read/write functions, then everywhere a group, dataset or attribute is written either use xml or hdf5

    def get_metadata(self, key):
        return self._attributes[key]

    def add_metadata(self, key, value):
        self._attributes[key] = value

    def get_s1xx_attr(self, s1xx_name):
        expected_items = self.get_standard_properties_mapping()
        return self.__getattribute__(expected_items[s1xx_name])

    def set_s1xx_attr(self, s1xx_name, val):
        # could make this expected_items a class variable
        # (- forget the function name __class__ or __new__?) or instance variable (compute in __init__)
        expected_items = self.get_standard_properties_mapping()
        setattr(self, expected_items[s1xx_name], val)

    def get_write_order(self):
        """ Override this method if the write order of attributes/groups/dataset items is important
        Returns
        -------
        A list of key names if order is important, None otherwise.
        """
        return None

    def get_all_keys(self):  # this includes custom keys
        current_keys = set(self._attributes.keys())
        total_keys = current_keys.update(self.get_standard_keys())
        return total_keys

    # def get_standard_keys(self):  # this is only things that have properties associated
    #     props = [p[0] for p in inspect.getmembers(self.__class__, lambda x: isinstance(x, property))]
    #     implemented_properties = [self.__getattribute__(p) for p in props if p.endswith(self._attr_name_suffix) and p[:-len(self._attr_name_suffix)] in props]
    #     return implemented_properties

    def get_standard_keys(self):  # this is only things that have properties associated
        """
        Returns
        -------
        list
            The S102 HDF5 group/attribute/dataset names from this object that will be written or read from an HDF5 file.
            e.g. BathymetryCoverage or westBoundLongitude etc.

            For the class "Root":
            ['BathymetryCoverage', 'Group_F', 'TrackingListCoverage']
        """
        return list(self.get_standard_properties_mapping().keys())
        # props = [p[0] for p in inspect.getmembers(self.__class__, lambda x: isinstance(x, property))]
        # implemented_properties = [self.__getattribute__(p) for p in props if p.endswith(self._attr_name_suffix) and p[:-len(self._attr_name_suffix)] in props]
        # return implemented_properties

    def get_standard_list_properties(self):
        """

        Returns
        -------

        """
        s100_to_property = self.get_standard_properties_mapping()
        s100_to_property_for_lists = {}
        for s100_attr, prop in s100_to_property.items():
            if issubclass(self.__getattribute__(prop + "_type"), S1XX_MetadataList_base):
                s100_to_property_for_lists[s100_attr] = prop
        return s100_to_property_for_lists

    def get_standard_properties_mapping(self):
        """ This function autodetermines the HDF5 or xml names and their associated property names.

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
            s100_to_property[self.__getattribute__(prop + self._attr_name_suffix)] = prop
        return s100_to_property

    def initialize_properties(self, fill_empty=False):
        """ Calls the create function for all the properties of the class.

        Returns
        -------
        None
        """
        for prop in self.get_standard_properties():
            exec("self.{}_create()".format(prop))
            o = eval("self.{}".format(prop))
            if fill_empty and isinstance(o, S1XX_Attributes_base):
                o.initialize_properties(fill_empty)

    def get_standard_properties(self):
        """  This function autodetermines the properties implemented (which have get/set and _attribute_name methods associated)

        Returns
        -------
        list
            Names of the properties implemented.

            For class "Root":
            ['bathymetry_coverage', 'feature_information', 'tracking_list_coverage']
        """
        props = [p[0] for p in inspect.getmembers(self.__class__, lambda x: isinstance(x, property))]
        implemented_properties = [p[:-len(self._attr_name_suffix)] for p in props if
                                  p.endswith(self._attr_name_suffix) and p[:-len(self._attr_name_suffix)] in props]
        return implemented_properties

    def set_enum_attribute(self, val, attribute_name, enum_type):
        if isinstance(val, str):
            val = enum_type[val]
        if isinstance(val , (int, numpy.integer)):
            val = enum_type(val)
        self._attributes[attribute_name] = val

class S1XX_WritesOwnGroup_base(S1XX_Attributes_base):
    """ Derive things that either create a dataset and have to combine data into it or make multiple sub groups that the parent can't predict
    The S1XX_Attributes_base will call the derived class' writer without pre-making group for it.
    i.e. the derived class can specify its own group name or dataset and apply specialized logic as needed.
    """
    pass


class S1XX_MetadataList_base(list, S1XX_WritesOwnGroup_base):
    """ This class represents arrays (noted in UML as *, 1..*, 0..* etc) which is not really part of HDF5.
    The S100 spec is using a atttribute.NNN to repreent this type of record.
    This class takes the supplied name and type and will make it act like a list in python and read/write the data in HDF5 like S102 wants.

    """
    read_re_pattern = r"[_\.](\d+)"
    write_format_str = "_%03d"

    def __init__(self, *args, **opts):
        # initialize the list in case data was passed in.
        super().__init__(*args, **opts)  # standard init for lists
        S1XX_Attributes_base.__init__(self)  # initialize the s100 class

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

    def read(self, group_object, indent=0):
        indentstr = "    " * indent
        logging.debug(indentstr + "Reading " + str(self))

        # keys are HDF5 groups or datasets
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
            obj.read(group_object[data_key], indent=indent + 1)
            self.append(obj)

    def write(self, group_object, indent=0):
        # Iterate through the values in the list and write to hdf5
        # Write each item as self.metadata_name + ".%03d" % index
        # They should all be the same type but we aren't sure what type they are.
        # Most likely to be another S100_Attribute type.
        # if a value is numpy.array then convert to h5py.Dataset
        # if a value is a S1XX_Attributes_base instance then create a group and call it's write function
        # if a value is a date, time - convert to character string per S100, section 10C-7 table 10C-1
        # otherwise write as a simple attribute and simple type
        try:
            logging.debug(indent * "  " + "Writing" + " " + str(self))
            # create N new group objects named as metadata_name.NNN
            for index, val in enumerate(self):
                name = self.metadata_name + ".%03d" % (index + 1)
                if isinstance(val, s1xx_sequence.__args__):  # this looks inside the typing.Union to see what arrays should be treated like this
                    raise NotImplementedError()
                # elif isinstance(val, S1XX_MetadataList_base):
                #     raise NotImplementedError("Nested Lists")
                elif isinstance(val, S1XX_Attributes_base):  # Attributes, List and Datasets all work the same (for now)
                    new_group = group_object.create_group(name)
                    val.write(new_group, indent + 1)
                elif isinstance(val, (datetime.date, datetime.datetime, datetime.time)):
                    logging.debug(name + " datetime: {}", val)
                    # @TODO: figure out how to write datetimes
                    raise NotImplementedError("DateTimes not supported yet")
                else:
                    logging.debug(name + " simple type: " + str(val))
                    group_object[name] = val
        except Exception as e:
            raise e

    def write_as_xml(self, etree_object):
        raise NotImplementedError("flesh this out if we want an xml representation of S102 bathy file")


class S1XX_Dataset_base(list, S1XX_WritesOwnGroup_base):
    """ The S102 spec stores some things as attributes that could (or should) be stored as attributes.
    This class reads/writes datasets but stores/accesses them as a list of class instances.
    Data access should then be used as object[index].attribute
    So for the FeatureInformation class that would be feat[0].name = "depth" and feat[1].name = "uncertainty"
    """

    def __init__(self, *args, **opts):
        # initialize the list in case data was passed in.
        super().__init__(*args, **opts)  # standard init for lists
        S1XX_Attributes_base.__init__(self)  # initialize the s102 class

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
        s = S1XX_Attributes_base.__repr__(self)
        for data_object in self:
            s += str(data_object)
        return s

    def read(self, group_object, indent=0):
        self.read_hdf5_attributes(group_object, indent)  # put any attributes from the dataset obect into the overall _attributes

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

    def write(self, group_object, indent=0):
        # @todo - is there a bug here if some instances are missing attributes leading to a mismatched array?
        """ Write out the dataset using order specified with any extra values as unordered but named at the end.

        Parameters
        ----------
        group_object
            HDF5 object to write into
        indent

        Returns
        -------
        HDF5 dataset created during the write method

        """
        try:
            # First determine the write order of the keys
            logging.debug(indent * "  " + "Writing" + " " + str(self))

            dataset = None
            if len(self) > 0:
                val = self[0]
                write_keys = []
                if self.get_write_order():  # @todo I think bathycoverage and trackingcoverage in the feature information may want to be ordered
                    write_keys.extend(val.get_write_order())

                # to preserve order of other keys - iterate instead of using set logic
                for key in val._attributes:
                    if key not in write_keys:
                        write_keys.append(key)
                # write_keys.extend(set(self._attributes.keys()).difference(write_keys))
                write_array = []
                for val in self:
                    list_vals=[]
                    for key in write_keys:
                        v = val._attributes[key]
                        if isinstance(v, str):  # convert unicode strings into ascii since HDF5 doesn't like the unicode strings that numpy will produce
                            v = v.encode("utf-8")
                        elif isinstance(v, Enum):  # convert Enums to integars
                            v = v.value
                        list_vals.append(v)

                    # list_vals = [v if not isinstance(v, str) else v.encode("utf-8") for v in list_vals]
                    write_array.append(list_vals)

                # hdf5 needs names to the columns which is done in a record array or structured array.
                # but to create that without specifying type we need to transpose first then call 'fromarrays'
                transposed_array = list(map(list, zip(*write_array)))
                rec_array = numpy.core.records.fromarrays(transposed_array, names=write_keys)
                dataset = group_object.create_dataset(self.metadata_name, data=rec_array)
            return dataset
        except Exception as e:
            raise e

    def write_as_xml(self, etree_object):
        raise NotImplementedError("flesh this out if we want an xml representation of S102 bathy file")



class S1XXFile(h5py.File):
    """
    hdf5 files have primary creation methods of
    create_dataset  to insert array data
    attrs           a dictionary-like to add/read metadata about the current group
    create_group    to make a group containing datasets and/or metadata
    """
    # these keys allow backward compatibility, the first key is current at time of writing
    top_level_keys = ('BathymetryCoverage', 'S102_Grid', 'S102_BathymetryCoverage')
    tracking_list_top_level = ("TrackingListCoverage",)
    tracking_list_second_level = ("TrackingListCoverage.01",)
    tracking_list_group_level = ("Group.001",)
    second_level_keys = (
        'BathymetryCoverage.01', 'S102_Grid.01', 'S102_BathymetryCoverage.01', 'BathymetryCoverage_01', 'S102_Grid_01', 'S102_BathymetryCoverage_01',)
    group_level_keys = ('Group.001', 'Group_001',)
    value_level_keys = ("values",)
    depth_keys = ("depth", "depths", 'elevation', "elevations", "S102_Elevation")

    def __init__(self, *args, **kywrds):
        # @TODO: This is the NAVO default setting, have to decide if that is best and handle other options too.
        kywrds.setdefault('root', None)
        self.root_type = kywrds.pop('root')
        if "driver" in kywrds:
            if kywrds['driver'] == 'family':  # @todo @fixme -- this is from the NAVO files, figure how to set memb_size automatically.
                kywrds.setdefault('memb_size', 681574400)
        super().__init__(*args, **kywrds)
        # initialize with the s102 data if the file already exists.
        # if this is an empty file or opening for write then this is essentially a no-op
        if self.root_type:
            self.read()

    def read(self):
        self.root = self.root_type()
        self.root.read(self)

    def write(self):
        self.root.write(self)

    def create_empty_metadata(self):
        self.root = self.root_type(True)

    def print_overview(self, display_nodes=10):
        depths = self.get_depths()
        print("shape of grid is", depths.shape, "of type", depths.dtype)
        with numpy.printoptions(precision=2, suppress=True, linewidth=200):
            x, y = depths.shape
            r = max(x, y)
            step = int(r / display_nodes)
            print(depths[::step, ::step])

    def print_depth_attributes(self):
        hdf5 = self.get_depth_dataset()
        print(hdf5.attrs)

    def get_depth_dataset(self):
        for k in self.top_level_keys:
            if k in self:
                d = self[k]
                break
        try:
            d
        except NameError:
            raise KeyError(str(self.top_level_keys) + " were not found in " + str(list(self.keys())))

        for k in self.second_level_keys:
            if k in d:
                g = d[k]
                break

        try:
            g
        except NameError:
            raise KeyError(str(self.second_level_keys) + " were not found in " + str(list(d.keys())))

        for k in self.group_level_keys:
            if k in g:
                gp = g[k]
                break

        try:
            gp
        except NameError:
            raise KeyError(str(self.group_level_keys) + " were not found in " + str(list(g.keys())))

        for k in self.value_level_keys:
            if k in gp:
                v = gp[k]
                break
        try:
            v
        except NameError:
            raise KeyError(str(self.value_level_keys) + " were not found in " + str(list(gp.keys())))
        return v

    def get_depths(self):
        v = self.get_depth_dataset()
        # v.dtype
        # dtype([('S102_Elevation', '<f4'), ('S102_Uncertainty', '<f4')])
        for k in self.depth_keys:
            if k in v.dtype.names:
                return v[k]
        raise KeyError(str(self.depth_keys) + " were not found in " + str(list(v.dtype.names)))

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
