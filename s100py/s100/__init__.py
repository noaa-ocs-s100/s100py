"""
S-100 Python Utilities

Tools for converting various datasets to S-100 compliant formats.
"""

import pathlib

from .v5_2.api import *
import s100py

def open(filename: (str, pathlib.Path), mode: str = "r") -> S1XXFile:
    # sends back the S100 version if it can't figure out which exact type to use
    file_object = S100File(filename, mode)
    # @todo use the new case statement when python 3.10 is required
    spec = file_object.root.product_specification

    for product in [s100py.s102.v3_0, s100py.s102.v2_2, s100py.s102.v2_1, s100py.s102.v2_0,
                    s100py.s111.v1_0, s100py.s111.v1_2, s100py.s111.v2_0,
                    s100py.s104.v1_0, s100py.s104.v1_1, s100py.s104.v2_0]:
        if h5py_string_comp(spec, product.PRODUCT_SPECIFICATION):
            if 'S-102' in spec:
                file_object = product.api.S102File(filename, mode)
            elif 'S-111' in spec:
                file_object = product.api.S111File(filename, mode)
            elif 'S-104' in spec:
                file_object = product.api.S104File(filename, mode)
            else:
                print("Warning: unrecognized s100 product specification, using generic S100File object")
            break

    else:
        print("Warning: unrecognized s100 product specification, using generic S100File object")
    return file_object

